import math
import numpy as np
from tqdm import tqdm
from os.path import join
from os import listdir
import pandas as pd
import sys
import subprocess

def get_read_data(all_data):
    '''
    Given a line from an NGS file, return the name, sequence, and quality scores of the read
    Args:
        all_data: string containing name, read, +, quality
    Returns:
        name: identifier for the read
        read: DNA sequences read
        quality: quality scores corresponding to each position in the read
    '''

    all_data = all_data.split('\n')
    name = all_data[0] 
    try:
        read = all_data[1]
    except:
        print(all_data)
    quality = all_data[3]
    return name, read, quality


def process_read_pair(fwd_all_data, rev_all_data):
    '''
    Take a pair of reads (including their metadata), align them and return a compiled read
    Args:
        fwd_all_data: raw string contain name, read, +, quality for the forward read
        rev_all_data: raw string contain name, read, +, quality for the reverse read
    Returns:
        A processed read: aligned reads, select nucleotide with highest quality score for each position
    '''

    fwd_name, fwd_read, fwd_quality = get_read_data(fwd_all_data)
    rev_name, rev_read, rev_quality = get_read_data(rev_all_data)
    dna_seqs_final = []
    complement = {'C':'G', 'G':'C', 'A':'T', 'T':'A', 'N':'N'}
    if len(fwd_read) != len(fwd_quality) or len(rev_read) != len(rev_quality):
        print('read and quality score lengths not the same!')
    # take reverse complement
    fwd_read = ''.join([complement[fwd_read[len(fwd_read) - i - 1]] for i in range(len(fwd_read))])
    # convert quality scores to int
    rev_quality = [ord(q) - 33 for q in rev_quality]
    fwd_quality = [ord(fwd_quality[len(fwd_quality) - i - 1]) - 33 for i in range(len(fwd_quality))]
    # merge reads
    sequence = [None]*270
    for i in range(270):
        if i < 151:
            base = rev_read[i]
        if i >= 151 or (i >= 119 and i < 151 and rev_quality[i] < fwd_quality[i-119]):
            base = fwd_read[i-119]
        sequence[i] = base
    return ''.join(sequence)


if __name__ == '__main__':
    '''
    inputs:
        input_directory: directory that holds all the fastq files (make sure to unzip them)
        date: date the NGS was performed on, used for naming outputs
        merged_reads_directory: place to store merged read files, intermediate files
        input_seq_df: filed that should contain a column called 'dna_seq' which has the nucleotide seq we want to match
        output_counts_file: where to output the counts for each NGS file
    '''
    _, input_directory, date, merged_reads_directory, input_seq_df, output_counts_file = sys.argv
    seq_df = pd.read_csv(input_seq_df)

    # read file names
    paired_file_paths = pd.read_csv(join(input_directory, 'sra_file_pairs.csv'))

    # remove date variable from dataframe columns - makes count file more interchangeable
    merged_reads_output_names = ['both_reads'.join(file.split('R1')) for file in paired_file_paths.R1]
    read_descriptors = [file.split('R1')[0] for file in paired_file_paths.R1]

    # save output stats here to print report when done
    output_stats = []

    # explicitly identify r1 and r2 files while keeping with `file_name_pair` naming scheme below
    for r1_file, r2_file, merged_reads_output_name in tqdm(zip(paired_file_paths.R1, paired_file_paths.R2, merged_reads_output_names), total=len(paired_file_paths), ncols=100, leave=True, desc='File'):
        file_name_pair = [r1_file, r2_file]
        read_count = 0
        fwd_all_data = ''
        rev_all_data = ''

        # open read files
        f1 = open(join(input_directory, file_name_pair[0]), 'r')
        f2 = open(join(input_directory, file_name_pair[1]), 'r')
        
        # write empty merged read file
        with open(join(merged_reads_directory, merged_reads_output_name), 'w') as f:
            f.write('')

        # open in appending mode
        out_file = open(join(merged_reads_directory, merged_reads_output_name), 'a')
        for l1, l2 in tqdm(zip(f1, f2), leave=False, desc='Processing reads'):
            # if a new read description is identified and compiled fwd/rev_all_data is not '',
            # then fwd/rev_all_data is complete read. Process and identify read.
            if l1[0] == '@' and l2[0] == '@' and fwd_all_data != '' and rev_all_data != '':
                final_sequence = process_read_pair(fwd_all_data, rev_all_data)
                out_file.write(final_sequence + '\n')
                read_count += 1

                # reset fwd/rev_all_data to ''
                fwd_all_data = ''
                rev_all_data = ''
            
            # append fastq data to fwd/rev_all_data until complete read data is compiled
            # i.e. read descriptor, sequence, description, and quality scores.
            fwd_all_data += l1
            rev_all_data += l2

        # for loop does not process last read in file - do that here
        final_sequence = process_read_pair(fwd_all_data, rev_all_data)
        out_file.write(final_sequence + '\n')
        read_count += 1

        # save output stats
        output_stats.append((*file_name_pair, read_count))

        # close files
        f1.close()
        f2.close()
        out_file.close()

    # print output_stats
    print('All reads filtered for quality scores, final read counts:')
    for r1_file, r2_file, read_cnt in output_stats:
        print(f"Experiment: {r1_file.split('R1')[0]}  --  Total reads: {read_cnt}")

    print('Identifying reads in filtered fastq files...')
    unsorted_counts = None
    for merged_reads_output_name, read_descriptor in tqdm(zip(merged_reads_output_names, read_descriptors), total=len(paired_file_paths), ncols=100, leave=True, desc='File'):
        # open merged reads again
        with open(join(merged_reads_directory, merged_reads_output_name), 'r') as f:
            sequences = f.read().split('\n')
        possible_sequences = seq_df['dna_seq'].values.tolist()

        # calculate counts (add a count to a possible sequence if the read matches exactly)
        # counts index matches that of sequence index
        counts = [0 for _ in possible_sequences]
        print('Determining counts for {} reads from {} dataset'.format(len(sequences), read_descriptor))
        for i, sequence in enumerate(sequences):
            sequence = sequence[58:223]  # this is the unique protein sequence in the read
            if sequence in possible_sequences:
                counts[possible_sequences.index(sequence)] += 1
        seq_df[read_descriptor + 'count'] = counts
    
    # add date here to distinguish different counts file outputs
    seq_df.to_csv(date+'_'+output_counts_file, index=False)
