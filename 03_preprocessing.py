import math
import numpy as np
from os.path import join
from os import listdir
import pandas as pd
import sys

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

    file_names = listdir(input_directory)
    file_name_pairs = []
    for file_name in file_names:
        if '.fastq' in file_name and 'R1' in file_name:
            file_name_match = 'R2'.join(file_name.split('R1'))
            if file_name in file_names:
                file_name_pairs.append([file_name, file_name_match])
                print('Found pair of file reads {}, {}'.format(file_name, file_name_match))
            else:
                print('Unable to find match for {}'.format(file_name))

    merged_reads_output_names = ['both_reads'.join(file_name_pair[0].split('R1')) for file_name_pair in file_name_pairs]
    read_descriptors = [file_name_pair[0].split('R1')[0] + date for file_name_pair in file_name_pairs]

    for file_name_pair, merged_reads_output_name in zip(file_name_pairs, merged_reads_output_names):
        print('Processing {}, {}'.format(*file_name_pair))
        read_count = 0
        fwd_all_data = ''
        rev_all_data = ''
        # open each pair of files
        f1 = open(join(input_directory, file_name_pair[0]), 'r')
        f2 = open(join(input_directory, file_name_pair[1]), 'r')
        # overwrite existing merged_read_file
        with open(join(merged_reads_directory, merged_reads_output_name), 'w') as f:
            f.write('')
        # open in appending mode
        out_file = open(join(merged_reads_directory, merged_reads_output_name), 'a')
        for l1, l2 in zip(f1, f2):
            # process and write the read
            if l1[0] == '@' and l2[0] == '@' and fwd_all_data != '' and rev_all_data != '':
                final_sequence = process_read_pair(fwd_all_data, rev_all_data)
                out_file.write(final_sequence + '\n')
                read_count += 1
                fwd_all_data = ''
                rev_all_data = ''
            # append to the new line
            fwd_all_data += l1
            rev_all_data += l2
        # process and write final read
        final_sequence = process_read_pair(fwd_all_data, rev_all_data)
        out_file.write(final_sequence + '\n')
        read_count += 1
        print('Merged {} reads'.format(read_count))

    unsorted_counts = None
    for merged_reads_output_name, read_descriptor in zip(merged_reads_output_names, read_descriptors):
        # open merged reads again
        with open(join(merged_reads_directory, merged_reads_output_name), 'r') as f:
            sequences = f.read().split('\n')
        seq_ids = seq_df['seq_ID'].values.tolist()
        possible_sequences = seq_df['dna_seq'].values.tolist()
        # calculate counts (add a count to a possible sequence iff the read matches exactly)
        counts = [0 for _ in possible_sequences]
        print('Determining counts for {} reads from {} dataset'.format(len(sequences), read_descriptor))
        for i, sequence in enumerate(sequences):
            sequence = sequence[58:223] # this portion will also be unique to each protein library
            if sequence in possible_sequences:
                counts[possible_sequences.index(sequence)] += 1
        seq_df[read_descriptor + '_count'] = counts

    seq_df.to_csv(output_counts_file, index=False)
