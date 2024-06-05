"""Download SRA sequencing data for project PRJNA1117877."""

# imports
import subprocess
import shutil
import os
from tqdm import tqdm

import pandas as pd


output_dir = '.'
run_info_file = 'PRJNA1117877_run_info.csv'

run_info = pd.read_csv(run_info_file)

# download fastq files from the SRA
print('Downloading files from SRA...')
for index, row in tqdm(run_info.iterrows(), total=len(run_info), ncols=100, desc='Progress'):
    subprocess.run(
        f'fasterq-dump {row.Run} --outfile {row.LibraryName} --outdir {output_dir}',
        shell=True,
        check=True
    )

    # rename output file to use R1 and R2 to indicate read direction
    file_base = os.path.join(output_dir, row.LibraryName)

    for r in [1, 2]:
        shutil.move(f'{file_base}_{r}.fastq', f'{file_base}_R{r}.fastq')