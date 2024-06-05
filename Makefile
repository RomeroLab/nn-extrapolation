DATE=$(shell date +%m-%d-%Y)

env:
	conda env create -f environment.yml
	conda activate gb1_inf

extrapolation:
	python 01_extrapolation_predictions.py

process_sequencing:# sra_download
	echo "${DATE}"
	mkdir merged_reads
	python 03_preprocessing.py fastq_files ${DATE} merged_reads/ designs.csv designs_counts.csv

sra_download: 
	# Download fastq files from NCBI SRA
	echo "Total download size for SRA data is ~56 GB"
	cd fastq_files; \
		python download_sra_data.py

