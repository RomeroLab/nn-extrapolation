DATE=$(date +'%m%d%Y')

all: preprocess
	echo "Done"

extrapolation:
	python 01_extrapolation_predictions.py

preprocess: download # check to see if this is the best way to specify to download if folders do not exist
	mkdir merged_reads
	for d in fastq_files/ ; do
	python preprocess.py fastq_files/${d} ${d:0:6} merged_reads/ designs.csv designs_counts.csv

download: 
	mkdir fastq_files
	# Download fastq files from NCBI SRA
	echo "Not implemented yet"

