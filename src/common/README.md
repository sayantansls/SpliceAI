### This directory contains the scripts which are used commonly or post generation processing

#### The script 'utils.py' has various functions which are used in other scripts


#### The script 'multiseq_processing.py' consumes a raw data TSV file containing sequences for a particular splice predictor and converts it into a FASTA format which can be used for bulk upload on the respective splice predictor sites

	DATA=../../data
	python multiseq_processing.py \
	$DATA/processed_files/splice_predictor_raw_data/Splice_Variants_Matrix\ -\ MaxEntScan\ Raw\ Data.tsv \
	$DATA/processed_files/fasta_upload_files/maxentscan_ref.fasta \
	$DATA/processed_files/fasta_upload_files/maxentscan_var.fasta

	DATA=../../data
	python multiseq_processing.py \
	$DATA/processed_files/splice_predictor_raw_data/Splice_Variants_Matrix\ -\ NNSplice\ Raw\ Data.tsv \
	$DATA/processed_files/fasta_upload_files/maxentscan_ref.fasta \
	$DATA/processed_files/fasta_upload_files/maxentscan_var.fasta

	DATA=../../data
	python multiseq_processing.py \
	$DATA/processed_files/splice_predictor_raw_data/Splice_Variants_Matrix\ -\ ASSP\ Raw\ Data.tsv \
	$DATA/processed_files/fasta_upload_files/maxentscan_ref.fasta \
	$DATA/processed_files/fasta_upload_files/maxentscan_var.fasta

