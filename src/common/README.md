### This directory contains the scripts which are used commonly or post generation processing

#### The script 'utils.py' has various functions which are used in other scripts


#### The script 'multiseq_processing.py' consumes a raw data TSV file containing sequences for a particular splice predictor and converts it into a FASTA format which can be used for bulk upload on the respective splice predictor sites

	DATA=../../data/processed_files
	python multiseq_processing.py \
	$DATA/splice_predictor_raw_data/Splice_Variants_Matrix\ -\ MaxEntScan\ Raw\ Data.tsv \
	$DATA/fasta_upload_files/splice_variants/maxentscan_ref.fasta \
	$DATA/fasta_upload_files/splice_variants/maxentscan_var.fasta

	DATA=../../data/processed_files
	python multiseq_processing.py \
	$DATA/splice_predictor_raw_data/Splice_Variants_Matrix\ -\ NNSplice\ Raw\ Data.tsv \
	$DATA/fasta_upload_files/splice_variants/nnsplice_ref.fasta \
	$DATA/fasta_upload_files/splice_variants/nnsplice_var.fasta

	DATA=../../data/processed_files
	python multiseq_processing.py \
	$DATA/splice_predictor_raw_data/Splice_Variants_Matrix\ -\ ASSP\ Raw\ Data.tsv \
	$DATA/fasta_upload_files/splice_variants/assp_ref.fasta \
	$DATA/fasta_upload_files/splice_variants/assp_var.fasta

#### The script 'exon_multiseq_processing.py' consumes a raw data TSV file containing sequences of ESS sites for each of the splice predictors and converts it into a FASTA format which can be used for bulk upload on the respective splice predictor web tools

	DATA=../../data/processed_files
	python exon_multiseq_processing.py \
	$DATA/exon_distribution_data/exon_distribution\ -\ Strong_Exons_MES.tsv \
	$DATA/fasta_upload_files/exon_data/mes_strong_exon_donor.fasta \
	$DATA/fasta_upload_files/exon_data/mes_strong_exon_acceptor.fasta

	DATA=../../data/processed_files
	python exon_multiseq_processing.py \
	$DATA/exon_distribution_data/exon_distribution\ -\ Strong_Exons_NNSplice.tsv \
	$DATA/fasta_upload_files/exon_data/nnsplice_strong_exon_donor.fasta \
	$DATA/fasta_upload_files/exon_data/nnsplice_strong_exon_acceptor.fasta

	DATA=../../data/processed_files
	python exon_multiseq_processing.py \
	$DATA/exon_distribution_data/exon_distribution\ -\ Strong_Exons_ASSP.tsv \
	$DATA/fasta_upload_files/exon_data/assp_strong_exon_donor.fasta \
	$DATA/fasta_upload_files/exon_data/assp_strong_exon_acceptor.fasta

	DATA=../../data/processed_files
	python exon_multiseq_processing.py \
	$DATA/exon_distribution_data/exon_distribution\ -\ Rare_Exons_MES.tsv \
	$DATA/fasta_upload_files/exon_data/mes_rare_exon_donor.fasta \
	$DATA/fasta_upload_files/exon_data/mes_rare_exon_acceptor.fasta

	DATA=../../data/processed_files
	python exon_multiseq_processing.py \
	$DATA/exon_distribution_data/exon_distribution\ -\ Rare_Exons_NNSplice.tsv \
	$DATA/fasta_upload_files/exon_data/nnsplice_rare_exon_donor.fasta \
	$DATA/fasta_upload_files/exon_data/nnsplice_rare_exon_acceptor.fasta

	DATA=../../data/processed_files
	python exon_multiseq_processing.py \
	$DATA/exon_distribution_data/exon_distribution\ -\ Rare_Exons_ASSP.tsv \
	$DATA/fasta_upload_files/exon_data/assp_rare_exon_donor.fasta \
	$DATA/fasta_upload_files/exon_data/assp_rare_exon_acceptor.fasta

