### Processing of Exons and their ranges

#### The script "exon_processing.py" lists out all the strong and exons for the splice variants genes:
	DATA=../../data
	OUT=$DATA/output_data/exon_processed_files
	PYTHONPATH=../common python exon_processing.py \
	$DATA/strandomics_input_data/transcripts.tsv \
	$DATA/variants_input_data/splice_variants.tsv \
	$DATA/strandomics_input_data/genes.tsv \
	$OUT/strong_exons.tsv \
	$OUT/rare_exons.tsv

#### The script "get_sequences_for_exon_distribution.py" consumes exons file and gets the sequences for them:
	DATA=../../data/output_data/
	PYTHONPATH=../common:$PYTHONPATH
	PYTHONPATH=../create_splice_matrices:$PYTHONPATH
	export PYTHONPATH
	python get_sequences_for_exon_distribution.py \
	$DATA/exon_processed_files/strong_exons.tsv \
	$DATA/exon_with_sequences/strong_exons_sequences.tsv

	DATA=../../data
	PYTHONPATH=../common:$PYTHONPATH
	PYTHONPATH=../create_splice_matrices:$PYTHONPATH
	export PYHTONPATH
	python get_sequences_for_exon_distribution.py \
	$DATA/exon_processed_files/rare_exons.tsv \
	$DATA/exon_with_sequences/rare_exons_sequences.tsv
