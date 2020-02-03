### The script "get_exon_distribution_data.py" lists out all the strong exons for the splice variants genes:
	DATA=../../data
	PYTHONPATH=../common python get_exon_distribution_data.py \
	$DATA/variants_input_data/splice_variants.tsv \
	$DATA/strandomics_input_data/transcripts.tsv
