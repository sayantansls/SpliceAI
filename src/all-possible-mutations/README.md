#### The script 'all_possible_snv_mutations.py' creates all possible substitution mutations for the essential splice sites of the strong and rare exons

	DATA=../../data
	PYTHONPATH=../common:$PYTHONPATH
	PYTHONPATH=../create_splice_matrices:$PYTHONPATH
	export PYTHONPATH
	python all_possible_snv_mutations.py \
	$DATA/output_data/exon_processed_files/strong_exons.tsv \
	$DATA/output_data/exon_processed_files/rare_exons.tsv \
	str_out.tsv rare_out.tsv
