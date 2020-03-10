echo "Starting..."
DATA=../../data
PYTHONPATH=../common python ../get_reported_variants_statistics/get_reported_variants_statistics.py \
$DATA/variants_input_data/positive.tsv \
$DATA/variants_input_data/negative.tsv \
$DATA/variants_input_data/cohort_tests.txt \
$DATA/variants_input_data/cohort_gene_lists.txt