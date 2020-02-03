## Getting the statistical analysis of the reported variants files - positive.tsv and negative.tsv

### The format of the statistics is as follows:

####P ositive File statistics:
 #Cases in the file: 5979
 #Cases in the clinical exome panel: 2758
 #Cases in the clinical exome panel (LP + P + VUSD): 1637
 #Cases in the clinical exome + neurodevelopmental cohort: 307
 #Cases in the clinical exome + neurodevelopmental cohort (LP + P + VUSD): 201
 Gender Stats: {'Male': 117, 'Female': 84}
 Age Stats: {'below 10 years': 116, 'above 20 years': 41, 'between 10 to 20 years': 44}
 Variant Label Stats: {'Pathogenic': 123, 'VUSD': 39, 'VUS': 0, 'Likely Benign': 0, 'Likely Pathogenic': 39}

#### Negative File statistics:
 #Cases in the file: 6237
 #Cases in the clinical exome panel: 913
 #Cases in the clinical exome + neurodevelopmental cohort: 125
 Gender Stats: {'Male': 74, 'Female': 51}
 Age Stats: {'below 10 years': 72, 'above 20 years': 18, 'between 10 to 20 years': 35}

### The script "get_reported_variants_statistics.py" gives us the reported variants statitics:
	DATA=../../data
	PYTHONPATH=../common python get_reported_variants_statistics.py \
	$DATA/variants_input_data/positive.tsv \
	$DATA/variants_input_data/negative.tsv \
	$DATA/variants_input_data/cohort_tests.txt \
	$DATA/variants_input_data/cohort_gene_lists.txt
