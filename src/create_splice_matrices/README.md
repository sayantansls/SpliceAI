### Creating splice matrices for 5 different splice predictors
The splice predictors being used in StrandOmics are as follows:
a. MaxEntScan (MES)
b. SplicePort
c. NNSplice
d. Human Splicing Finder (HSF)
e. ASSP

### The script "create_splice_matrices.py" creates the desired sequences for the splice variants
	DATA=../../data
	PYTHONPATH= ../common python create_splice_matrices.py \
	$DATA/variants_input_data/splice_variants.tsv \
	$DATA/output_data/splice_matrices.tsv \
	$DATA/strandomics_clinical_data/genes.tsv

### The script "maxentscan_processing.py" creates the sequences for MaxEntScan and processes it into a FASTA file
	DATA=../../data
	OUT=../../data/processed_files/
	PYTHONPATH=../common python maxentscan_processing.py \
	$DATA/output_data/splice_matrices_maxentscan.tsv \
	$OUT/mes_donor_ref.tsv \
	$OUT/mes_donor_var.tsv \
	$OUT/mes_acceptor_ref.tsv \
	$OUT/mes_acceptor_var.tsv

### The script "multiseq_processing.py" creates the sequences for NNSplice and ASSP into FASTA format sequences
	DATA=../../data
	OUT=../../data/processed_files/
	PYTHONPATH=../common python multiseq_processing.py \
	$DATA/output_data/splice_variants_nnsplice.tsv \
	$OUT/nnsplice_ref_file.tsv \
	$OUT/nnsplice_var_file.tsv

	DATA=../../data
	OUT=../../data/processed_files/
	PYTHONPATH=../common python multiseq_processing.py \
	$DATA/output_data/splice_variants_assp.tsv \
	$OUT/assp_ref_file.tsv \
	$OUT/assp_var_file.tsv           
