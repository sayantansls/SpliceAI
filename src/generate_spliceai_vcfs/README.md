### Consume a file of variants and convert into a VCF file format (recognized by SpliceAI)

#### The format of the 'input.vcf' for Splice AI is as follows:

	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	1	25000	.	A	C,G,T	.	.	.
	2	152389953	.	T	A,C,G	.	.	.
	2	179415988	.	C	CA	.	.	.
	2	179446218	.	ATACT	A	.	.	.
	2	179446218	.	ATACT	AT,ATA	.	.	.

#### The script 'spliceai_vcf_generator.py' consumes an input TSV file of the following format:

	Gene name	genomicHGVS
	PLP1	g.103041656G>A
	MFSD8	g.128859994C>T
	CEP290	g.88512260C>T
	PNKP	g.50365794C>A 

#### The following command runs the script and creates the VCF file for the variants in the input TSV file:

	DATA=../../data
	OUT=../../data/processed_files/
	PYTHONPATH=../common python spliceai_vcf_generator.py \
	$DATA/variant_input_data/splice_variants_for_vcf.tsv \
	$OUT/splice_variants_vcf.vcf \
	$DATA/human_genome_data/hg19.2bit \
	$DATA/spliceai_templates/input_template.vcf \
	$DATA/strandomics_input_data/genes.tsv
