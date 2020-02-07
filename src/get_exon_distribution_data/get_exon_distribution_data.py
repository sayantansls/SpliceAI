# author -- sayantan
# The purpose of this code is as follows:
# 1 -- This code parses through and obtains the genes in the splice_variants.tsv file
# 2 -- Then for those genes it parses through the file transcripts.tsv and obatins the strong exons (present in all transcripts) 
#      for those corresponding genes

from __future__ import print_function
import re
import utils
import time as tm 

# The headers in splice_variants.tsv file are as follows:
# 1 -- Source	
# 2 -- Case ID	
# 3 -- UK1	
# 4 -- UK2	
# 5 -- Panel	
# 6 -- Test	
# 7 -- Time and Date	
# 8 -- ClinicalManifestation	
# 9 -- Gender	
# 10 -- Ethnicity	
# 11 -- Age	
# 12 -- ReportWiseGeneListInfo
# 13 -- UK3	
# 14 -- UK4	
# 15 -- UK5	
# 16 -- Chromosome	
# 17 -- Gene	
# 18 -- pHGVS	
# 19 -- Genomic HGVS	
# 20 -- Transcript	
# 21 -- VariantType	
# 22 -- AlleleFrequency	
# 23 -- Global PPDB	
# 24 -- Local PPDB	
# 25 -- Zygosity	
# 26 -- Clinvar IDs	
# 27 -- Variant Label Reason	
# 28 -- UK6	
# 29 -- rsID	
# 30 -- EVS	
# 31 -- ExAc	
# 32 -- dbSNP	
# 33 -- 1000 Genomes	
# 34 -- HGMD ID	
# 35 -- Bioinfo Summary	
# 36 -- Literature Summary

# The headers in the transcripts.tsv file are as follows:
# 1 -- Strand_gene_id	
# 2 -- Transcript_name	
# 3 -- Transcript_start	
# 4 -- Transcript_end	
# 5 -- CDS_start	
# 6 -- CDS_end	
# 7 -- Exon_count	
# 8 -- Exon_start	
# 9 -- Exon_end	
# 10 -- Str_validity	
# 11 -- Genomic_override	
# 12 -- Non_standard_codon	
# 13 -- Version	
# 14 -- is_deleted	
# 15 -- invalid_reason	
# 16 -- positional_override	
# 17 -- Comments

splice_variant_genes = set()

def create_splice_variants_genes_set(splice_variants_file):
	for variant in utils.records_iterator(splice_variants_file):
		splice_variant_genes.add(variant['Gene'])

gene_exons = dict()

def get_strong_exons(transcripts_file):
	for entry in utils.records_iterator(transcripts_file):
		if '-'.join(entry['Strand_gene_id'].split('-')[:-2]) in splice_variant_genes:
			if entry['Strand_gene_id'] not in gene_exons:
				gene_exons[entry['Strand_gene_id']] = {'exon_start': list(),
														'exon_end': list()}
			gene_exons[entry['Strand_gene_id']]['exon_start'].append(entry['Exon_start'])
			gene_exons[entry['Strand_gene_id']]['exon_end'].append(entry['Exon_end'])

def main(splice_variants_file, transcripts_file):
	print("Start of code:", tm.ctime(tm.time()))

	create_splice_variants_genes_set(splice_variants_file)
	get_strong_exons(transcripts_file)
	print(gene_exons)
	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2])