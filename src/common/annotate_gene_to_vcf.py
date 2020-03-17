"""
@author : sayantan
This file consumes a vcf file, maps and annotates the gene name to another column
based on the position of the variant
"""
import time as tm
import utils
import sys

"""
The headers in the genes.tsv file are as follows:
1 -- Strand_gene_id
2 -- ChrName
3 -- Strand
4 -- Gene_start
5 -- Gene_end
6 -- Symbol
7 -- Entrez_id
"""
def create_gene_position_map(genesfile):
	global gene_positions_map
	gene_positions_map = dict()
	for gene in utils.records_iterator(genesfile):
		gene_positions_map[gene['Symbol']] = {'start': gene['Gene_start'], 'end': gene['Gene_end']}

def annotate_entries(vcffile):
	pass

def main(vcffile, genesfile):
	print("Start of code:", tm.ctime(tm.time()))

	create_gene_position_map(genesfile)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])