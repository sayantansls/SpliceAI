"""
@author : sayantan
This file consumes a vcf file, maps and annotates the gene name to another column
based on the position of the variant
"""
from __future__ import print_function
import time as tm
import utils, copy
import sys, csv, os

HEADERS = ['#CHROM',
		   'POS',
		   'ID',
		   'REF',
		   'ALT',
		   'QUAL',
		   'FILTER',
		   'INFO',
		   'FORMAT',
		   'DETAILS',
		   'GENES']

ENTRY_T = {'#CHROM': '',
		   'POS': '',
		   'ID': '',
		   'REF': '',
		   'ALT': '',
		   'QUAL': '',
		   'FILTER': '',
		   'INFO': '',
		   'FORMAT': '',
		   'DETAILS': '',
		   'GENES': ''}

sep = '\t'
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
		gene_positions_map[gene['Symbol']] = {'chromosome': gene['ChrName'],
				                              'start': int(gene['Gene_start']),
				                              'end': int(gene['Gene_end'])}

def read_vcf(vcffile):
	f = open(vcffile, 'r')
	for line in f.readlines():
		if line.startswith('##'):
			continue
		yield line

def find_gene(position, chromosome):
	genes = list()
	for key, value in gene_positions_map.items():
		if chromosome == value['chromosome']:
			if position >= value['start'] and position <= value['end']:
				genes.append(key)
	return genes
"""
The headers in the vcf file are as follows:
1 -- #CHROM
2 -- POS
3 -- ID
4 -- REF
5 -- ALT
6 -- QUAL
7 -- FILTER
8 -- INFO
9 -- FORMAT
10 -- DETAILS
"""
def annotate_entries(vcffile, output):
	data = csv.DictReader(read_vcf(vcffile), delimiter='\t')
	count = 0
	for entry in data:
		count = count + 1
		ENTRY = copy.deepcopy(ENTRY_T)
		position, chromosome = [int(entry['POS']), entry['#CHROM']]
		genes = find_gene(position, chromosome)
		ENTRY['#CHROM'] = chromosome
		ENTRY['POS'] = entry['POS']
		ENTRY['ID'] = entry['ID']
		ENTRY['REF'] = entry['REF']
		ENTRY['ALT'] = entry['ALT']
		ENTRY['QUAL'] = entry['QUAL']
		ENTRY['FILTER'] = entry['FILTER']
		ENTRY['INFO'] = entry['INFO']
		ENTRY['FORMAT'] = entry['FORMAT']
		ENTRY['DETAILS'] = entry['DETAILS']

		for gene in genes:
			ENTRY['GENES'] = gene

			field_values = [str(ENTRY[i]) for i in HEADERS]
			output.write(sep.join(field_values))
			output.write('\n')
	return count

def main(vcf_folder_path, genesfile):

	create_gene_position_map(genesfile)

	all_files = os.listdir(vcf_folder_path)
	file_count = 0
	for file in all_files:
		start_time = tm.ctime(tm.time())
		print("FILE NAME:", file)
		file_count = file_count + 1
		filename = os.path.basename(file)
		output_file = filename.replace('.vcf', '_output.vcf')
		output = open(vcf_folder_path+output_file, 'w')
		output.write(sep.join(HEADERS))
		output.write('\n')

		variants_count = annotate_entries(vcf_folder_path+file, output)

		end_time = tm.ctime(tm.time())
		print("NUMBER OF VARIANTS PROCESSED:", variants_count)
		print("Start time:", start_time)
		print("End time:", end_time)

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
