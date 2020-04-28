"""
@author : sayantan
This script consumes the genelist dump from StrandOmics and adds certain annotations to it.
"""

from __future__ import print_function
import utils
import time as tm
import copy, sys, os

sep = '\t'

HEADERS = ['gene_list_id',
		   'gene_list_name',
		   'gene_list_display_name',
		   'gene_list_description',
		   'gene_list_status',
		   'gene_symbol',
		   'chromosome',
		   'strand']

ENTRY_T = {'gene_list_id': '',
		   'gene_list_name': '',
		   'gene_list_display_name': '',
		   'gene_list_description': '',
		   'gene_list_status': '',
		   'gene_symbol': '',
		   'chromosome': '',
		   'strand': ''}

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
def create_gene_details_map(genesfile):
	global gene_details
	gene_details = dict()
	for gene in utils.records_iterator(genesfile):
		if not gene['Strand_gene_id'] in gene_details:
			gene_details[gene['Strand_gene_id']] = {'gene_name': '',
												    'chromosome': '',
					                                'strand': ''}
		gene_details[gene['Strand_gene_id']] = {'gene_name': gene['Symbol'],
												'chromosome': gene['ChrName'],
												'strand': gene['Strand']}

"""
The headers in the gene_list.tsv file are as follows:
1 -- gene_list_id	
2 -- gene_list_name	
3 -- gene_list_display_name	
4 -- gene_list_description	
5 -- gene_list_status	
6 -- creation_time	
7 -- strand_gene_id
"""

def annotate_gene_lists(inputfile, output):
	for entry in utils.records_iterator(inputfile):
		ENTRY = copy.deepcopy(ENTRY_T)

		strand_gene_id = entry['strand_gene_id']

		ENTRY['gene_list_id'] = entry['gene_list_id']
		ENTRY['gene_list_name'] = entry['gene_list_name']
		ENTRY['gene_list_display_name'] = entry['gene_list_display_name']
		ENTRY['gene_list_description'] = entry['gene_list_description']
		ENTRY['gene_list_status'] = entry['gene_list_status']

		try:
			ENTRY['gene_symbol'] = gene_details[strand_gene_id]['gene_name']
			ENTRY['chromosome'] = gene_details[strand_gene_id]['chromosome']
			ENTRY['strand'] = gene_details[strand_gene_id]['strand']
		except:
			ENTRY['gene_symbol'] = "Could not map gene symbol:" + strand_gene_id

		field_values = [str(ENTRY[i]) for i in HEADERS]
		output.write(sep.join(field_values))
		output.write('\n')

def main(inputfile, genesfile):
	print("Start of code:", tm.ctime(tm.time()))

	inputfile_name = os.path.basename(inputfile)
	inputfile_path = "/home/sayantan/Desktop/SpliceAI/data/strandomics_input_data"

	outputfile_name = inputfile_name.replace('.tsv','_annotated.tsv')

	output = open(os.path.join(inputfile_path,outputfile_name), 'w')
	output.write(sep.join(HEADERS))
	output.write('\n') 

	create_gene_details_map(genesfile)
	annotate_gene_lists(inputfile, output)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
