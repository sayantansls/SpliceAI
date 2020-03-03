"""
@author : sayantan
This script consumes the output files from the 'all_possible_snv_mutations.py' and converts it into a medium
for feeding as input to the vcf_generator_script
"""

from __future__ import print_function
import time as tm
import utils
import sys
import copy
import re
import all_possible_snv_mutations

sep = '\t'

HEADERS = ['Gene name',
		   'genomicHGVS',
		   'Site',
		   'Strand']

ENTRY_T = {'Gene name': '',
           'genomicHGVS': '',
           'Site': '',
           'Strand': ''}

"""
The headers in the inputfile are as follows:
1 -- Gene	
2 -- Chromosome	
3 -- Strand	
4 -- Exon_start	
5 -- Exon_start_site_present	
6 -- Exon_start_site_positions	
7 -- Exon_end	
8 -- Exon_end_site_present
9 -- Exon_end_site_positions	
10 -- Comments
"""
def get_positions(string):
	tmp_list = string.split(',')
	string_list = [''.join(re.findall(r'[0-9]+', i)) for i in tmp_list]
	positions = [int(i) for i in string_list]

	return positions

def process_entries(inputfile, output):
	for entry in utils.records_iterator(inputfile):
		if entry['Comments'] == '':
			ENTRY = copy.deepcopy(ENTRY_T)
			ENTRY['Gene name'] = entry['Gene']
			ENTRY['Strand'] = entry['Strand']
			positions_start = get_positions(entry['Exon_start_site_positions'])
			positions_end = get_positions(entry['Exon_end_site_positions'])
			chrom = entry['Chromosome']

			for position in positions_start:
				ENTRY['genomicHGVS'] = str(position) + all_possible_snv_mutations.get_bases(position, position, chrom) + ">N"
				ENTRY['Site'] = entry['Exon_start_site_present']

				field_values = [str(ENTRY[i]) for i in HEADERS]
				output.write(sep.join(field_values))
				output.write('\n') 

			for position in positions_end:
				ENTRY['genomicHGVS'] = str(position) + all_possible_snv_mutations.get_bases(position, position, chrom) + ">N"
				ENTRY['Site'] = entry['Exon_end_site_present']

				field_values = [str(ENTRY[i]) for i in HEADERS]
				output.write(sep.join(field_values))
				output.write('\n')

def main(inputfile, outputfile):
	print("Start of code:", tm.ctime(tm.time()))

	all_possible_snv_mutations.load_genome_sequence()
	output = open(outputfile, 'w')
	output.write(sep.join(HEADERS))
	output.write('\n')
	process_entries(inputfile, output)

	print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])