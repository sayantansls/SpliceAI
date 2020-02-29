"""
@author : sayantan
This script will consume genomic positions of essential splice sites (donor and acceptor) and give all possible mutations 
of the GT (donor) and AG(acceptor) consensus bases of these sites 
"""

from __future__ import print_function
import time as tm
import sys
import utils
import copy

"""
The headers in the strong_exons.tsv and rare_exons.tsv file are as follows:
1 -- Gene	
2 -- Chromosome	
3 -- Strand	
4 -- Exon_start	
5 -- Exon_end
"""

HEADERS = ['Gene',
		   'Chromosome',
		   'Strand',
		   'Exon_start',
		   'Exon_start_site_present',
		   'Exon_start_site_positions',
		   'Exon_end',
		   'Exon_end_site_present',
		   'Exon_end_site_positions']

ENTRY_T = {'Gene': '',
		   'Chromosome': '',
		   'Strand': '',
		   'Exon_start' : '',
		   'Exon_start_site_present': '',
		   'Exon_start_site_positions': '',
		   'Exon_end': '',
		   'Exon_end_site_present': '',
		   'Exon_end_site_positions': ''}

def processing_entries(inputfile, output):
	for entry in utils.records_iterator(inputfile):
		ENTRY = copy.deepcopy(ENTRY_T)
		ENTRY['Gene'] = entry['Gene']
		ENTRY['Chromosome'] = entry['Chromosome']
		strand = entry['Strand']
		ENTRY['Strand'] = strand
		exon_start, exon_end = [int(entry['Exon_start']), int(entry['Exon_end'])]
		ENTRY['Exon_start'] = exon_start
		ENTRY['Exon_end'] = exon_end
		if strand == '+':
			ENTRY['Exon_start_site_present'] = 'acceptor'
			ENTRY['Exon_end_site_present'] = 'donor'
		else:
			ENTRY['Exon_start_site_present'] = 'donor'
			ENTRY['Exon_end_site_present'] = 'acceptor'



def main(strong_exons_inputfile, rare_exons_inputfile, strong_exons_outputfile, rare_exons_outputfile):
	print("Start of code:", tm.ctime(tm.time()))

	strong_exons_output = open(strong_exons_outputfile, 'w')
	rare_exons_output = open(rare_exons_outputfile, 'w')

	strong_exons_output.write(HEADERS)
	strong_exons_output.write('\n')
	rare_exons_output.write(HEADERS)
	rare_exons_output.write('\n')

	processing_entries(strong_exons_inputfile, strong_exons_output)
	processing_entries(rare_exons_inputfile, rare_exons_output)

	print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])