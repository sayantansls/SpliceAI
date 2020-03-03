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
import twobitreader
import create_splice_matrices

"""
The headers in the strong_exons.tsv and rare_exons.tsv file are as follows:
1 -- Gene	
2 -- Chromosome	
3 -- Strand	
4 -- Exon_start	
5 -- Exon_end
"""
sep = '\t'

HEADERS = ['Gene',
		   'Chromosome',
		   'Strand',
		   'Exon_start',
		   'Exon_start_site_present',
		   'Exon_start_site_positions',
		   'Exon_end',
		   'Exon_end_site_present',
		   'Exon_end_site_positions',
		   'Comments']

ENTRY_T = {'Gene': '',
		   'Chromosome': '',
		   'Strand': '',
		   'Exon_start' : '',
		   'Exon_start_site_present': '',
		   'Exon_start_site_positions': '',
		   'Exon_end': '',
		   'Exon_end_site_present': '',
		   'Exon_end_site_positions': '',
		   'Comments': ''}

def load_genome_sequence():
	global genome
	genome = twobitreader.TwoBitFile('/home/sayantan/Desktop/hg19.2bit')

def get_bases(pos1, pos2, chrom):
	return genome[chrom][pos1-1:pos2]

def processing_entries(inputfile, output):
	for entry in utils.records_iterator(inputfile):
		ENTRY = copy.deepcopy(ENTRY_T)
		ENTRY['Gene'] = entry['Gene']
		chrom = entry['Chromosome']
		ENTRY['Chromosome'] = chrom
		strand = entry['Strand']
		ENTRY['Strand'] = strand
		exon_start, exon_end = [int(entry['Exon_start']), int(entry['Exon_end'])]
		ENTRY['Exon_start'] = exon_start
		ENTRY['Exon_end'] = exon_end
		if strand == '+':
			ENTRY['Exon_start_site_present'] = 'acceptor'
			ENTRY['Exon_end_site_present'] = 'donor'
			(start_site_A, start_site_G)  = (exon_start-2, exon_start-1)
			(end_site_G, end_site_T) = (exon_end+1, exon_end+2)
			(AG, GT) = (get_bases(start_site_A, start_site_G, chrom), get_bases(end_site_G, end_site_T, chrom)) 
			ENTRY['Exon_start_site_positions'] = (start_site_A, start_site_G)
			ENTRY['Exon_end_site_positions'] = (end_site_G, end_site_T)
		else:
			ENTRY['Exon_start_site_present'] = 'donor'
			ENTRY['Exon_end_site_present'] = 'acceptor'
			(start_site_G, start_site_T) = (exon_start-2, exon_start-1)
			(end_site_A, end_site_G) = (exon_end+1, exon_end+2)
			(AG_tmp, GT_tmp) = (get_bases(end_site_A, end_site_G, chrom), get_bases(start_site_G, start_site_T, chrom))
			(AG, GT) = (create_splice_matrices.create_reverse_complementary_sequence(AG_tmp), \
			create_splice_matrices.create_reverse_complementary_sequence(GT_tmp))
			ENTRY['Exon_start_site_positions'] = (start_site_G, start_site_T)
			ENTRY['Exon_end_site_positions'] = (end_site_A, end_site_G)

		if AG.upper() != 'AG':
			ENTRY['Comments'] = 'No acceptor splice (AG) site'
		elif GT.upper() != 'GT':
			ENTRY['Comments'] = 'No donor splice (GT) site'
		elif AG.upper() != 'AG' and GT.upper() != 'GT':
			ENTRY['Comments'] = 'Neither acceptor (AG) nor donor (GT) splice site present'

		#print(strand, AG, GT)
		field_values = [str(ENTRY[i]) for i in HEADERS]
		output.write(sep.join(field_values))
		output.write('\n')

def main(strong_exons_inputfile, rare_exons_inputfile, strong_exons_outputfile, rare_exons_outputfile):
	print("Start of code:", tm.ctime(tm.time()))

	strong_exons_output = open(strong_exons_outputfile, 'w')
	rare_exons_output = open(rare_exons_outputfile, 'w')

	load_genome_sequence()

	strong_exons_output.write(sep.join(HEADERS))
	strong_exons_output.write('\n')
	rare_exons_output.write(sep.join(HEADERS))
	rare_exons_output.write('\n')

	processing_entries(strong_exons_inputfile, strong_exons_output)
	processing_entries(rare_exons_inputfile, rare_exons_output)

	print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])