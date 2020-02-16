"""
@author : sayantan
The purpose of this code is to consume a file containing the start and end positions of strong
exons and get the donor and acceptor splice sequences for them
"""

from __future__ import print_function
import time as tm 
import utils
import re
import sys
import copy
import create_splice_matrices

sep = '\t'

HEADERS = ['Gene',
		   'Chromosome',
		   'Strand',
		   'Exon_start',
		   'MES_acceptor_seq',
		   'NNSplice_acceptor_seq',
		   'ASSP/HSF_acceptor_seq',
		   'Exon_end',
		   'MES_donor_seq',
		   'NNSplice_donor_seq',
		   'ASSP/HSF_donor_seq']

ENTRY_T = {'Gene': '',
		   'Chromosome': '',
		   'Strand': '',
		   'Exon_start': '',
		   'MES_acceptor_seq': '',
		   'NNSplice_acceptor_seq': '',
		   'ASSP/HSF_acceptor_seq': '',
		   'Exon_end': '',
		   'MES_donor_seq': '',
		   'NNSplice_donor_seq': '',
		   'ASSP/HSF_donor_seq': ''}

# The headers in the exons_distribution.tsv file are as follows:
# 1 -- Gene
# 2 -- Chromosome	
# 3 -- Strand	
# 4 -- Exon Start	
# 5 -- Exon End

def process_input_file(input_file, output):
	for gene in utils.records_iterator(input_file):
		ENTRY = copy.deepcopy(ENTRY_T)

		ENTRY['Gene'] = gene['Gene']
		chrom = gene['Chromosome']
		ENTRY['Chromosome'] = chrom
		strand = gene['Strand']
		ENTRY['Strand'] = strand

		exon_start, exon_end = [int(gene['Exon_start']), int(gene['Exon_end'])]
		ENTRY['Exon_start'], ENTRY['Exon_end'] = [exon_start, exon_end] 
		
		if strand == '+':
			mes_donor = create_splice_matrices.create_mes_donor_range(exon_end, strand)
			mes_acceptor = create_splice_matrices.create_mes_acceptor_range(exon_start, strand)
			nnplice_donor = create_splice_matrices.create_nnsplice_donor_range(exon_end)
			nnsplice_acceptor = create_splice_matrices.create_nnsplice_acceptor_range(exon_start)
			assp_hsf_donor = create_splice_matrices.create_assp_hsf_range(exon_end)
			assp_hsf_acceptor = create_splice_matrices.create_assp_hsf_range(exon_start)
		else:
			mes_donor = create_splice_matrices.create_mes_donor_range(exon_start, strand)
			mes_acceptor = create_splice_matrices.create_mes_acceptor_range(exon_end, strand)
			nnplice_donor = create_splice_matrices.create_nnsplice_donor_range(exon_start)
			nnsplice_acceptor = create_splice_matrices.create_nnsplice_acceptor_range(exon_end)
			assp_hsf_donor = create_splice_matrices.create_assp_hsf_range(exon_start)
			assp_hsf_acceptor = create_splice_matrices.create_assp_hsf_range(exon_end)


		mes_donor_seq = create_splice_matrices.get_original_sequence(mes_donor, chrom)
		mes_acceptor_seq = create_splice_matrices.get_original_sequence(mes_acceptor, chrom)
		nnsplice_donor_seq = create_splice_matrices.get_original_sequence(nnplice_donor, chrom)
		nnsplice_acceptor_seq = create_splice_matrices.get_original_sequence(nnsplice_acceptor, chrom)
		assp_hsf_donor_seq = create_splice_matrices.get_original_sequence(assp_hsf_donor, chrom)
		assp_hsf_acceptor_seq = create_splice_matrices.get_original_sequence(assp_hsf_acceptor, chrom)

		if strand == '+':
			ENTRY['MES_donor_seq'] = mes_donor_seq
			ENTRY['MES_acceptor_seq'] = mes_acceptor_seq
			ENTRY['NNSplice_donor_seq'] = nnsplice_donor_seq
			ENTRY['NNSplice_acceptor_seq'] = nnsplice_acceptor_seq
			ENTRY['ASSP/HSF_donor_seq'] = assp_hsf_donor_seq
			ENTRY['ASSP/HSF_acceptor_seq'] = assp_hsf_acceptor_seq

		elif strand == '-':
			ENTRY['MES_donor_seq'] = create_splice_matrices.create_reverse_complementary_sequence(mes_donor_seq)
			ENTRY['MES_acceptor_seq'] = create_splice_matrices.create_reverse_complementary_sequence(mes_acceptor_seq)
			ENTRY['NNSplice_donor_seq'] = create_splice_matrices.create_reverse_complementary_sequence(nnsplice_donor_seq)
			ENTRY['NNSplice_acceptor_seq'] = create_splice_matrices.create_reverse_complementary_sequence(nnsplice_acceptor_seq)
			ENTRY['ASSP/HSF_donor_seq'] = create_splice_matrices.create_reverse_complementary_sequence(assp_hsf_donor_seq)
			ENTRY['ASSP/HSF_acceptor_seq'] = create_splice_matrices.create_reverse_complementary_sequence(assp_hsf_acceptor_seq)

		field_values = [str(ENTRY[i]) for i in HEADERS]
		output.write(sep.join(field_values))
		output.write('\n')

def main(input_file, output_file):
	print("Start of code:", tm.ctime(tm.time()))
	output = open(output_file, 'w')
	output.write(sep.join(HEADERS))
	output.write('\n')

	process_input_file(input_file, output)

	print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
