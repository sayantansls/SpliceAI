"""
@author: sayantan
This code consumes the files of exon_distribution.tsv files and converts them into 
fasta format for bulk upload and scoring in case of MaxEntScan, NNsplice and ASSP.
"""
from __future__ import print_function
import time as tm 
import utils
import os
import sys

def exon_processing_donor(inputfile, donor):
	for line in utils.records_iterator(inputfile):
		fasta_line = ">" + line['Gene'] + line['Exon_end']
		donor.write(fasta_line)
		donor.write('\n')
		donor.write(line['Donor_seq'])
		donor.write('\n')

def exon_processing_acceptor(inputfile, acceptor):
	for line in utils.records_iterator(inputfile):
		fasta_line = ">" + line['Gene'] + line['Exon_start']
		acceptor.write(fasta_line)
		acceptor.write('\n')
		acceptor.write(line['Acceptor_seq'])
		acceptor.write('\n')

def main(inputfile, donor_file, acceptor_file):
	print("Start of code:", tm.ctime(tm.time()))
	donor = open(donor_file, 'w')
	acceptor = open(acceptor_file, 'w')

	exon_processing_donor(inputfile, donor)
	exon_processing_acceptor(inputfile, acceptor)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3])