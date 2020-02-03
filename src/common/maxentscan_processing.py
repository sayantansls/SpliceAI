# author -- sayantan
# The purpose of this code is as follows:
# 1 -- It consumes the output.tsv file and loads the sequences for MaxEntScan
# 2 -- Converts the sequences into fasta format to be used for bulk uploading on MaxEntScan

from __future__ import print_function
import time as tm 
import utils

def mes_processing(inputfile, donor_ref, donor_var, acceptor_ref, acceptor_var):
	for variant in utils.records_iterator(inputfile):
		fasta_line = "> " + variant['genomicHGVS']
		
		donor_ref.write(fasta_line)
		donor_ref.write('\n')
		donor_ref.write(variant['RefSeq_MES_donor'])
		donor_ref.write('\n')

		donor_var.write(fasta_line)
		donor_var.write('\n')
		donor_var.write(variant['VarSeq_MES_donor'])
		donor_var.write('\n')

		acceptor_ref.write(fasta_line)
		acceptor_ref.write('\n')
		acceptor_ref.write(variant['RefSeq_MES_acceptor'])
		acceptor_ref.write('\n')

		acceptor_var.write(fasta_line)
		acceptor_var.write('\n')
		acceptor_var.write(variant['VarSeq_MES_acceptor'])
		acceptor_var.write('\n')

def main(inputfile, donor_ref_file, donor_var_file, acceptor_ref_file, acceptor_var_file):
	print("Start of code:", tm.ctime(tm.time()))

	donor_ref = open(donor_ref_file, 'w')
	donor_var = open(donor_var_file, 'w')
	acceptor_ref = open(acceptor_ref_file, 'w')
	acceptor_var = open(acceptor_var_file, 'w')

	mes_processing(inputfile, donor_ref, donor_var, acceptor_ref, acceptor_var)
	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
