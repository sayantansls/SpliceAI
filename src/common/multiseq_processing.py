# author -- sayantan
# The purpose of this code is as follows:
# 1 -- It consumes the output.tsv file and loads the sequences for NNSplice and ASSP
# 2 -- Converts the sequences into fasta format to be used for bulk uploading on NNSplice and ASSP

from __future__ import print_function
import time as tm 
import utils

def multiseq_processing(inputfile, ref_output, var_output):
	for variant in utils.records_iterator(inputfile):
		fasta_line = "> " + variant['genomicHGVS']
		
		donor_ref.write(fasta_line)
		donor_ref.write('\n')
		donor_ref.write(variant['RefSeq'])
		donor_ref.write('\n')

		donor_var.write(fasta_line)
		donor_var.write('\n')
		donor_var.write(variant['VarSeq'])
		donor_var.write('\n')


def main(inputfile, ref_output_file, var_output_file):
	print("Start of code:", tm.ctime(tm.time()))

	ref_output = open(ref_output_file, 'w')
	var_output = open(var_output_file, 'w')

	mes_processing(inputfile, ref_output, var_output)
	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2], sys.argv[3])