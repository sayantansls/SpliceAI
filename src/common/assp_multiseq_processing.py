# author -- sayantan
# The purpose of this code is as follows:
# 1 -- This code consumes assp_splice_variants.tsv file and converts the sequences into a FASTA format

from __future__ import print_function
import time as tm 
import utils
import re

def assp_sequence_processing(infile, ref_output, var_output):
	for variant in utils.records_iterator(infile):
		pos = ''.join(re.findall(r'[0-9]+', variant['genomicHGVS']))

		ref_output.write("> Variant- " + pos)
		ref_output.write('\n')
		ref_output.write(variant['RefSeq_ASSP_donor'])
		ref_output.write('\n')

		var_output.write("> Variant- " + pos)
		var_output.write('\n')
		var_output.write(variant['VarSeq_ASSP_donor'])
		var_output.write('\n')


def main(infile, ref_outfile, var_outfile):
	print("Start of code:", tm.ctime(tm.time()))
	ref_output = open(ref_outfile, 'w')
	var_output = open(var_outfile, 'w')
	assp_sequence_processing(infile, ref_output, var_output)

	print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2], sys.argv[3])