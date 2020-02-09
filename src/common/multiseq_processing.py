# author -- sayantan
# The purpose of this code is as follows:
# 1 -- It consumes the output.tsv file and loads the sequences for NNSplice, ASSP and MaxEntScan
# 2 -- Converts the sequences into fasta format to be used for bulk uploading on the respective_splice_predictors

from __future__ import print_function
import time as tm 
import utils
import os

def multiseq_processing(inputfile, ref_output, var_output):
	for variant in utils.records_iterator(inputfile):
		fasta_line = "> " + variant['gene'] + variant['genomicHGVS']
		
		ref_output.write(fasta_line)
		ref_output.write('\n')
		ref_output.write(variant['RefSeq'])
		ref_output.write('\n')

		var_output.write(fasta_line)
		var_output.write('\n')
		var_output.write(variant['VarSeq'])
		var_output.write('\n')

def mes_processing(inputfile, ref_output, var_output):
	donor_list, acceptor_list = [list(), list()]
	for variant in utils.records_iterator(inputfile):
		if variant['MES_Seq_used'] == 'Donor':
			donor_list.append(variant)
		else:
			acceptor_list.append(variant)

	ref_output.write("Donor sequences:")
	ref_output.write('\n')
	var_output.write("Donor sequences:")
	var_output.write('\n')

	for entry in donor_list:
		fasta_line = "> " + entry['gene'] + " " + entry['genomicHGVS']
		ref_output.write(fasta_line)
		ref_output.write('\n')
		ref_output.write(entry['RefSeq'])
		ref_output.write('\n')

		var_output.write(fasta_line)
		var_output.write('\n')
		var_output.write(entry['VarSeq'])
		var_output.write('\n')

	ref_output.write('\n')
	ref_output.write("Acceptor sequences:")
	ref_output.write('\n')

	var_output.write('\n')
	var_output.write("Acceptor sequences:")
	var_output.write('\n')

	for entry in acceptor_list:
		fasta_line = ">" + entry['gene'] + " " + entry['genomicHGVS']
		ref_output.write(fasta_line)
		ref_output.write('\n')
		ref_output.write(entry['RefSeq'])
		ref_output.write('\n')

		var_output.write(fasta_line)
		var_output.write('\n')
		var_output.write(entry['VarSeq'])
		var_output.write('\n')

def main(inputfile, ref_output_file, var_output_file):
	print("Start of code:", tm.ctime(tm.time()))

	ref_output = open(ref_output_file, 'w')
	var_output = open(var_output_file, 'w')

	filename = os.path.basename(inputfile)

	if "MaxEntScan" in filename:
		mes_processing(inputfile, ref_output, var_output)
	else:
		multiseq_processing(inputfile, ref_output, var_output)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	import sys
	main(sys.argv[1], sys.argv[2], sys.argv[3])