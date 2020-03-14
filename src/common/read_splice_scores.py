"""
@author : sayantan
This script reads txt files containing the splice scores from the three predictors - MaxEntSxan(MES), 
SplicePort and NNSplice. It annotates the scores to the sequences file.
"""
from __future__ import print_function
import time as tm
import sys, os
import re, copy

HEADERS = ['Gene_name',
	       'genomic_HGVS',
	       'sequence',
	       'predictor',
	       'splice_site_type',
	       'score']

ENTRY_T = {'Gene_name': '',
	       'genomic_HGVS': '',
	       'sequence': '',
	       'predictor': '',
	       'splice_site_type': '',
	       'score': ''}
"""
The MaxEntScan scores file is as follows:
> D2HGDH g.242689707A>G
aAGGTACTG	MAXENT: 8.56	
> SCN1A g.166895930T>A
TTGGTaAAA	MAXENT: 3.23	
> PLA2G6 g.38528838C>T
TCgGTGAGC	MAXENT: 9.10	
> ALS2 g.202590081G>A
cAGGTAACA	MAXENT: 8.88	
"""
def read_maxentscan_scores(file, filename):
	for line in file.readlines():
		if line.startswith(">"):
			ENTRY = copy.deepcopy(ENTRY_T)
			gene, gHGVS, seq, score = ['','','','']

		if line.startswith(">"):
			gene_gHGVS = line.replace("> ", "")
			gene, gHGVS = gene_gHGVS.split()
		else:
			seq_score = line.replace("MAXENT:","")
			seq, score = seq_score.split()
			
		ENTRY['Gene_name'] = gene
		ENTRY['genomic_HGVS'] = gHGVS
		ENTRY['sequence'] = seq
		ENTRY['predictor'] = 'MaxEntScan'
		ENTRY['score'] = score

		if 'donor' in filename:
			ENTRY['splice_site_type'] = 'donor'
		elif 'acceptor' in filename:
			ENTRY['splice_site_type'] = 'acceptor'
		else:
			print("Wrong file name")

		field_values = [ENTRY[i] for i in HEADERS]
		if not '' in field_values:
			print(field_values)

"""
The NNSplice scores file is as follows:
Acceptor site predictions for TPP1 :
Start   End    Score     Intron               Exon
    1    41     0.94     ccctgacccctgaccctacaggctgcccccaggctgggtgt

Donor site predictions for LAMA2 :
Start   End    Score     Exon   Intron
   25    39     0.86     tctcgctgtgagtgt

Acceptor site predictions for LAMA2 :
Start   End    Score     Intron               Exon
    1    41     0.95     ttcctccctctttttgactagaaatctcgctgtgagtgtga

Donor site predictions for LAMA2 :
Start   End    Score     Exon   Intron
    1    15     0.99     ttgtaaggtgagtga
"""
def read_nnsplice_scores(file):
	for line in file.readlines():
		ENTRY = copy.deepcopy(ENTRY_T)
		site_type, gene, score, seq = ['','','','']

		ENTRY['predictor'] = 'NNSplice'
		line_elements = line.split()
		if line_elements[0] == 'Donor' or line_elements[0] == 'Acceptor':
			ENTRY['splice_site_type'] = line_elements[0]
			ENTRY['Gene_name'] = line_elements[-2]
		elif line_elements[0] == 'Start':
			pass
		elif int(line_elements[0]).isdigit():
			ENTRY['sequence'] = line_elements[-1]
			ENTRY['score'] = line_elements[-2]
		elif line_elements == []:
			pass

		field_values = [ENTRY[i] for i in HEADERS]
		print(field_values)

"""
Sequence: G
  Position|Putative splice site|SequenceScore|Intron GC|Activation alt./cryptic|Activation constitutive|Confidence
  38|Alt. isoform/cryptic donor|TTGATCCTGGgtctgtggga|4.859|0.543|0.927|0.053|0.943
  42|Alt. isoform/cryptic donor|TCCTGGGTCTgtgggattct|5.371|0.543|0.933|0.048|0.949
  71|Constitutive acceptor|tgtgttgcagATCACAGGGA|7.033|0.543|0.309|0.676|0.543
  78|Alt. isoform/cryptic acceptor|cagatcacagGGATGACCTG|4.124|0.514|0.565|0.416|0.265
Sequence: A
  Position|Putative splice site|SequenceScore|Intron GC|Activation alt./cryptic|Activation constitutive|Confidence
  70|Alt. isoform/cryptic donor|TATTTTCACAgtatgttgtt|8.326|0.400|0.735|0.199|0.729
  72|Alt. isoform/cryptic acceptor|attttcacagTATGTTGTTG|9.939|0.300|0.565|0.425|0.248
  100|Alt. isoform/cryptic acceptor|tgtctgaaagCTGGTTTGGA|3.503|0.300|0.803|0.191|0.763

"""
def read_assp_scores(file):
	for line in file.readlines():
		print(line)


def main(inputfile):
	print("Start of code:", tm.ctime(tm.time()))

	file = open(inputfile, 'r')

	filename = os.path.basename(inputfile)
	if 'mes' in filename:
		read_maxentscan_scores(file, filename)
	elif 'nnsplice' in filename:
		read_nnsplice_scores(file)
	elif 'assp' in filename:
		read_assp_scores(file)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1])