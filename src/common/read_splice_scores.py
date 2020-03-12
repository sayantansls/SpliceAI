"""
@author : sayantan
This script reads txt files containing the splice scores from the three predictors - MaxEntSxan(MES), 
SplicePort and NNSplice. It annotates the scores to the sequences file.
"""

import time as tm
import sys
import os

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
def read_maxentscan_scores(inputfile):


"""
The NNsplice scores file is as follows:
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
def read_nnsplice_scores(inputfile):


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
def read_assp_scores(inputfile):


def main(inputfile):
	print("Start of code:", tm.ctime(tm.time()))

	filename = os.path.basename(inputfile)
	if 'mes' in filename:
		read_maxentscan_scores(inputfile)
	elif 'nnsplice' in filename:
		read_nnsplice_scores(inputfile)
	elif 'assp' in filename:
		read_assp_scores(inputfile)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1])