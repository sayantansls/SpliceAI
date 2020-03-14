"""
@author : sayantan
This code will consume the splice_variants.tsv file and search for and get the cryptic splice sites
near the variant on a search window range of +80 to -80 bases
"""

from __future__ import print_function
import utils
import time as tm
import re, sys
import create_splice_matrices
from maxentpy import maxent
from pprint import pprint

"""
The headers in the splice_variants.tsv file are as follows:
1 -- Source	
2 -- Case ID	
3 -- UK1	
4 -- UK2	
5 -- Panel	
6 -- Test	
7 -- Time and Date	
8 -- ClinicalManifestation	
9 -- Gender	
10 -- Ethnicity	
11 -- Age	
12 -- ReportWiseGeneListInfo
13 -- UK3	
14 -- UK4	
15 -- UK5	
16 -- Chromosome	
17 -- Gene	
18 -- pHGVS	
19 -- Genomic HGVS	
20 -- Transcript	
21 -- VariantType
22 -- AlleleFrequency	
23 -- Global PPDB	
24 -- Local PPDB	
25 -- Zygosity	
26 -- Clinvar IDs	
27 -- Variant Label Reason	
28 -- UK6	
29 -- rsID	
30 -- EVS	
31 -- ExAc	
32 -- dbSNP	
33 -- 1000 Genomes	
34 -- HGMD ID	
35 -- Bioinfo Summary	
36 -- Literature Summary
"""
END_CUTOFF = 6

def get_neighbourhood_region(position):
	region = {'start': (position - 80),
			  'end': (position + 80)}
	return region

def process_entries(inputfile):
	for entry in utils.records_iterator(inputfile):
		position = int(''.join(re.findall(r'[0-9]+', entry['Genomic HGVS'])))
		search_window = get_neighbourhood_region(position)
		chrom = entry['Chromosome']
		#strand = create_splice_matrices.get_strand_sense(entry['Gene'])
		neighbourhood_sequence = create_splice_matrices.get_original_sequence(search_window, chrom)

		donor_seq_list, acceptor_seq_list = [list(), list()]
		# Handling built only for positive strand
		# TODO -- build handling for negative strand
		for i in range(len(neighbourhood_sequence) - END_CUTOFF):
			if neighbourhood_sequence[i] == 'G' and neighbourhood_sequence[i+1] == 'T':
				donor_seq_list.append(neighbourhood_sequence[i-3:i+6])
			if neighbourhood_sequence[i] == 'A' and neighbourhood_sequence[i+1] == 'G':
				acceptor_seq_list.append(neighbourhood_sequence[i-18:i+5])

		donor_seq_list = list(filter(None, donor_seq_list))
		acceptor_seq_list = list(filter(None, acceptor_seq_list))

		donor_sites_with_scores, acceptor_sites_with_scores = [dict(), dict()]
		for donor in donor_seq_list:
			donor_sites_with_scores[donor] = maxent.score5(donor)
		for acceptor in acceptor_seq_list:
			acceptor_sites_with_scores[acceptor] = maxent.score3(acceptor)

		pprint(donor_sites_with_scores)
		pprint(acceptor_sites_with_scores)


def main(inputfile):
	print("Start of code:", tm.ctime(tm.time()))

	process_entries(inputfile)
	
	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1])