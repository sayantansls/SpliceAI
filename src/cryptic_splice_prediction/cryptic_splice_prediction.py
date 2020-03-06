"""
@author : sayantan
This code will consume the splice_variants.tsv file and search for and get the cryptic splice sites
near the variant on a search window range of +80 to -80 bases
"""

from __future__ import print_function
import utils
import time as tm
import re

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

def get_neighbourhood_region(position):
	region = {'start': (position - 80),
			  'end': (position + 80)}
	return region



def main():
	print("Start of code:", tm.ctime(tm.time()))

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main()