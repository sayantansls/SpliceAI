"""
@author : sayantan
This code will consume the splice_variants.tsv file and search for and get the
cryptic splice sites near the variant on a search window range of +80 to -80 bases
"""

from __future__ import print_function
import utils
import time as tm
import re, sys
import create_splice_matrices
from maxentpy import maxent
from pprint import pprint

"""
The headers in the genes.tsv file are as follows:
1 -- Strand_gene_id
2 -- ChrName
3 -- Strand
4 -- Gene_start
5 -- Gene_end
6 -- Symbol
7 -- Entrez_id
"""
def create_genes_strand_map(genesfile):
	global genes_strand_map
	genes_strand_map = dict()
	for gene in utils.records_iterator(genesfile):
		genes_strand_map[gene['Symbol']] = gene['Strand']

def get_neighbourhood_region(position):
	region = {'start': (position - 80),
			  'end': (position + 80)}
	return region

def get_mes_donor_sequence(position, strand, neighbourhood_sequence):
	if strand == '+':
		sequence = neighbourhood_sequence[position - 3: position + 6]
	elif strand == '-':
		sequence = neighbourhood_sequence[position - 4: position + 5]
	return sequence

def get_mes_acceptor_sequence(position, strand, neighbourhood_sequence):
	if strand == '+':
		sequence = neighbourhood_sequence[position - 18: position + 5]
	elif strand == '-':
		sequence = neighbourhood_sequence[position - 3: position + 20]
	return sequence

def populate_sequences(neighbourhood_sequence, strand):
	donor_seq_list, acceptor_seq_list = [list(), list()]
	if strand == '+':
		CUTOFF = 6
		donor_base1, donor_base2 = ['G', 'T']
		acceptor_base1, acceptor_base2 = ['A', 'G']
	elif strand == '-':
		CUTOFF = 20
		donor_base1, donor_base2 = ['A', 'C']
		acceptor_base1, acceptor_base2 = ['C', 'T']

	for i in range(len(neighbourhood_sequence) - CUTOFF):
			if neighbourhood_sequence[i] == donor_base1 and neighbourhood_sequence[i+1] == donor_base2:
				donor_seq_list.append(get_mes_donor_sequence(i, strand, neighbourhood_sequence))
			if neighbourhood_sequence[i] == acceptor_base1 and neighbourhood_sequence[i+1] == acceptor_base2:
				acceptor_seq_list.append(get_mes_acceptor_sequence(i, strand, neighbourhood_sequence))

	donor_seq_list = list(filter(None, donor_seq_list))
	acceptor_seq_list = list(filter(None, acceptor_seq_list))

	if strand == '-':
		donor_seq_list = [create_splice_matrices.create_reverse_complementary_sequence(i) for i in donor_seq_list]
		acceptor_seq_list = [create_splice_matrices.create_reverse_complementary_sequence(i) for i in acceptor_seq_list]

	return donor_seq_list, acceptor_seq_list

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
def process_entries(inputfile):
	for entry in utils.records_iterator(inputfile):
		position = int(''.join(re.findall(r'[0-9]+', entry['Genomic HGVS'])))
		search_window = get_neighbourhood_region(position)
		chrom = entry['Chromosome']
		strand = genes_strand_map[entry['Gene']]
		(ref, alt) = create_splice_matrices.get_ref_alt(entry['Genomic HGVS'])
		print(entry['Genomic HGVS'])
		neighbourhood_sequence = create_splice_matrices.get_original_sequence(search_window, chrom)
		ref_neighbourhood_sequence = create_splice_matrices.get_reference_sequence(search_window, position, chrom)
		var_neighbourhhod_sequence = create_splice_matrices.get_reference_sequence(search_window, position, chrom)

		print(neighbourhood_sequence == ref_neighbourhood_sequence)
		print(neighbourhood_sequence == var_neighbourhhod_sequence)
		(donor_ref_list, acceptor_ref_list) = populate_sequences(neighbourhood_sequence, strand)

		donor_sites_with_scores, acceptor_sites_with_scores = [dict(), dict()]
		for donor in donor_ref_list:
			donor_sites_with_scores[donor] = maxent.score5(donor)
		for acceptor in acceptor_ref_list:
			acceptor_sites_with_scores[acceptor] = maxent.score3(acceptor)

		pprint(donor_sites_with_scores)
		pprint(acceptor_sites_with_scores)


def main(inputfile, genesfile):
	print("Start of code:", tm.ctime(tm.time()))

	create_genes_strand_map(genesfile)
	process_entries(inputfile)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
