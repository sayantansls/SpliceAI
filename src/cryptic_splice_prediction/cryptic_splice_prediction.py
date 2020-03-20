"""
@author : sayantan
This code will consume the splice_variants.tsv file and search for and get the
cryptic splice sites near the variant on a search window range of +80 to -80 bases
"""

from __future__ import print_function
import deepdiff, json
import time as tm
import re, sys, utils
import twobitreader
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

genome = twobitreader.TwoBitFile('/home/sayantan/Desktop/hg19.2bit')

class SpliceSite(object):
	"""SpliceSite defines a site to be investigated for splicing,
	it has start and end coordinates and the sequence"""

	def __init__(self, start, end, sequence):
		sups SpliceSite, self).__init__()
		self.start = start
		self.end = end
		self.sequence = sequence

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

def get_reference_sequence(window, position, chrom):
    sequence = genome[chrom][window['start'] - 1:position - 1] + genome[chrom][position - 1] + genome[chrom][position:window['end']]
    return sequence

def get_variant_sequence(window, position, chrom, alt):
	sequence = genome[chrom][window['start'] - 1:position - 1] + alt + genome[chrom][position:window['end']]
	return sequence

def annotate_scores(donor_list, acceptor_list):
	donor_with_scores, acceptor_with_scores = [dict(), dict()]
	for donor in donor_list:
		donor_with_scores[donor] = maxent.score5(donor)
	for acceptor in acceptor_list:
		acceptor_with_scores[acceptor] = maxent.score3(acceptor)
	return donor_with_scores, acceptor_with_scores

def dict_differ(ref_dict, var_dict):
	ref_keys, var_keys = [ref_dict.keys(), var_dict.keys()]
	matched_keys = set(ref_keys).intersection(set(var_keys))
	unmatched_keys = set(ref_keys).symmetric_difference(set(var_keys))

	for key in matched_keys:
		if ref_dict[key] != var_dict[key]:
			unmatched_keys.add(key)
	diff_dict = dict()
	for key in unmatched_keys:
		if key in ref_dict and key not in var_dict:
			diff_dict[key] = ref_dict[key]
		elif key in var_dict and key not in ref_dict:
			diff_dict[key] = var_dict[key]
		elif key in ref_dict and key in var_dict:
			diff_dict[key] = {'reference' : ref_dict[key], 'variant': var_dict[key]}
	return diff_dict

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
		ref_neighbourhood_sequence = get_reference_sequence(search_window, position, chrom)
		var_neighbourhhod_sequence = get_variant_sequence(search_window, position, chrom, alt)

		(donor_ref_list, acceptor_ref_list) = populate_sequences(ref_neighbourhood_sequence, strand)
		(donor_var_list, acceptor_var_list) = populate_sequences(var_neighbourhhod_sequence, strand)

		(donor_ref_sites_with_scores, acceptor_ref_sites_with_scores) = annotate_scores(donor_ref_list, acceptor_ref_list)
		(donor_var_sites_with_scores, acceptor_var_sites_with_scores) = annotate_scores(donor_var_list, acceptor_var_list)

		donor_diff = dict_differ(donor_ref_sites_with_scores, donor_var_sites_with_scores)
		acceptor_diff = dict_differ(acceptor_ref_sites_with_scores, acceptor_var_sites_with_scores)

		print(entry['Genomic HGVS'])
		#pprint(donor_ref_sites_with_scores)
		#pprint(acceptor_ref_sites_with_scores)
		#print('\n')
		#pprint(donor_var_sites_with_scores)
		#pprint(acceptor_var_sites_with_scores)
		print("Donor diff")
		pprint(donor_diff)
		print("Acceptor diff")
		pprint(acceptor_diff)

def main(inputfile, genesfile):
	print("Start of code:", tm.ctime(tm.time()))

	create_genes_strand_map(genesfile)
	process_entries(inputfile)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2])
