# author -- sayantan
# This code consumes a file containing gene name and genomic HGVS and converts it into a VCF file 
# according to the format specified by Illumina/SpliceAI

import time as tm 
import twobitreader
import re
import copy
import sys
import utils
import vcf_generator_ver2

headers = ['#CHROM',
		   'POS',
		   'ID',
		   'REF',
		   'ALT',
		   'QUAL',
		   'FILTER',
		   'INFO']

ENTRY_T = {'#CHROM': '',
           'POS': '',
           'ID': '.',
           'REF': '',
           'ALT': '',
           'QUAL': '.',
           'FILTER': '.',
           'INFO': '.'}

sep = '\t'

def create_substitution_entries(output):
	for variant in subs:
		ENTRY = copy.deepcopy(ENTRY_T)
		gene, gHGVS = [variant[0], variant[1]]
		ref_alt = re.findall(r'[A-Z]', gHGVS)
		position = ''.join(re.findall(r'[0-9]+', gHGVS))
		ref, alt = [ref_alt[0], ref_alt[1]]

		ENTRY['#CHROM'] = vcf_generator_ver2.get_chromosome(gene).replace('chr','')
		ENTRY['POS'] = position
		ENTRY['REF'] = ref
		ENTRY['ALT'] = alt

		field_values = [ENTRY[i] for i in headers]
		output.write(sep.join(field_values))
		output.write('\n')

def create_non_substitution_entries(output):
	for variant in others:
		ENTRY = copy.deepcopy(ENTRY_T)
		gene, gHGVS = [variant[0], variant[1]]
		chrom = get_chromosome(gene)

		if "del" in gHGVS:
			(pos, ref, alt) = vcf_generator_ver2.del_handling(gHGVS, chrom)
		else:
			(pos, ref, alt) = vcf_generator_ver2.other_handling(gHGVS, chrom)

		ENTRY['#CHROM'] = vcf_generator_ver2.get_chromosome(gene).replace('chr','')
		ENTRY['POS'] = pos
		ENTRY['REF'] = ref
		ENTRY['ALT'] = alt

		field_values = [str(ENTRY[i]) for i in headers]
		output.write(sep.join(field_values))
		output.write('\n')

def main(input_file, output_file, genome_file, input_template, genesfile):
	print("Start of code:", tm.ctime(tm.time()))
	
	vcf_generator_ver2.load_human_genome_sequence(genome_file)
	vcf_generator_ver2.load_metadata(input_template)

	output = open(output_file, 'w')
	for line in metadata:
		output.write(line)

	output.write(sep.join(headers))
	output.write('\n')

	vcf_generator_ver2.create_gene_chromosome_map(genesfile)
	vcf_generator_ver2.process_input_file(input_file)

	print("Total variants in the file:", len(variants_list))
	
	vcf_generator_ver2.segregate_variants()
	
	print("#Substitution variants:", len(subs))
	print("#INDEL variants:", len(others))

	create_substitution_entries(output)
	create_non_substitution_entries(output)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])