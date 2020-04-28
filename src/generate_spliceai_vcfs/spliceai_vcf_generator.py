"""
@author -- sayantan
This code consumes a file containing gene name and genomic HGVS and converts it into a VCF file
according to the format specified by Illumina/SpliceAI
"""
from __future__ import print_function
import time as tm
import twobitreader
import re, copy, 	sys
import utils, os
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

def get_all_alt(base):
	org_bases = ['A', 'T', 'G', 'C']
	bases = copy.deepcopy(org_bases)
	if base in bases:
		bases.remove(base)
	return ','.join(bases)

def create_substitution_entries(output, subs):
	for variant in subs:
		ENTRY = copy.deepcopy(ENTRY_T)
		gene, gHGVS = [variant[0], variant[1]]
		ref_alt = re.findall(r'[A-Za-z]', gHGVS)
		position = ''.join(re.findall(r'[0-9]+', gHGVS))
		ref, alt = [ref_alt[0], ref_alt[1]]

		ENTRY['#CHROM'] = vcf_generator_ver2.get_chromosome(gene).replace('chr','')
		ENTRY['POS'] = position
		ENTRY['REF'] = ref
		ENTRY['ALT'] = alt
		ENTRY['ALT'] = get_all_alt(ref.upper())

		field_values = [ENTRY[i] for i in headers]
		output.write(sep.join(field_values))
		output.write('\n')

def create_non_substitution_entries(output, others):
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

def process_entries(vcf_file, output):
    for variant in utils.records_iterator(vcf_file):
        ENTRY = copy.deepcopy(ENTRY_T)

        ENTRY['#CHROM'] = variant['CHROM'].replace('chr','')
        ENTRY['POS'] = variant['POS']
        ENTRY['REF'] = variant['REF']
        ENTRY['ALT'] = variant['ALT']

        field_values = [str(ENTRY[i]) for i in headers]
        output.write(sep.join(field_values))
        output.write('\n')

def main(input_files_path, genome_file, input_template, genesfile):
	print("Start of code:", tm.ctime(tm.time()))

	vcf_generator_ver2.load_human_genome_sequence(genome_file)
	metadata = vcf_generator_ver2.load_metadata(input_template)
	files = os.listdir(input_files_path)

	for file in files:
		filename = os.path.basename(file)
		if '.tsv' in filename:
			output_file_name = filename.replace('.tsv', '_spliceai.vcf')
			output = open(input_files_path+output_file_name, 'w')
			for line in metadata:
				output.write(line)
			output.write(sep.join(headers))
			output.write('\n')

			vcf_generator_ver2.create_gene_chromosome_map(genesfile)
			variants_list = vcf_generator_ver2.process_input_file(input_files_path+file)

			print("Total variants in the file:", len(variants_list))

			(subs, others) = vcf_generator_ver2.segregate_variants()

			print("#Substitution variants:", len(subs))
			print("#INDEL variants:", len(others))

			create_substitution_entries(output, subs)
			create_non_substitution_entries(output, others)

		elif '.vcf' in filename:
			output_file_name = filename.replace('.vcf', '_spliceai.vcf')
			output = open(input_files_path+output_file_name, 'w')
			for line in metadata:
				output.write(line)
			output.write(sep.join(headers))
			output.write('\n')

			print("INPUT FILE:", filename)
	        process_entries(input_files_path+file, output)
	        print("OUTPUT FILE:", output_file_name)

	print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
