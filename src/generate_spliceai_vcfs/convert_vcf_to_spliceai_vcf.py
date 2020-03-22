"""
@author : sayantan
This code converts StrandOmics VCF to Splice AI input VCF format
"""
import utils
import sys, os, copy
import time as tm

"""
The headers in StrandOmics VCF format are as follows:
1 -- #CHROM
2 -- POS
3 -- ID
4 -- REF
5 -- ALT
6 -- QUAL
7 -- FILTER
8 -- INFO
9 -- FORMAT
10 -- DETAILS
"""
HEADERS = ['#CHROM',
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
"""
The headers in the SpliceAI VCF file are as follows:
1 -- #CHROM
2 -- POS
3 -- ID
4 -- REF
5 -- ALT
6 -- QUAL
7 -- FILTER
8 -- INFO
"""
def process_entries(vcf_file, output):
    for variant in utils.records_iterator(vcf_file):
        ENTRY = copy.deepcopy(ENTRY_T)

        ENTRY['#CHROM'] = variant['CHROM'].replace('chr','')
        ENTRY['POS'] = variant['POS']
        ENTRY['REF'] = variant['REF']
        ENTRY['ALT'] = variant['ALT']

        field_values = [str(ENTRY[i]) for i in HEADERS]
        output.write(sep.join(field_values))
        output.write('\n')

def main(vcf_folder):
    print("Start of code:", tm.ctime(tm.time()))
    files = os.listdir(vcf_folder)
    for file in files:
        filename = os.path.basename(file)
        output_file_name = filename.replace('.vcf','_spliceai.vcf')
        output = open(vcf_folder+output_file_name, 'w')
        output.write(sep.join(HEADERS))
        output.write('\n')
        print("INPUT FILE:", filename)
        process_entries(vcf_folder+file, output)
        print("OUTPUT FILE:", output_file_name)

    print("End of code:", tm.ctime(tm.time()))

if __name__ == '__main__':
    main(sys.argv[1])
