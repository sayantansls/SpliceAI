"""
@author - sayantan
The purpose of this code is as follows:
1 -- Consume a tsv file "SpliceVariants.tsv" and get the genomicHGVS
2 -- Create a sliding window range for the particular splice prediction software
  a. SplicePort
  b. NNSplice
  c. ASSP
  d. MaxEntScan
  e. Human Splicing Finder (HSF)
Parse the human genomic sequence and get the sequences for each of the sliding ranges
Create a tsv file with the particular variant, details, sequences fdr each predictor
"""

from __future__ import print_function
import time as tm
import copy
import re
import utils
import twobitreader

sep = '\t'

predictor_keywords = ['NNSPLICE', 'ASSP', 'MaxEntScan', 'SplicePort', 'Human Splicing Finder']

HEADER = ['caseID',
          'panel',
          'test',
          'clinical_manifestation',
          'gender',
          'ethnicity',
          'age',
          'report_wise_gene_list_info',
          'chromosome',
          'gene',
          'strand',
          'c./p.HGVS',
          'genomicHGVS',
          'transcript',
          'predictors',
          'RefSeq_MES_donor',
          'VarSeq_MES_donor',
          'RefSeq_NNSplice_donor',
          'VarSeq_NNSplice_donor',
          'RefSeq_SplicePort_donor',
          'VarSeq_SplicePort_donor',
          'RefSeq_HSF_donor',
          'VarSeq_HSF_donor',
          'RefSeq_ASSP_donor',
          'VarSeq_ASSP_donor',
          'RefSeq_MES_acceptor',
          'VarSeq_MES_acceptor',
          'RefSeq_NNSplice_acceptor',
          'VarSeq_NNSplice_acceptor',
          'RefSeq_SplicePort_acceptor',
          'VarSeq_SplicePort_acceptor',
          'RefSeq_HSF_acceptor',
          'VarSeq_HSF_acceptor',
          'RefSeq_ASSP_acceptor',
          'VarSeq_ASSP_acceptor',
          'variant_type',
          'bioinfo_summary',
          'literature_summary']

ENTRY_T = {'caseID': '',
           'panel': '',
           'test': '',
           'clinical_manifestation': '',
           'gender': '',
           'ethnicity': '',
           'age': '',
           'report_wise_gene_list_info': '',
           'chromosome': '',
           'gene': '',
           'strand': '',
           'c./p.HGVS': '',
           'genomicHGVS': '',
           'transcript': '',
           'predictors': '',
           'RefSeq_MES_donor': '',
           'VarSeq_MES_donor': '',
           'RefSeq_NNSplice_donor': '',
           'VarSeq_NNSplice_donor': '',
           'RefSeq_SplicePort_donor': '',
           'VarSeq_SplicePort_donor': '',
           'RefSeq_HSF_donor': '',
           'VarSeq_HSF_donor': '',
           'RefSeq_ASSP_donor': '',
           'VarSeq_ASSP_donor': '',
           'RefSeq_MES_acceptor': '',
           'VarSeq_MES_acceptor': '',
           'RefSeq_NNSplice_acceptor': '',
           'VarSeq_NNSplice_acceptor': '',
           'RefSeq_SplicePort_acceptor': '',
           'VarSeq_SplicePort_acceptor': '',
           'RefSeq_HSF_acceptor': '',
           'VarSeq_HSF_acceptor': '',
           'RefSeq_ASSP_acceptor': '',
           'VarSeq_ASSP_acceptor': '',
           'variant_type': '',
           'bioinfo_summary': '',
           'literature_summary': ''}

genome = twobitreader.TwoBitFile('/home/sayantan/Desktop/hg19.2bit')


def get_genomic_position(genomicHGVS):
    if ">" in genomicHGVS:  # substitution variants
        position = ''.join(re.findall(r'[0-9]+', genomicHGVS))
    elif "del" in genomicHGVS:  # deletion variants
        positions = re.findall(r'[0-9]+', genomicHGVS)
        if len(positions) == 1:
            position = int(positions[0]) - 1
        else:
            start, end = [positions[0], positions[1]]
            position = int(start) - 1
    else:  # insertion and delins variants
        positions = re.findall(r'[0-9]+', genomicHGVS)
        position = positions[0]
    return position


def get_ref_alt(position):
    try:
        ref, alt = re.findall(r'[A-Z]', position)
        return ref, alt
    except:
        pass

def create_gene_strand_map(genes_file):
  global gene_strand_map
  gene_strand_map = dict()
  for gene in utils.records_iterator(genes_file):
    gene_strand_map[gene['Symbol']] = gene['Strand']

def get_strand_sense(gene):
  return gene_strand_map[gene]

def create_mes_ranges(position, strand):
    if strand == '+':
        donor_splice = {'start': (position - 2), 'end': (position + 6)}
        acceptor_splice = {'start': (position - 20), 'end': (position + 2)}
    else:
        donor_splice = {'start': (position - 6), 'end': (position + 2)}
        acceptor_splice = {'start': (position - 2), 'end': (position + 20)}

    return donor_splice, acceptor_splice


def create_mes_donor_range(position, strand):
    if strand == '+':
        donor_splice = {'start': (position - 2), 'end': (position + 6)}
    elif strand == '-':
        donor_splice = {'start': (position - 6), 'end': (position + 2)}
    else:
        print(position, "Wrong strand sense")

    return donor_splice


def create_mes_acceptor_range(position, strand):
    if strand == '+':
        acceptor_splice = {'start': (position - 20), 'end': (position + 2)}
    elif strand == '-':
        acceptor_splice = {'start': (position - 2), 'end': (position + 20)}
    else:
        print(position, "Wrong strand sense")

    return acceptor_splice


def create_hsf_ranges(position):
    # Human Splicing Finder (HSF)
    # Donor splice  site = -40 to +40
    # Acceptor splice site = -40 to +40
    donor_splice = {'start': (position - 40), 'end': (position + 41)}
    acceptor_splice = {'start': (position - 40), 'end': (position + 41)}
    return donor_splice, acceptor_splice


def create_nnsplice_ranges(position):
    # NN Splice
    # Donor splice  site = -7 to +8
    # Acceptor splice site = -20 to +20
    # Build handling for negative strand?????
    donor_splice = {'start': (position - 7), 'end': (position + 8)}
    acceptor_splice = {'start': (position - 20), 'end': (position + 20)}
    return donor_splice, acceptor_splice


def create_nnsplice_donor_range(position, strand):
    if strand == '+':
        donor_splice = {'start': (position - 6), 'end': (position + 8)}
    else:
        donor_splice = {'start': (position - 8), 'end': (position + 6)}
    return donor_splice


def create_nnsplice_acceptor_range(position, strand):
    if strand == '+':
        acceptor_splice = {'start': (position - 21), 'end': (position + 19)}
    else:
        acceptor_splice = {'start': (position - 19), 'end': (position + 21)}
    return acceptor_splice


def create_spliceport_ranges(position):
    # Spliceport
    # Donor splice  site = -80 to +80
    # Acceptor splice site = -80 to +80
    donor_splice = {'start': (position - 80), 'end': (position + 81)}
    acceptor_splice = {'start': (position - 80), 'end': (position + 81)}
    return donor_splice, acceptor_splice


def create_assp_ranges(position):
    # ASSP
    # Donor splice site = -70 to +70
    # Acceptor splice site = -70 to +70
    donor_splice = {'start': (position - 70), 'end': (position + 70)}
    acceptor_splice = {'start': (position - 70), 'end': (position + 70)}
    return donor_splice, acceptor_splice


def create_assp_hsf_range(position):
    splice_range = {'start': (position - 70), 'end': (position + 70)}
    return splice_range


def create_reverse_complementary_sequence(sequence):
    complementary = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    c_sequence = str()
    for base in sequence:
        c_sequence += complementary[base]
    rc_sequence = c_sequence[::-1]
    return rc_sequence


def get_original_sequence(predictor, chrom):
    sequence = genome[chrom][predictor['start'] - 1:predictor['end']]

    return sequence


def get_sequence(predictor, base, position, chrom):
    sequence = genome[chrom][predictor['start'] - 1:position] + base + genome[chrom][
                                                                           position+1:predictor['end']]
    return sequence

def get_reference_sequence(predictor, position, chrom):
    sequence = genome[chrom][predictor['start'] -1:position] + genome[chrom][position] + genome[chrom][position+1:predictor['end']]
    return sequence


"""
The headers in the "SpliceVariants.tsv" file are as follows:
1 -- Source
2 -- Case ID
3 -- UK1 (Unknown)
4 -- UK2 (Unknown)
5 -- Panel
6 -- Test
7 -- Time and Date
8 -- ClinicalManifestation
9 -- Gender
10 -- Ethnicity
11 -- Age
12 -- ReportWiseGeneListInfo
13 -- UK3 (Unknown)
14 -- UK4 (Unknown)
15 -- UK5 (Unknown)
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
28 -- UK6 (Unknown)
29 -- rsID
30 -- EVS
31 -- ExAc
32 -- dbSNP
33 -- 1000 Genomes
34 -- HGMD ID
35 -- Bioinfo Summary
36 -- Literature Summary
"""


def create_splice_matrix(inputfile, outfile):
    for variant in utils.records_iterator(inputfile):
        ENTRY = copy.deepcopy(ENTRY_T)

        ENTRY['caseID'] = variant['Case ID']
        ENTRY['panel'] = variant['Panel']
        ENTRY['test'] = variant['Test']
        ENTRY['clinical_manifestation'] = variant['ClinicalManifestation']
        ENTRY['gender'] = variant['Gender'].split('-')[-1]
        ENTRY['ethnicity'] = variant['Ethnicity'].split('-')[-1]
        ENTRY['age'] = variant['Age'].split('-')[-1]
        ENTRY['report_wise_gene_list_info'] = variant['ReportWiseGeneListInfo']
        ENTRY['chromosome'] = variant['Chromosome']
        ENTRY['gene'] = variant['Gene']
        strand = get_strand_sense(variant['Gene'])
        ENTRY['strand'] = strand
        ENTRY['c./p.HGVS'] = variant['pHGVS']
        ENTRY['genomicHGVS'] = variant['Genomic HGVS']
        ENTRY['transcript'] = variant['Transcript']
        ENTRY['variant_type'] = variant['VariantType']
        bioinfo = variant['Bioinfo Summary']
        ENTRY['bioinfo_summary'] = bioinfo

        predictors = ''

        for keyword in predictor_keywords:
            if keyword in bioinfo:
                predictors += keyword + ','

        ENTRY['predictors'] = predictors[:-1]

        ENTRY['literature_summary'] = variant['Literature Summary']

        position = int(get_genomic_position(variant['Genomic HGVS']))

        try:
            reference, alternate = get_ref_alt(variant['Genomic HGVS'])
        except:
            pass

        (mes_donor, mes_acceptor) = create_mes_ranges(position, strand)
        (nnsplice_donor, nnsplice_acceptor) = create_nnsplice_ranges(position)
        (spliceport_donor, spliceport_acceptor) = create_spliceport_ranges(position)
        (hsf_donor, hsf_acceptor) = create_hsf_ranges(position)
        (assp_donor, assp_acceptor) = create_assp_ranges(position)

        chrom = variant['Chromosome']

        refSeq_MES_donor = get_sequence(mes_donor, reference, position, chrom)
        refSeq_NNSplice_donor = get_sequence(nnsplice_donor, reference, position, chrom)
        refSeq_SplicePort_donor = get_sequence(spliceport_donor, reference, position, chrom)
        refSeq_HSF_donor = get_sequence(hsf_donor, reference, position, chrom)
        refSeq_ASSP_donor = get_sequence(assp_donor, reference, position, chrom)

        refSeq_MES_acceptor = get_sequence(mes_acceptor, reference, position, chrom)
        refSeq_NNSplice_acceptor = get_sequence(nnsplice_acceptor, reference, position, chrom)
        refSeq_SplicePort_acceptor = get_sequence(spliceport_acceptor, reference, position, chrom)
        refSeq_HSF_acceptor = get_sequence(hsf_acceptor, reference, position, chrom)
        refSeq_ASSP_acceptor = get_sequence(assp_acceptor, reference, position, chrom)

        varSeq_MES_donor = get_sequence(mes_donor, alternate, position, chrom)
        varSeq_NNSplice_donor = get_sequence(nnsplice_donor, alternate, position, chrom)
        varSeq_SplicePort_donor = get_sequence(spliceport_donor, alternate, position, chrom)
        varSeq_HSF_donor = get_sequence(hsf_donor, alternate, position, chrom)
        varSeq_ASSP_donor = get_sequence(assp_donor, alternate, position, chrom)

        varSeq_MES_acceptor = get_sequence(mes_acceptor, alternate, position, chrom)
        varSeq_NNSplice_acceptor = get_sequence(nnsplice_acceptor, alternate, position, chrom)
        varSeq_SplicePort_acceptor = get_sequence(spliceport_acceptor, alternate, position, chrom)
        varSeq_HSF_acceptor = get_sequence(hsf_acceptor, alternate, position, chrom)
        varSeq_ASSP_acceptor = get_sequence(assp_acceptor, alternate, position, chrom)

        if strand == '+':
            ENTRY['RefSeq_MES_donor'] = refSeq_MES_donor
            ENTRY['RefSeq_NNSplice_donor'] = refSeq_NNSplice_donor
            ENTRY['RefSeq_SplicePort_donor'] = refSeq_SplicePort_donor
            ENTRY['RefSeq_HSF_donor'] = refSeq_HSF_donor
            ENTRY['RefSeq_ASSP_donor'] = refSeq_ASSP_donor

            ENTRY['RefSeq_MES_acceptor'] = refSeq_MES_acceptor
            ENTRY['RefSeq_NNSplice_acceptor'] = refSeq_NNSplice_acceptor
            ENTRY['RefSeq_SplicePort_acceptor'] = refSeq_SplicePort_acceptor
            ENTRY['RefSeq_HSF_acceptor'] = refSeq_HSF_acceptor
            ENTRY['RefSeq_ASSP_acceptor'] = refSeq_ASSP_acceptor

            ENTRY['VarSeq_MES_donor'] = varSeq_MES_donor
            ENTRY['VarSeq_NNSplice_donor'] = varSeq_NNSplice_donor
            ENTRY['VarSeq_SplicePort_donor'] = varSeq_SplicePort_donor
            ENTRY['VarSeq_HSF_donor'] = varSeq_HSF_donor
            ENTRY['VarSeq_ASSP_donor'] = varSeq_ASSP_donor

            ENTRY['VarSeq_MES_acceptor'] = varSeq_MES_acceptor
            ENTRY['VarSeq_NNSplice_acceptor'] = varSeq_NNSplice_acceptor
            ENTRY['VarSeq_SplicePort_acceptor'] = varSeq_SplicePort_acceptor
            ENTRY['VarSeq_HSF_acceptor'] = varSeq_HSF_acceptor
            ENTRY['VarSeq_ASSP_acceptor'] = varSeq_ASSP_acceptor

        else:
            ENTRY['RefSeq_MES_donor'] = create_reverse_complementary_sequence(refSeq_MES_donor)
            ENTRY['RefSeq_NNSplice_donor'] = create_reverse_complementary_sequence(refSeq_NNSplice_donor)
            ENTRY['RefSeq_SplicePort_donor'] = create_reverse_complementary_sequence(refSeq_SplicePort_donor)
            ENTRY['RefSeq_HSF_donor'] = create_reverse_complementary_sequence(refSeq_HSF_donor)
            ENTRY['RefSeq_ASSP_donor'] = create_reverse_complementary_sequence(refSeq_ASSP_donor)

            ENTRY['RefSeq_MES_acceptor'] = create_reverse_complementary_sequence(refSeq_MES_acceptor)
            ENTRY['RefSeq_NNSplice_acceptor'] = create_reverse_complementary_sequence(refSeq_NNSplice_acceptor)
            ENTRY['RefSeq_SplicePort_acceptor'] = create_reverse_complementary_sequence(refSeq_SplicePort_acceptor)
            ENTRY['RefSeq_HSF_acceptor'] = create_reverse_complementary_sequence(refSeq_HSF_acceptor)
            ENTRY['RefSeq_ASSP_acceptor'] = create_reverse_complementary_sequence(refSeq_ASSP_acceptor)

            ENTRY['VarSeq_MES_donor'] = create_reverse_complementary_sequence(varSeq_MES_donor)
            ENTRY['VarSeq_NNSplice_donor'] = create_reverse_complementary_sequence(varSeq_NNSplice_donor)
            ENTRY['VarSeq_SplicePort_donor'] = create_reverse_complementary_sequence(varSeq_SplicePort_donor)
            ENTRY['VarSeq_HSF_donor'] = create_reverse_complementary_sequence(varSeq_HSF_donor)
            ENTRY['VarSeq_ASSP_donor'] = create_reverse_complementary_sequence(varSeq_ASSP_donor)

            ENTRY['VarSeq_MES_acceptor'] = create_reverse_complementary_sequence(varSeq_MES_acceptor)
            ENTRY['VarSeq_NNSplice_acceptor'] = create_reverse_complementary_sequence(varSeq_NNSplice_acceptor)
            ENTRY['VarSeq_SplicePort_acceptor'] = create_reverse_complementary_sequence(varSeq_SplicePort_acceptor)
            ENTRY['VarSeq_HSF_acceptor'] = create_reverse_complementary_sequence(varSeq_HSF_acceptor)
            ENTRY['VarSeq_ASSP_acceptor'] = create_reverse_complementary_sequence(varSeq_ASSP_acceptor)

        field_values = [ENTRY[i] for i in HEADER]
        outfile.write(sep.join(field_values))
        outfile.write('\n')


def main(input_file, output_file, genes_file):
    print("Start of code:", tm.ctime(tm.time()))
    outfile = open(output_file, 'w')
    outfile.write(sep.join(HEADER))
    outfile.write('\n')

    create_gene_strand_map(genes_file)
    create_splice_matrix(input_file, outfile)
    print("End of code:", tm.ctime(tm.time()))


if __name__ == "__main__":
    import sys

    main(sys.argv[1], sys.argv[2], sys.argv[3])
