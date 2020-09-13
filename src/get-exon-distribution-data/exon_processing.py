"""
@author -- sayantan 
The purpose of this code is as follows: 1 -- This code parses through and obtains the genes in 
the splice_variants.tsv file 2 -- Then for those genes it parses through the file transcripts.tsv and obtains the 
strong exons (present in all transcripts) for those corresponding genes 
"""

from __future__ import print_function
import utils
import time as tm
import copy
from collections import OrderedDict, defaultdict

"""
The headers in the transcripts.tsv file are as follows:
1 -- Strand_gene_id	
2 -- Transcript_name	
3 -- Transcript_start	
4 -- Transcript_end	
5 -- CDS_start	
6 -- CDS_end	
7 -- Exon_count	
8 -- Exon_start	
9 -- Exon_end	
10 -- Str_validity	
11 -- Genomic_override	
12 -- Non_standard_codon	
13 -- Version	
14 -- is_deleted	
15 -- invalid_reason	
16 -- positional_override	
17 -- Comments
"""
sep = '\t'

def create_gene_transcripts_map(transcripts_file):
    """
    gene_transcripts creates a dictionary in the following format:
    {'EGFR-7-1':['NM_201283', 'NM_005228', 'NM_201284', 'NM_201282']}
    """
    gene_transcripts = dict()

    for gene in utils.records_iterator(transcripts_file):
        if gene['Strand_gene_id'] not in gene_transcripts:
            gene_transcripts[gene['Strand_gene_id']] = list()
        gene_transcripts[gene['Strand_gene_id']].append(gene['Transcript_name'])

    """
    gene_transcripts_count creates a dictionary of gene id and number of transcripts in the following format:
    {'EGFR-7-1': 4}
    """
    global gene_transcripts_count
    gene_transcripts_count = dict()

    for key in gene_transcripts.keys():
        gene_transcripts_count[key] = len(gene_transcripts[key])

def create_file_genes_list_with_three_or_more_transcripts():
    gene_list = list()

    for key in gene_transcripts_count.keys():
        if gene_transcripts_count[key] >= 3:
            gene_lista.append(key)

def create_gene_transcript_to_exons_map(transcripts_file):
    """
    gene_exons creates a dictionary with the 'gene_id|transcript' as the key and 
    the exon_start and exon_end for all the corresponding transcript as the value as follows:
    {'EGFR-7-1|NM_005228': {'exon_end': ['55087058','55210130','55211181','55214433',.....],
                            'exon_start': ['55086725','55209979','55210998','55214299',....]}}
    """
    genes_exons = OrderedDict()

    for line in utils.records_iterator(transcripts_file):
        key = line['Strand_gene_id'] + "|" + line['Transcript_name']
        if key not in genes_exons:
            genes_exons[key] = dict()
        exon_start = line['Exon_start'].split(',')
        exon_end = line['Exon_end'].split(',')

        exon_start = list(filter(None, exon_start))
        exon_end = list(filter(None, exon_end))

        genes_exons[key] = {'exon_start': exon_start,
                            'exon_end': exon_end}

    """
    genes_exons_limits creates a dictionary with the 'gene_id|transcript' as the key and
    combines each exon_start and exon_end as a tuple and the list of tuples is the value as follows:
    {'EGFR-7-1|NM_005228':[('55086725', '55087058'),
                           ('55209979', '55210130'),
                           ('55210998', '55211181'),
                           ('55214299', '55214433'),
                           ('55218987', '55219055'),...]}
    """
    genes_exons_limits = OrderedDict()

    for key in genes_exons.keys():
        exon_limits = list()
        for i in range(len(genes_exons[key]['exon_start'])):
            exon_limit = (genes_exons[key]['exon_start'][i], genes_exons[key]['exon_end'][i])
            exon_limits.append(exon_limit)
        genes_exons_limits[key] = exon_limits

    """
    exons_dict combines the tuples for all the respective transcripts for a particular gene_id
    {'EGFR-7-1':[('55086725', '55087058'),
                 ('55209979', '55210130'),
                 ('55210998', '55211181'),
                 ('55214299', '55214433'),
                 ('55218987', '55219055'),...]}
    """
    exons_dict = defaultdict(list)

    for key, value in genes_exons_limits.items():
        new_key = key.split('|')[0]
        exons_dict[new_key].extend(value)

    """
    exons_counts_dict creates a dictionary for a particular gene_id, the frequency of each tuple 
    in the list of tuples (each tuple is a combination of exon_start and exon_end)
    {'EGFR-7-1': {[('55086725', '55087058'): 4,
                   ('55209979', '55210130'): 4,
                   ('55210998', '55211181'): 4,
                   ('55214299', '55214433'): 4,
                   ('55224452', '55224525'): 3,
                   ('55224452', '55224644'): 1,
                   ('55225356', '55225446'): 3,...]}
    """
    global exons_counts_dict
    exons_counts_dict = OrderedDict()
    for key in exons_dict.keys():
        freq = dict()
        for item in exons_dict[key]:
            if item in freq:
                freq[item] += 1
            else:
                freq[item] = 1
        exons_counts_dict[key] = freq

def get_all_strong_exons():
    global all_strong_exons
    all_strong_exons = copy.deepcopy(exons_counts_dict)
    for key in all_strong_exons.keys():
        for sub_key in all_strong_exons[key].keys():
            if all_strong_exons[key][sub_key] != gene_transcripts_count[key]:
                del all_strong_exons[key][sub_key]

def get_all_rare_exons():
    global all_rare_exons
    all_rare_exons = copy.deepcopy(exons_counts_dict)
    for key in all_rare_exons.keys():
        for sub_key in all_rare_exons[key].keys():
            if all_rare_exons[key][sub_key] == gene_transcripts_count[key]:
                del all_rare_exons[key][sub_key]

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

def create_gene_details_map(genes_file):
    global gene_details_dict 
    gene_details_dict = dict()
    for line in utils.records_iterator(genes_file):
        gene_details_dict[line['Symbol']] = {'gene_id': line['Strand_gene_id'],
                                             'chromosome': line['ChrName'],
                                             'strand': line['Strand']}

HEADERS = ['Gene',
           'Chromosome',
           'Strand',
           'Exon_start',
           'Exon_end']

ENTRY_T = {'Gene': '',
           'Chromosome': '',
           'Strand': '',
           'Exon_start': '',
           'Exon_end': ''}

"""
The headers in splice_variants.tsv file are as follows:
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

def get_file_genes_list_from_file(splice_variants_file):
    global file_genes_list 
    file_genes_list = set()
    for entry in utils.records_iterator(splice_variants_file):
        file_genes_list.add(entry['Gene'])

def process_strong_exon_entries(output_strong_exons):
    for gene_entry in file_genes_list:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene_id = gene_details_dict[gene_entry]['gene_id']

        if gene_transcripts_count[gene_id] >= 3:
            ENTRY['Gene'] = gene_entry
            ENTRY['Chromosome'] = gene_details_dict[gene_entry]['chromosome']
            ENTRY['Strand'] = gene_details_dict[gene_entry]['strand']
            exon_ranges = all_strong_exons[gene_id].keys()
            for item in exon_ranges:
                (ENTRY['Exon_start'],ENTRY['Exon_end']) = item
                field_values = [ENTRY[i] for i in HEADERS]
                output_strong_exons.write(sep.join(field_values))
                output_strong_exons.write('\n')
        else:
            continue

def process_rare_exon_entries(output_rare_exons):
    for gene_entry in file_genes_list:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene_id = gene_details_dict[gene_entry]['gene_id']
        ENTRY['Gene'] = gene_entry
        ENTRY['Chromosome'] = gene_details_dict[gene_entry]['chromosome']
        ENTRY['Strand'] = gene_details_dict[gene_entry]['strand']
        exon_ranges = all_rare_exons[gene_id].keys()
        for item in exon_ranges:
            (ENTRY['Exon_start'],ENTRY['Exon_end']) = item

            field_values = [ENTRY[i] for i in HEADERS]
            output_rare_exons.write(sep.join(field_values))
            output_rare_exons.write('\n')

def main(transcripts_file, splice_variants_file, genes_file, output_strong_exons_file, output_rare_exons_file):
    print("Start of code:", tm.ctime(tm.time()))

    create_gene_transcripts_map(transcripts_file)
    create_gene_transcript_to_exons_map(transcripts_file)
    create_gene_details_map(genes_file)
    get_all_strong_exons()
    get_all_rare_exons()
    get_file_genes_list_from_file(splice_variants_file)

    output_strong_exons = open(output_strong_exons_file, 'w')
    output_strong_exons.write(sep.join(HEADERS))
    output_strong_exons.write('\n')
    process_strong_exon_entries(output_strong_exons)

    output_rare_exons = open(output_rare_exons_file, 'w')
    output_rare_exons.write(sep.join(HEADERS))
    output_rare_exons.write('\n')
    process_rare_exon_entries(output_rare_exons)

    print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
