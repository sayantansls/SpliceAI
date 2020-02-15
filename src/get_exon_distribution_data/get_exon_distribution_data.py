# author -- sayantan 
# The purpose of this code is as follows: 1 -- This code parses through and obtains the genes in 
# the splice_variants.tsv file 2 -- Then for those genes it parses through the file transcripts.tsv and obtains the 
# strong exons (present in all transcripts) for those corresponding genes 

from __future__ import print_function
import utils
import time as tm
from collections import OrderedDict
from collections import defaultdict


# The headers in splice_variants.tsv file are as follows:
# 1 -- Source	
# 2 -- Case ID	
# 3 -- UK1	
# 4 -- UK2	
# 5 -- Panel	
# 6 -- Test	
# 7 -- Time and Date	
# 8 -- ClinicalManifestation	
# 9 -- Gender	
# 10 -- Ethnicity	
# 11 -- Age	
# 12 -- ReportWiseGeneListInfo
# 13 -- UK3	
# 14 -- UK4	
# 15 -- UK5	
# 16 -- Chromosome	
# 17 -- Gene	
# 18 -- pHGVS	
# 19 -- Genomic HGVS	
# 20 -- Transcript	
# 21 -- VariantType	
# 22 -- AlleleFrequency	
# 23 -- Global PPDB	
# 24 -- Local PPDB	
# 25 -- Zygosity	
# 26 -- Clinvar IDs	
# 27 -- Variant Label Reason	
# 28 -- UK6	
# 29 -- rsID	
# 30 -- EVS	
# 31 -- ExAc	
# 32 -- dbSNP	
# 33 -- 1000 Genomes	
# 34 -- HGMD ID	
# 35 -- Bioinfo Summary	
# 36 -- Literature Summary

# The headers in the transcripts.tsv file are as follows:
# 1 -- Strand_gene_id	
# 2 -- Transcript_name	
# 3 -- Transcript_start	
# 4 -- Transcript_end	
# 5 -- CDS_start	
# 6 -- CDS_end	
# 7 -- Exon_count	
# 8 -- Exon_start	
# 9 -- Exon_end	
# 10 -- Str_validity	
# 11 -- Genomic_override	
# 12 -- Non_standard_codon	
# 13 -- Version	
# 14 -- is_deleted	
# 15 -- invalid_reason	
# 16 -- positional_override	
# 17 -- Comments

def create_gene_transcript_to_exons_map(transcripts_file):
    genes_exons = OrderedDict()

    for line in utils.records_iterator(transcripts_file):
        key = line['Strand_gene_id'] + "|" + line['Transcript_name']
        if key not in genes_exons:
            genes_exons[key] = dict()
        exon_start = line['Exon_start'].split(',')
        exon_end = line['Exon_end'].split(',')

        if '' in exon_start:
            exon_start.remove('')
        if '' in exon_end:
            exon_end.remove('')

        genes_exons[key] = {'exon_start': exon_start,
                            'exon_end': exon_end}

    genes_exons_limits = OrderedDict()

    for key in genes_exons.keys():
        exon_limits = list()
        for i in range(len(genes_exons[key]['exon_start'])):
            exon_limit = (genes_exons[key]['exon_start'][i], genes_exons[key]['exon_end'][i])
            exon_limits.append(exon_limit)
        genes_exons_limits[key] = exon_limits

    global exons_dict
    exons_dict = defaultdict(list)

    for key, value in genes_exons_limits.items():
        new_key = key.split('|')[0]
        exons_dict[new_key].extend(value)


def main(transcripts_file):
    print("Start of code:", tm.ctime(tm.time()))

    create_gene_transcript_to_exons_map(transcripts_file)

    print(exons_dict)
    print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
