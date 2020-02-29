"""
@author: sayantan
The purpose of this code is as follows:
1 -- This code consumes two files - positive.tsv and negative.tsv and gets some statistics
2 -- It obtains the splice variants from the positive.tsv file and displays splice statistics
"""

from __future__ import print_function
import time as tm
import copy
import utils
import re
import os

"""
The headers in the positive.tsv file are as follows:
1 -- CaseID
2 -- UK1
3 -- UK2
4 -- Panel
5 -- Test
6 -- Date/Time
7 -- ClinicalManifestations
8 -- Gender
9 -- Ethnicity
10 -- Age
11 -- ReportWiseGeneListInfo
12 -- UK3
13 -- UK4
14 -- UK5
15 -- Chromosome
16 -- Gene
17 -- p./c.HGVS
18 -- GenomicHGVS
19 -- Transcript
20 -- VariantType
21 -- AlleleFrequency
22 -- GlobalPPDB
23 -- LocalPPDB
24 -- Zygosity
25 -- ClinvarIDs
26 -- VariantLabelReason
27 -- UK6
28 -- rsID
29 -- EVS
30 -- ExAc
31 -- dbSNP
32 -- 1000genomes
33 -- HGMD ID
34 -- Bioinfo Summary
35 -- Literature Summary
"""

"""
The headers in the negative.tsv file are as follows:
1 -- CaseID
2 -- UK1
3 -- UK2
4 -- Panel
5 -- Test
6 -- Date/Time
7 -- ClinicalManifestations
8 -- Gender
9 -- Ethnicity
10 -- Age
11 -- ReportWiseGeneListInfo
12 -- UK3
13 -- UK4
14 -- UK5
"""

def create_age_distribution(age_list):
    below_ten, ten_to_twenty, above_twenty = [0, 0, 0]

    for case_age in age_list:
        if 'Not Available' not in case_age:
            str_age = case_age.split('-')[-1].strip()
            try:
                int_age = int(''.join(re.findall(r'[0-9]+', str_age)))
            except:
                pass

            if 'years' in str_age or 'year' in str_age:
                age = int_age
            elif 'months' in str_age or 'month' in str_age:
                age = int_age / 12
            elif 'days' in str_age or 'day' in str_age:
                age = int_age / 365
            else:
                pass

            if age < 10:
                below_ten += 1
            elif 10 <= age <= 20:
                ten_to_twenty += 1
            elif age > 20:
                above_twenty += 1

    age_stats = {'below 10 years': below_ten,
                 'between 10 to 20 years': ten_to_twenty,
                 'above 20 years': above_twenty}
    return age_stats


def create_gender_distribution(gender_list):
    male_variants, female_variants = [0, 0]

    for case_gender in gender_list:
        gender = case_gender.split('-')[-1].strip()

        if gender == 'Male':
            male_variants += 1
        elif gender == 'Female':
            female_variants += 1

    gender_stats = {'Male': male_variants,
                    'Female': female_variants}
    return gender_stats


def create_variant_label_distribution(variant_label_list):
    pathogenic, likely_pathogenic, likely_benign, vus, vusd = [0, 0, 0, 0, 0]

    for label in variant_label_list:
        if label == 'Pathogenic':
            pathogenic += 1
        elif label == 'LikelyPathogenic':
            likely_pathogenic += 1
        elif label == 'LikelyBenign':
            likely_benign += 1
        elif label == 'VUS':
            vus += 1
        elif label == 'VUSD':
            vusd += 1
        else:
            print("undefined label format", label)

    variant_label_stats = {'Pathogenic': pathogenic,
                           'Likely Pathogenic': likely_pathogenic,
                           'Likely Benign': likely_benign,
                           'VUS': vus,
                           'VUSD': vusd}
    return variant_label_stats


test_cohorts, genelist_cohorts = [set(), set()]


def create_cohorts_list(test_file, genelist_file):
    f_test = open(test_file, 'r')
    for line in f_test.readlines():
        test_cohorts.add(line.strip('\n'))

    f_genelist = open(genelist_file, 'r')
    for line in f_genelist.readlines():
        genelist_cohorts.add(line.split('-')[-1].strip('\n'))

    return test_cohorts, genelist_cohorts


def is_test_cohort(test):
    if test in test_cohorts:
        return 1
    else:
        return 0


def is_genelist_cohort(genelist):
    if genelist in genelist_cohorts:
        return 1
    else:
        return 0


def is_pathogenic(variant_label):
    pathogenic_labels = ['Pathogenic', 'LikelyPathogenic', 'VUSD']
    if variant_label in pathogenic_labels:
        return 1
    else:
        return 0


def create_reported_cases_statistics(inputfile):
    cases, clinical_exome_cases, clinical_exome_pathogenic_cases = [0, 0, 0]
    clinical_exome_cohort_cases, clinical_exome_cohort_pathogenic_cases = [0, 0]
    gender_list, age_list, variant_label_list = [list(), list(), list()]

    for case in utils.records_iterator(inputfile):
        cases += 1

        filename = os.path.basename(inputfile)

        if case['Panel'] == 'Trusight one clinical exome(4805 genes)':
            clinical_exome_cases += 1

            if filename != 'negative.tsv':
                if is_pathogenic(case['VariantLabelReason']):
                    clinical_exome_pathogenic_cases += 1

            if is_test_cohort(case['Test'].strip()) == 1 or is_genelist_cohort(case['ReportWiseGeneListInfo'].strip()) == 1:
                clinical_exome_cohort_cases += 1
                #if filename == 'negative.tsv':
                 #   print(case['CaseID'] + "\t" + case['UK1'] + "\t" + case['UK2'] + "\t" + case['Panel'] + "\t" + case['Test'] + "\t" + case['ClinicalManifestations'])

                if filename == 'positive.tsv':
                    if is_pathogenic(case['VariantLabelReason']):
                        clinical_exome_cohort_pathogenic_cases += 1

                        age_list.append(case['Age'])
                        gender_list.append(case['Gender'])
                        variant_label_list.append(case['VariantLabelReason'])
                else:
                    age_list.append(case['Age'])
                    gender_list.append(case['Gender'])

    gender_stats = create_gender_distribution(gender_list)
    age_stats = create_age_distribution(age_list)
    variant_label_stats = create_variant_label_distribution(variant_label_list)

    return cases, clinical_exome_cases, clinical_exome_pathogenic_cases, clinical_exome_cohort_cases, clinical_exome_cohort_pathogenic_cases, \
           gender_stats, age_stats, variant_label_stats


def create_printable_data(file, cases, clinical_exome_cases, clinical_exome_pathogenic_cases, clinical_exome_cohort_cases,
                          clinical_exome_cohort_pathogenic_cases, gender_stats, age_stats, variant_label_stats):
    filename = os.path.basename(file)

    print("#Cases in the file:", cases)
    print("#Cases in the clinical exome panel:", clinical_exome_cases)

    if filename == 'positive.tsv':
        print("#Cases in the clinical exome panel (LP + P + VUSD):", clinical_exome_pathogenic_cases)
    print("#Cases in the clinical exome + neurodevelopmental cohort:", clinical_exome_cohort_cases)

    if filename == 'positive.tsv':
        print("#Cases in the clinical exome + neurodevelopmental cohort (LP + P + VUSD):", clinical_exome_cohort_pathogenic_cases)
    print("Gender Stats:", gender_stats)
    print("Age Stats:", age_stats)

    if filename == 'positive.tsv':
        print("Variant Label Stats:", variant_label_stats)


def main(positive_file, negative_file, test_file, genelist_file):
    print("Start of code:", tm.ctime(tm.time()))

    create_cohorts_list(test_file, genelist_file)

    (positive_cases, positive_clinical_exome_cases, positive_clinical_exome_pathogenic_cases,
     positive_clinical_exome_cohort_cases, positive_clinical_exome_cohort_pathogenic_cases, positive_gender_stats, positive_age_stats,
     positive_variant_label_stats) = create_reported_cases_statistics(positive_file)

    (negative_cases, negative_clinical_exome_cases, negative_clinical_exome_pathogenic,
     negative_cohort_cases, negative_cohort_case_pathogenic, negative_gender_stats, negative_age_stats,
     negative_variant_label_stats) = create_reported_cases_statistics(negative_file)

    print("Positive File statistics:")
    create_printable_data(positive_file, positive_cases, positive_clinical_exome_cases, positive_clinical_exome_pathogenic_cases,
                          positive_clinical_exome_cohort_cases, positive_clinical_exome_cohort_pathogenic_cases, positive_gender_stats,
                          positive_age_stats, positive_variant_label_stats)

    print("Negative File statistics:")
    create_printable_data(negative_file, negative_cases, negative_clinical_exome_cases, negative_clinical_exome_pathogenic,
                          negative_cohort_cases, negative_cohort_case_pathogenic, negative_gender_stats,
                          negative_age_stats, negative_variant_label_stats)

    print("End of code:", tm.ctime(tm.time()))


if __name__ == '__main__':
    import sys
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
