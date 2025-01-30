###
# New script based on some code from filterClinicalVariantsSV_v0.1.wdl, 
# tested in test_sv_annotations.ipynb. 
# Created 1/30/2025.

## CHANGE LOG:
'''
1/30/2025:
- added dominant_gt and recessive_gt annotations
'''
###

import datetime
import pandas as pd
import hail as hl
import numpy as np
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='vcf_file', help='Input VCF file')
parser.add_argument('--ped', dest='ped_uri', help='Input ped file')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--build', dest='build', help='Genome build')

args = parser.parse_args()

sv_vcf = args.vcf_file
ped_uri = args.ped_uri
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
genome_build = args.build

def get_transmission(phased_tm):
    phased_tm = phased_tm.annotate_entries(transmission=hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|0'), 'uninherited',
            hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|1'), 'inherited_from_mother',
                        hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|0'), 'inherited_from_father',
                                hl.or_missing(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|1'), 'inherited_from_both'))))
    )
    return phased_tm

def filter_and_annotate_tm(tm, variant_category, filter_type='trio'):
    '''
    Filter entries by presence in trio --> annotate by variant_category
    filter_type can be 'trio' or 'proband'
    '''
    if filter_type=='trio':
        tm = tm.filter_entries((tm.proband_entry.GT.is_non_ref()) | 
                                    (tm.mother_entry.GT.is_non_ref()) |
                                    (tm.father_entry.GT.is_non_ref()))
    if filter_type=='proband':
        tm = tm.filter_entries(tm.proband_entry.GT.is_non_ref())
    tm = tm.annotate_rows(variant_category=variant_category)
    return tm

sv_mt = hl.import_vcf(sv_vcf, force_bgz=sv_vcf.split('.')[-1] in ['gz', 'bgz'], 
    reference_genome=genome_build, array_elements_required=False, call_fields=[])
header = hl.get_vcf_metadata(sv_vcf)

# Annotate affected status/phenotype from pedigree
ped_ht = hl.import_table(ped_uri, delimiter='\t').key_by('sample_id')
sv_mt = sv_mt.annotate_cols(phenotype=ped_ht[sv_mt.s].phenotype)

# Get cohort unaffected/affected het and homvar counts
sv_mt = sv_mt.annotate_rows(**{
    "n_het_unaffected": hl.agg.filter(sv_mt.phenotype=='1', hl.agg.sum(sv_mt.GT.is_het())),
    "n_hom_var_unaffected": hl.agg.filter(sv_mt.phenotype=='1', hl.agg.sum(sv_mt.GT.is_hom_var())),
    "n_het_affected": hl.agg.filter(sv_mt.phenotype=='2', hl.agg.sum(sv_mt.GT.is_het())),
    "n_hom_var_affected": hl.agg.filter(sv_mt.phenotype=='2', hl.agg.sum(sv_mt.GT.is_hom_var()))
})

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
cropped_ped_uri = f"{os.path.basename(ped_uri).split('.ped')[0]}_crop.ped"
tmp_ped.to_csv(cropped_ped_uri, sep='\t', index=False)
pedigree = hl.Pedigree.read(cropped_ped_uri, delimiter='\t')

sv_tm = hl.trio_matrix(sv_mt, pedigree, complete_trios=False)
phased_sv_tm = hl.experimental.phase_trio_matrix_by_transmission(sv_tm, call_field='GT', phased_call_field='PBT_GT')

# Mendel errors
all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(sv_mt['GT'], pedigree)
all_errors_sv_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
phased_sv_tm = phased_sv_tm.annotate_entries(mendel_code=all_errors_sv_mt[phased_sv_tm.row_key, phased_sv_tm.col_key].mendel_code)
phased_sv_tm = get_transmission(phased_sv_tm)

# Annotate if complete trio or not
complete_trio_probands = [trio.s for trio in pedigree.complete_trios()]
phased_sv_tm = phased_sv_tm.annotate_cols(complete_trio=hl.array(complete_trio_probands).contains(phased_sv_tm.proband.s))

# Annotate affected status/phenotype from pedigree
phased_sv_tm = phased_sv_tm.annotate_cols(
    proband=phased_sv_tm.proband.annotate(
        phenotype=ped_ht[phased_sv_tm.proband.s].phenotype),
mother=phased_sv_tm.mother.annotate(
        phenotype=ped_ht[phased_sv_tm.mother.s].phenotype),
father=phased_sv_tm.father.annotate(
        phenotype=ped_ht[phased_sv_tm.father.s].phenotype))

affected_cols = ['n_het_unaffected', 'n_hom_var_unaffected', 'n_het_affected', 'n_hom_var_affected']
phased_sv_tm = phased_sv_tm.annotate_rows(**{col: sv_mt.rows()[phased_sv_tm.row_key][col] 
                                             for col in affected_cols})

## Annotate dominant_gt and recessive_gt
# denovo
dom_trio_criteria = ((phased_sv_tm.complete_trio) &  
                     (phased_sv_tm.mendel_code==2))
# het absent in unaff
dom_non_trio_criteria = ((~phased_sv_tm.complete_trio) & 
                        (phased_sv_tm.n_hom_var_unaffected==0) &
                        (phased_sv_tm.n_het_unaffected==0) & 
                        (phased_sv_tm.proband_entry.GT.is_het())
                        )

has_mother = hl.is_defined(phased_sv_tm.mother)
has_father = hl.is_defined(phased_sv_tm.father)

# homozygous and het parents
rec_trio_criteria = ((phased_sv_tm.complete_trio) &  
                     (phased_sv_tm.proband_entry.GT.is_hom_var()) &
                     (phased_sv_tm.mother_entry.GT.is_het()) &
                     (phased_sv_tm.father_entry.GT.is_het())
                    )  
# hom and unaff are not hom
rec_non_trio_criteria = ((~phased_sv_tm.complete_trio) &  
                        (phased_sv_tm.n_hom_var_unaffected==0) &
                        (phased_sv_tm.proband_entry.GT.is_hom_var())
                        )

phased_sv_tm.annotate_entries(dominant_gt=((rec_trio_criteria) | (rec_non_trio_criteria)),
                              recessive_gt=((rec_trio_criteria) | (rec_non_trio_criteria)))

# Output 1: grab Pathogenic only
path_sv_tm = phased_sv_tm.filter_rows(hl.any(lambda x: x.matches('athogenic'), phased_sv_tm.info.clinical_interpretation) |  # ClinVar P/LP
            (hl.is_defined(phased_sv_tm.info.dbvar_pathogenic)))  # dbVar Pathogenic
path_sv_tm = filter_and_annotate_tm(path_sv_tm, 'P/LP')

# Output 2: grab GD only
gd_sv_tm = phased_sv_tm.filter_rows((hl.is_defined(phased_sv_tm.info.gd_sv_name)))  # GD region  
gd_sv_tm = filter_and_annotate_tm(gd_sv_tm, 'GD')

# Output 3: grab large regions
size_field = [x for x in list(sv_mt.info) if 'passes_SVLEN_filter_' in x][0]
large_sv_tm = phased_sv_tm.filter_rows(phased_sv_tm.info[size_field])
large_sv_tm = filter_and_annotate_tm(large_sv_tm, 'large_region')

# Output 4: grab OMIM AD and XLD
dom_sv_tm = phased_sv_tm.filter_rows((hl.any(lambda x:x.matches('1'), phased_sv_tm.info.OMIM_inheritance_code)) |  # OMIM AD code
                                    (hl.any(lambda x:x.matches('3'), phased_sv_tm.info.OMIM_inheritance_code)))  # OMIM XLD code
dom_sv_tm = filter_and_annotate_tm(dom_sv_tm, 'OMIM_dominant')

# Output 5: grab OMIM AR and XLR
rec_sv_tm = phased_sv_tm.filter_rows((hl.any(lambda x:x.matches('2'), phased_sv_tm.info.OMIM_inheritance_code)) |  # OMIM AD code
                                    (hl.any(lambda x:x.matches('4'), phased_sv_tm.info.OMIM_inheritance_code)))  # OMIM XLD code
rec_sv_tm = filter_and_annotate_tm(rec_sv_tm, 'OMIM_recessive')

# export P/LP TSV
path_sv_tm.entries().flatten().export(os.path.basename(sv_vcf).split('.vcf')[0] + '_path_variants.tsv.gz', delimiter='\t')

# export GD TSV
gd_sv_tm.entries().flatten().export(os.path.basename(sv_vcf).split('.vcf')[0] + '_GD_variants.tsv.gz', delimiter='\t')

# export large region TSV
large_sv_tm.entries().flatten().export(os.path.basename(sv_vcf).split('.vcf')[0] + '_large_regions_variants.tsv.gz', delimiter='\t')

# export OMIM dominant TSV
dom_sv_tm.entries().flatten().export(os.path.basename(sv_vcf).split('.vcf')[0] + '_dominant_variants.tsv.gz', delimiter='\t')

# export OMIM recessive TSV
rec_sv_tm.entries().flatten().export(os.path.basename(sv_vcf).split('.vcf')[0] + '_recessive_variants.tsv.gz', delimiter='\t')
