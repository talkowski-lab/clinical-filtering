###
# Merged hail_filter_clinical_variants_NIFS_v0.1.py and
# hail_filter_clinical_variants_v0.1.py (small variant version)
# to have a single script for both pipelines, on 4/1/2025.

## CHANGE LOG:
'''
4/7/2025:
- add in_non_par annotation
4/19/2025:
- annotate_trio_matrix function (includes get_mendel_errors, get_transmission)
- filter_by_in_gene_list=False for ClinVar output
'''
###

from pyspark.sql import SparkSession
from clinical_helper_functions import filter_mt, remove_parent_probands_trio_matrix, load_split_vep_consequences, annotate_trio_matrix
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

import argparse
import numpy as np
import ast

# Set up argument parser
parser = argparse.ArgumentParser(description="First-pass filtering clinical variants.")

# Define the arguments
parser.add_argument('--vcf_file', type=str, help="The VCF file path")
parser.add_argument('--prefix', type=str, help="Prefix for output files")
parser.add_argument('--cores', type=str, help="Number of cores to use (as string)")
parser.add_argument('--mem', type=float, help="Memory in GB (as integer, floor converted from float)")
parser.add_argument('--ped_uri', type=str, help="URI to the PED file")
parser.add_argument('--af_threshold', type=float, help="Allele frequency threshold")
parser.add_argument('--ac_threshold', type=int, help="Allele count threshold")
parser.add_argument('--gnomad_af_threshold', type=float, help="gnomAD allele frequency threshold")
parser.add_argument('--build', type=str, help="Genome build (e.g., GRCh38, hg19)")
parser.add_argument('--pass_filter', type=str, help="Whether to apply the PASS filter (True/False)")
parser.add_argument('--include_all_maternal_carrier_variants', type=str, help="Include all maternal carrier variants (True/False)")

# Parse arguments
args = parser.parse_args()

vcf_file = args.vcf_file
prefix = args.prefix
cores = args.cores
mem = np.floor(float(args.mem))  # Ensure it's floored as a float
ped_uri = args.ped_uri
af_threshold = args.af_threshold
ac_threshold = args.ac_threshold
gnomad_af_threshold = args.gnomad_af_threshold
build = args.build
pass_filter = ast.literal_eval(args.pass_filter.capitalize())
include_all_maternal_carrier_variants = ast.literal_eval(args.include_all_maternal_carrier_variants.capitalize())

hl.init(min_block_size=128, 
        spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": "2",
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

mt = load_split_vep_consequences(vcf_file, build)
header = hl.get_vcf_metadata(vcf_file)

# NEW 4/7/2025: add in_non_par annotation
mt = mt.annotate_rows(in_non_par=~(mt.locus.in_autosome_or_par()))

# NEW 1/15/2025: moved all_csqs and gnomad_popmax_af annotations to INFO field
gnomad_fields = [x for x in list(mt.vep.transcript_consequences[0]) if 'gnomAD' in x]
mt = mt.annotate_rows(info=mt.info.annotate(
    all_csqs=hl.array(hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence))),  
    gnomad_popmax_af=hl.max([hl.or_missing(mt.vep.transcript_consequences[0][gnomad_field]!='',
                                    hl.float(mt.vep.transcript_consequences[0][gnomad_field])) 
                             for gnomad_field in gnomad_fields])))
# add all_csqs and gnomad_popmax_af fields to INFO in VCF header
header['info']['all_csqs'] = {'Description': "All unique consequences in vep.transcript_consequences.Consequence for each variant.",
    'Number': '.', 'Type': 'String'}
header['info']['gnomad_popmax_af'] = {'Description': "GnomAD Popmax AF taken across all fields in vep.transcript_consequences containing the string 'gnomAD'.",
    'Number': '1', 'Type': 'Float'}

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
tm = remove_parent_probands_trio_matrix(tm)  # NEW 1/31/2025: Removes redundant "trios"  
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

# Load pedigree as HT for sample annotations
ped_ht = hl.import_table(ped_uri, delimiter='\t').key_by('sample_id')
# NEW 4/19/2025: annotate_trio_matrix function (includes get_mendel_errors, get_transmission)
phased_tm = annotate_trio_matrix(phased_tm, mt, pedigree, ped_ht)

# NEW 1/21/2025: NIFS-specific
# make new row-level annotation, similar to CA, but purely based on GT
phased_tm = phased_tm.annotate_entries(CA_from_GT=hl.if_else(
    (phased_tm.proband_entry.GT.is_het()) & (phased_tm.mother_entry.GT.is_hom_ref()), 0,
    hl.if_else(
        (phased_tm.proband_entry.GT.is_hom_ref()) & (phased_tm.mother_entry.GT.is_het()), 1,
        hl.if_else(
            (phased_tm.proband_entry.GT.is_het()) & (phased_tm.mother_entry.GT.is_het()), 2,
            hl.if_else(
                (phased_tm.proband_entry.GT.is_hom_var()) & (phased_tm.mother_entry.GT.is_het()), 3,
                hl.if_else(
                    (phased_tm.proband_entry.GT.is_het()) & (phased_tm.mother_entry.GT.is_hom_var()), 4,
                    hl.or_missing(
                        (phased_tm.proband_entry.GT.is_hom_var()) & (phased_tm.mother_entry.GT.is_hom_var()), 5
                    )
                )
            )
        )
    ) 
))
phased_tm = phased_tm.annotate_rows(CA_from_GT_list=hl.array(hl.set(hl.agg.collect(phased_tm.CA_from_GT).filter(hl.is_defined))))  # CA_from_GT_list as intermediate field
phased_tm = phased_tm.annotate_rows(info=phased_tm.info.annotate(
    CA_from_GT=hl.or_missing(phased_tm.CA_from_GT_list.size()>0, phased_tm.CA_from_GT_list[0])))
# NEW 2/25/2025: Drop CA_from_GT from entries after adding to INFO
phased_tm = phased_tm.drop('CA_from_GT_list', 'CA_from_GT')
# annotate mt from phased_mt
mt = mt.annotate_rows(info=phased_tm.rows()[mt.row_key].info)
# add CA_from_GT to INFO in VCF header
header['info']['CA_from_GT'] = {'Description': "Cluster assignment, CA, based on fetal/maternal GTs only.",
    'Number': '1', 'Type': 'Int'}

# Output 1: grab ClinVar only
# NEW 3/5/2025: Fix string matching for ClinVar P/LP output to exclude 'Conflicting'
# NEW 4/5/2025: TEST including ClinVar 1*+ P/LP in CLNSIGCONF
# NEW 7/28/2025: Require passing cohort AC and AF filters for CLNSIGCONF P/LP 1*+
clinvar_CLNSIG_P_LP_no_conflicting_mt_cond = hl.any(lambda x: (x.matches('athogenic')) & (~x.matches('Conflicting')), mt.info.CLNSIG)
clinvar_CLNSIGCONF_P_LP_mt_cond = ((hl.any(lambda x: x.matches('athogenic'), mt.info.CLNSIGCONF)) & 
                                   ((mt.info.cohort_AC<=ac_threshold) | (mt.info.cohort_AF<=af_threshold)))
clinvar_mt = mt.filter_rows(clinvar_CLNSIG_P_LP_no_conflicting_mt_cond | clinvar_CLNSIGCONF_P_LP_mt_cond)

clinvar_CLNSIG_P_LP_no_conflicting_tm_cond = hl.any(lambda x: (x.matches('athogenic')) & (~x.matches('Conflicting')), phased_tm.info.CLNSIG)
clinvar_CLNSIGCONF_P_LP_tm_cond = ((hl.any(lambda x: x.matches('athogenic'), phased_tm.info.CLNSIGCONF)) & 
                                   ((mt.info.cohort_AC<=ac_threshold) | (mt.info.cohort_AF<=af_threshold)))
clinvar_tm = phased_tm.filter_rows(clinvar_CLNSIG_P_LP_no_conflicting_tm_cond | clinvar_CLNSIGCONF_P_LP_tm_cond)

# NEW 1/9/2025: Keep 2*+ ClinVar only
clnrevstat_one_star_plus = [['practice_guideline'], ['reviewed_by_expert_panel'], ['criteria_provided','_multiple_submitters','_no_conflicts'], ['criteria_provided','_conflicting_classifications'], ['criteria_provided','_single_submitter']]
clnrevstat_two_star_plus = [['practice_guideline'], ['reviewed_by_expert_panel'], ['criteria_provided', '_multiple_submitters', '_no_conflicts']]
# NEW 2/27/2025: Revert to all ClinVar P/LP, NOT 2*+ only
# NEW 3/10/2025: Change to ClinVar 1*+ P/LP
clinvar_mt = clinvar_mt.filter_rows(hl.any([clinvar_mt.info.CLNREVSTAT==category for category in clnrevstat_one_star_plus]))
clinvar_tm = clinvar_tm.filter_rows(hl.any([clinvar_tm.info.CLNREVSTAT==category for category in clnrevstat_one_star_plus]))

clinvar_tm = clinvar_tm.filter_entries((clinvar_tm.proband_entry.GT.is_non_ref()) | 
                                   (clinvar_tm.mother_entry.GT.is_non_ref()) |
                                   (clinvar_tm.father_entry.GT.is_non_ref()))
clinvar_tm = clinvar_tm.annotate_rows(variant_category='ClinVar_P/LP')
clinvar_tm = clinvar_tm.explode_rows(clinvar_tm.vep.transcript_consequences)
clinvar_tm = filter_mt(clinvar_tm, filter_csq=False, filter_impact=False, filter_by_in_gene_list=False)  # filter to CANONICAL and/or MANE_PLUS_CLINICAL

# NEW 1/15/2025: liberal set of maternal carrier variants
# Output 2: grab all GenCC_OMIM variants
gencc_omim_tm = phased_tm.explode_rows(phased_tm.vep.transcript_consequences)
# grab anything in GenCC_OMIM gene list and in clusters 1-3 (maternal het clusters)
gencc_omim_tm = gencc_omim_tm.filter_rows((gencc_omim_tm.vep.transcript_consequences.inheritance_code!='') &
                       (hl.array([1, 2, 3]).contains(gencc_omim_tm.info.CA_from_GT)))
gencc_omim_tm = gencc_omim_tm.annotate_rows(variant_category='GenCC_OMIM')
gencc_omim_tm = filter_mt(gencc_omim_tm, filter_csq=False, filter_impact=False)  # filter to CANONICAL and/or MANE_PLUS_CLINICAL

# filter out ClinVar benign
# NEW 1/10/2025 filter out 2*+ benign only!
# NEW: 1/15/2025: fixed bug where empty CLNSIG/CLNREVSTAT (not in ClinVar) gets filtered out
is_clinvar_benign = ((hl.is_defined(mt.info.CLNSIG)) & (hl.any(lambda x: x.matches('enign'), mt.info.CLNSIG)))
is_clnrevstat_two_star_plus = ((hl.is_defined(mt.info.CLNREVSTAT)) & 
                (hl.any([mt.info.CLNREVSTAT==category for category in clnrevstat_two_star_plus])))
mt = mt.filter_rows(is_clinvar_benign & is_clnrevstat_two_star_plus, keep=False)

# filter PASS
if (pass_filter):
    mt = mt.filter_rows(mt.filters.size()==0)

# filter out variants containing only these consequences
exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']

mt = mt.filter_rows(hl.set(exclude_csqs).intersection(hl.set(mt.info.all_csqs)).size()!=mt.info.all_csqs.size())

# filter by cohort_AC, cohort AF and gnomAD AF
# NEW 3/7/2025: added cohort_AC filter
mt = mt.filter_rows((mt.info.cohort_AC<=ac_threshold) | (mt.info.cohort_AF<=af_threshold))
mt = mt.filter_rows((mt.info.gnomad_popmax_af<=gnomad_af_threshold) | (hl.is_missing(mt.info.gnomad_popmax_af)))

# export intermediate VCF
hl.export_vcf(mt, prefix+'_clinical.vcf.bgz', metadata=header)

# export ClinVar VCF
hl.export_vcf(clinvar_mt, prefix+'_clinvar_variants.vcf.bgz', metadata=header, tabix=True)

# NEW 1/17/2025: only include fetal sample in output (mother_entry will be filled)
# export ClinVar TSV
clinvar_tm.entries().flatten().export(prefix+'_clinvar_variants.tsv.gz', delimiter='\t')

# NEW 2/10/2025: added include_all_maternal_carrier_variants parameter
# export GenCC_OMIM TSV
if include_all_maternal_carrier_variants:
    gencc_omim_tm.entries().flatten().export(prefix+'_mat_carrier_variants.tsv.gz', delimiter='\t')