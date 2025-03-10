###
# Finally forked from hail_filter_clinical_variants_v0.1.py on 1/15/2025
# to output maternal carrier variants as mat_carrier_tsv.

## CHANGE LOG:
'''
1/15/2025:
- fixed bug where empty CLNSIG/CLNREVSTAT (not in ClinVar) gets filtered out
- removed get_transmission function (irrelevant for NIFS)
- added mat_carrier_tsv output (gencc_omim_tm)
- moved all_csqs and gnomad_popmax_af annotations to INFO field
1/17/2025: 
- only include fetal sample in output (mother_entry will be filled)
1/21/2025:
- added CA_from_GT annotation to INFO
2/10/2025:
- added include_all_maternal_carrier_variants parameter
2/25/2025:
- comment out CA_list code because deprecated for CA_from_GT in INFO
- drop CA_from_GT from entries after adding to INFO
2/27/2025:
- revert to all ClinVar P/LP, NOT 2*+ only
3/5/2025:
- fix string matching for ClinVar P/LP output to exclude 'Conflicting'
3/7/2025: 
- added cohort_AC filter
3/10/2025:
- change to ClinVar 1*+ P/LP for clinvar_tsv output
- filter by in gene list in filter_mt (in hail_clinical_helper_functions.py)
'''
###

from pyspark.sql import SparkSession
from clinical_helper_functions import filter_mt, load_split_vep_consequences
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

vcf_file = sys.argv[1]
prefix = sys.argv[2]
cores = sys.argv[3]  # string
mem = int(np.floor(float(sys.argv[4])))
ped_uri = sys.argv[5]
af_threshold = float(sys.argv[6])
ac_threshold = int(sys.argv[7])
gnomad_af_threshold = float(sys.argv[8])
build = sys.argv[9]
pass_filter = ast.literal_eval(sys.argv[10].capitalize())
include_all_maternal_carrier_variants = ast.literal_eval(sys.argv[11].capitalize())

hl.init(min_block_size=128, 
        spark_conf={"spark.executor.cores": cores, 
                    "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                    "spark.driver.cores": "2",
                    "spark.driver.memory": f"{int(np.floor(mem*0.4))}g",
        #             'spark.hadoop.fs.gs.requester.pays.mode': 'CUSTOM',
        #             'spark.hadoop.fs.gs.requester.pays.buckets': 'hail-datasets-us-central1',
        #             'spark.hadoop.fs.gs.requester.pays.project.id': gcp_project,
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

mt = load_split_vep_consequences(vcf_file, build)
header = hl.get_vcf_metadata(vcf_file)

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

# NEW 2/25/2025: Comment out CA_list code because deprecated for CA_from_GT in INFO
# # NEW 1/15/2025: NIFS-specific
# # move CA annotation from entry-level to row-level (FORMAT to INFO)
# mt = mt.annotate_rows(CA_list=hl.array(hl.set(hl.agg.collect(mt.CA).filter(hl.is_defined))))  # CA_list as intermediate field
# mt = mt.annotate_rows(info=mt.info.annotate(
#     CA=hl.or_missing(mt.CA_list.size()>0, mt.CA_list[0])))
# # add CA field to INFO in VCF header
# header['info']['CA'] = header['format']['CA']

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

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

# Mendel errors
all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
all_errors_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)

# Output 1: grab ClinVar only
# NEW 3/5/2025: Fix string matching for ClinVar P/LP output to exclude 'Conflicting'
clinvar_mt = mt.filter_rows(hl.any(lambda x: (x.matches('athogenic')) & (~x.matches('Conflicting')), mt.info.CLNSIG))
clinvar_tm = phased_tm.filter_rows(hl.any(lambda x: (x.matches('athogenic')) & (~x.matches('Conflicting')), phased_tm.info.CLNSIG))
# NEW 1/9/2025: Keep 2*+ ClinVar only
clnrevstat_one_star_plus = ["criteria_provided,_multiple_submitters,_no_conflicts", "criteria_provided,_single_submitter", "practice_guideline", "reviewed_by_expert_panel"]
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
clinvar_tm = filter_mt(clinvar_tm, filter_csq=False, filter_impact=False)  # filter to CANONICAL and/or MANE_PLUS_CLINICAL

# NEW 1/15/2025: liberal set of maternal carrier variants
# Output 2: grab all GenCC_OMIM variants
gencc_omim_tm = phased_tm.explode_rows(phased_tm.vep.transcript_consequences)
# grab anything in GenCC_OMIM gene list and in clusters 1-3 (maternal het clusters)
gencc_omim_tm = gencc_omim_tm.filter_rows((gencc_omim_tm.vep.transcript_consequences.OMIM_inheritance_code!='') &
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
clinvar_tm.filter_cols(clinvar_tm.proband.s.matches('_fetal')).entries().flatten().export(prefix+'_clinvar_variants.tsv.gz', delimiter='\t')

# NEW 2/10/2025: added include_all_maternal_carrier_variants parameter
# export GenCC_OMIM TSV
if include_all_maternal_carrier_variants:
    gencc_omim_tm.filter_cols(gencc_omim_tm.proband.s.matches('_fetal')).entries().flatten().export(prefix+'_mat_carrier_variants.tsv.gz', delimiter='\t')
