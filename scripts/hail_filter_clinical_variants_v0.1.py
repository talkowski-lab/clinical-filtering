###
# Copied over from hail_filter_clinical_variants_NIFS_v0.1.py on 1/27/2025
# to incorporate updates made to NIFS pipeline

## CHANGE LOG:
'''
1/27/2025:
- removed all code related to CA (NIFS-specific)
- edited gencc_omim_tm filtering to use mother_entry.GT instead of CA for mother het status, filter by non-ref GT in proband
- removed "# NEW 1/17/2025: only include fetal sample in output (mother_entry will be filled)" code logic for outputs (NIFS-specific)
1/31/2025:
- added remove_parent_probands_trio_matrix function --> removes redundant "trios"
2/10/2025:
- added include_all_maternal_carrier_variants parameter
2/27/2025:
- revert to all ClinVar P/LP, NOT 2*+ only
3/5/2025:
- fix string matching for ClinVar P/LP output to exclude 'Conflicting'
3/7/2025: 
- added cohort_AC filter
'''
###

from pyspark.sql import SparkSession
from clinical_helper_functions import filter_mt, remove_parent_probands_trio_matrix, load_split_vep_consequences
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

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)
pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
tm = remove_parent_probands_trio_matrix(tm)  # NEW 1/31/2025: Removes redundant "trios"  
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

# Mendel errors
all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
all_errors_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)

# Output 1: grab ClinVar only
# NEW 3/5/2025: Fix string matching for ClinVar P/LP output to exclude 'Conflicting'
clinvar_mt = mt.filter_rows(hl.any(lambda x: (x.matches('athogenic')) & (~x.matches('Conflicting')), mt.info.CLNSIG))
clinvar_tm = phased_tm.filter_rows(hl.any(lambda x: (x.matches('athogenic')) & (~x.matches('Conflicting')), phased_tm.info.CLNSIG))
# NEW 1/9/2025: keep 2*+ ClinVar only
clinvar_two_star_plus = [['practice_guideline'], ['reviewed_by_expert_panel'], ['criteria_provided', '_multiple_submitters', '_no_conflicts']]
# NEW 2/27/2025: Revert to all ClinVar P/LP, NOT 2*+ only
# clinvar_mt = clinvar_mt.filter_rows(hl.any([clinvar_mt.info.CLNREVSTAT==category for category in clinvar_two_star_plus]))
# clinvar_tm = clinvar_tm.filter_rows(hl.any([clinvar_tm.info.CLNREVSTAT==category for category in clinvar_two_star_plus]))

clinvar_tm = clinvar_tm.filter_entries((clinvar_tm.proband_entry.GT.is_non_ref()) | 
                                   (clinvar_tm.mother_entry.GT.is_non_ref()) |
                                   (clinvar_tm.father_entry.GT.is_non_ref()))
clinvar_tm = clinvar_tm.annotate_rows(variant_category='ClinVar_P/LP')
clinvar_tm = clinvar_tm.explode_rows(clinvar_tm.vep.transcript_consequences)
clinvar_tm = filter_mt(clinvar_tm, filter_csq=False, filter_impact=False)  # filter to CANONICAL and/or MANE_PLUS_CLINICAL

# NEW 1/15/2025: liberal set of maternal carrier variants
# Output 2: grab all GenCC_OMIM variants
gencc_omim_tm = phased_tm.explode_rows(phased_tm.vep.transcript_consequences)
# NEW 1/27/2025: grab anything in GenCC_OMIM gene list and het in mother
gencc_omim_tm = gencc_omim_tm.filter_rows(gencc_omim_tm.vep.transcript_consequences.inheritance_code!='')
gencc_omim_tm = gencc_omim_tm.filter_entries((gencc_omim_tm.mother_entry.GT.is_het()) & 
                                             (gencc_omim_tm.proband_entry.GT.is_non_ref()))
gencc_omim_tm = gencc_omim_tm.annotate_rows(variant_category='GenCC_OMIM')
gencc_omim_tm = filter_mt(gencc_omim_tm, filter_csq=False, filter_impact=False)  # filter to CANONICAL and/or MANE_PLUS_CLINICAL

# filter out ClinVar benign
# NEW 1/10/2025 filter out 2*+ benign only!
# NEW: 1/15/2025: fixed bug where empty CLNSIG/CLNREVSTAT (not in ClinVar) gets filtered out
is_clinvar_benign = ((hl.is_defined(mt.info.CLNSIG)) & (hl.any(lambda x: x.matches('enign'), mt.info.CLNSIG)))
is_clinvar_two_star_plus = ((hl.is_defined(mt.info.CLNREVSTAT)) & 
                (hl.any([mt.info.CLNREVSTAT==category for category in clinvar_two_star_plus])))
mt = mt.filter_rows(is_clinvar_benign & is_clinvar_two_star_plus, keep=False)

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

output_filename = f"{prefix}_{variant_types}_comp_hets_xlr_hom_var_mat_carrier.tsv.gz"
if len(output_filename)>(os.pathconf('/', 'PC_NAME_MAX')-len('/cromwell_root/.')):  # if filename too long
    output_filename = f"{variant_types}_comp_hets_xlr_hom_var_mat_carrier.tsv.gz"

# export intermediate VCF
hl.export_vcf(mt, prefix+'_clinical.vcf.bgz', metadata=header)

# export ClinVar VCF
hl.export_vcf(clinvar_mt, prefix+'_clinvar_variants.vcf.bgz', metadata=header, tabix=True)

# export ClinVar TSV
clinvar_tm.entries().flatten().export(prefix+'_clinvar_variants.tsv.gz', delimiter='\t')

# NEW 2/10/2025: added include_all_maternal_carrier_variants parameter
# export GenCC_OMIM TSV
if include_all_maternal_carrier_variants:
    gencc_omim_tm.entries().flatten().export(prefix+'_mat_carrier_variants.tsv.gz', delimiter='\t')
