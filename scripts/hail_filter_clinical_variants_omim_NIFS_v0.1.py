###
# Finally forked from hail_filter_clinical_variants_omim_v0.1.py on 1/15/2025
# to match other scripts.

## CHANGE LOG:
'''
1/15/2025:
- removed get_transmission function (irrelevant for NIFS)
- commented out all_csqs and gnomad_popmax_af annotations because now annotated (in INFO) in hail_filter_clinical_variants_NIFS_v0.1.py :)
- changed has_low_or_modifier_impact to is_moderate_or_high_impact inverse logic (should have the same result)
1/17/2025: 
- add recessive_tsv output
- only include fetal sample in output (mother_entry will be filled)
1/27/2025:
- filter by AD before outputting recessive_tsv
3/4/2025:
- change OMIM_recessive/OMIM_dominant to just recessive/dominant
3/10/2025: 
- cohort_AC OR cohort_AF filter
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

import argparse
import numpy as np
import ast

# Create the parser
parser = argparse.ArgumentParser(description='Process some genomic data.')

# Add arguments
parser.add_argument('--vcf_file', type=str, help='Input VCF file')
parser.add_argument('--prefix', type=str, help='Output prefix')
parser.add_argument('--cores', type=str, help='Number of cores')
parser.add_argument('--mem', type=float, help='Memory (as float)')
parser.add_argument('--ped_uri', type=str, help='PED file URI')
parser.add_argument('--ac_rec_threshold', type=int, help='AC rec threshold')
parser.add_argument('--af_rec_threshold', type=float, help='AF rec threshold')
parser.add_argument('--ac_dom_threshold', type=int, help='AC dom threshold')
parser.add_argument('--af_dom_threshold', type=float, help='AF dom threshold')
parser.add_argument('--am_rec_threshold', type=float, help='AM rec threshold')
parser.add_argument('--am_dom_threshold', type=float, help='AM dom threshold')
parser.add_argument('--mpc_rec_threshold', type=float, help='MPC rec threshold')
parser.add_argument('--mpc_dom_threshold', type=float, help='MPC dom threshold')
parser.add_argument('--gnomad_af_rec_threshold', type=float, help='gnomAD AF rec threshold')
parser.add_argument('--gnomad_af_dom_threshold', type=float, help='gnomAD AF dom threshold')
parser.add_argument('--loeuf_v2_threshold', type=float, help='LOEUF v2 threshold')
parser.add_argument('--loeuf_v4_threshold', type=float, help='LOEUF v4 threshold')
parser.add_argument('--build', type=str, help='Build version')
parser.add_argument('--ad_alt_threshold', type=int, help='AD ALT threshold')
parser.add_argument('--include_not_omim', type=str, help='Include not OMIM (True/False)')
parser.add_argument('--spliceAI_threshold', type=float, help='SpliceAI threshold')
parser.add_argument('--rec_gene_list_tsv', type=str, help='rec gene list TSV')
parser.add_argument('--dom_gene_list_tsv', type=str, help='dom gene list TSV')

args = parser.parse_args()

mem = int(np.floor(args.mem))
include_not_omim = ast.literal_eval(args.include_not_omim.capitalize())

vcf_file = args.vcf_file
prefix = args.prefix
cores = args.cores
ped_uri = args.ped_uri
ac_rec_threshold = args.ac_rec_threshold
af_rec_threshold = args.af_rec_threshold
ac_dom_threshold = args.ac_dom_threshold
af_dom_threshold = args.af_dom_threshold
am_rec_threshold = args.am_rec_threshold
am_dom_threshold = args.am_dom_threshold
mpc_rec_threshold = args.mpc_rec_threshold
mpc_dom_threshold = args.mpc_dom_threshold
gnomad_af_rec_threshold = args.gnomad_af_rec_threshold
gnomad_af_dom_threshold = args.gnomad_af_dom_threshold
loeuf_v2_threshold = args.loeuf_v2_threshold
loeuf_v4_threshold = args.loeuf_v4_threshold
build = args.build
ad_alt_threshold = args.ad_alt_threshold
spliceAI_threshold = args.spliceAI_threshold
rec_gene_list_tsv = args.rec_gene_list_tsv
dom_gene_list_tsv = args.dom_gene_list_tsv

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

# NEW 1/9/2025: moved gnomad_popmax_af
# NEW 1/15/2025: commented out because now annotated (in INFO) in hail_filter_clinical_variants_NIFS_v0.1.py :)
# mt = mt.annotate_rows(all_csqs=hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence)))

# Phasing
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
cropped_ped_uri = f"{os.path.basename(ped_uri).split('.ped')[0]}_crop.ped"
tmp_ped.to_csv(cropped_ped_uri, sep='\t', index=False)
pedigree = hl.Pedigree.read(cropped_ped_uri, delimiter='\t')

tm = hl.trio_matrix(mt, pedigree, complete_trios=False)
phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')

# Mendel errors
all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(mt['GT'], pedigree)
all_errors_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)

gene_phased_tm = phased_tm.explode_rows(phased_tm.vep.transcript_consequences)

# NEW 1/9/2025: moved gnomad_popmax_af to after exploding by transcript
# NEW 1/15/2025: commented out because now annotated (in INFO) in hail_filter_clinical_variants_NIFS_v0.1.py :)
# gnomad_fields = [x for x in list(gene_phased_tm.vep.transcript_consequences) if 'gnomAD' in x]
# gene_phased_tm = gene_phased_tm.annotate_rows(
#     gnomad_popmax_af=hl.max([hl.or_missing(gene_phased_tm.vep.transcript_consequences[gnomad_field]!='',
#                                     hl.float(gene_phased_tm.vep.transcript_consequences[gnomad_field])) 
#                              for gnomad_field in gnomad_fields]))

gene_phased_tm = filter_mt(gene_phased_tm)

# annotate spliceAI score if missing
if 'spliceAI_score' not in list(gene_phased_tm.vep.transcript_consequences):
    fields = ['DS_AG','DS_AL','DS_DG','DS_DL']
    gene_phased_tm = gene_phased_tm.annotate_rows(vep=gene_phased_tm.vep.annotate(
        transcript_consequences=(gene_phased_tm.vep.transcript_consequences.annotate(
            spliceAI_score=hl.str(hl.max([
                hl.or_missing(gene_phased_tm.vep.transcript_consequences[field]!='', 
                            hl.float(gene_phased_tm.vep.transcript_consequences[field])) 
                for field in fields])))))
    )    
    # replace Hail missing values with empty string
    gene_phased_tm = gene_phased_tm.annotate_rows(vep=gene_phased_tm.vep.annotate(
        transcript_consequences=(gene_phased_tm.vep.transcript_consequences.annotate(
            spliceAI_score=hl.if_else(hl.is_missing(gene_phased_tm.vep.transcript_consequences.spliceAI_score),
                                    '', gene_phased_tm.vep.transcript_consequences.spliceAI_score
            )))))

# filter by spliceAI score
# NEW 1/7/2025 only apply on splice variants --> updated below
# NEW 1/10/2025 only apply on splice variants with no other HIGH/MODERATE impact consequences
splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
has_splice_var = (hl.set(splice_vars).intersection(
            hl.set(gene_phased_tm.vep.transcript_consequences.Consequence)).size()>0)
is_splice_var_only = (hl.set(splice_vars).intersection(
            hl.set(gene_phased_tm.vep.transcript_consequences.Consequence)).size()==
                hl.set(gene_phased_tm.vep.transcript_consequences.Consequence).size())

fails_spliceAI_score = (hl.if_else(gene_phased_tm.vep.transcript_consequences.spliceAI_score=='', 1, 
                hl.float(gene_phased_tm.vep.transcript_consequences.spliceAI_score))<spliceAI_threshold)

# NEW 1/15/2025 changed has_low_or_modifier_impact to is_moderate_or_high_impact inverse logic
# has_low_or_modifier_impact = (hl.array(['LOW','MODIFIER']).contains(gene_phased_tm.vep.transcript_consequences.IMPACT))
# gene_phased_tm = gene_phased_tm.filter_rows((is_splice_var_only | 
#                                              (has_splice_var & has_low_or_modifier_impact)) & fails_spliceAI_score, keep=False)

is_moderate_or_high_impact = (hl.array(['HIGH','MODERATE']).contains(gene_phased_tm.vep.transcript_consequences.IMPACT))
gene_phased_tm = gene_phased_tm.filter_rows((is_splice_var_only | 
                                             (has_splice_var & ~is_moderate_or_high_impact)) & fails_spliceAI_score, keep=False)

# Output 2: OMIM Recessive
# Filter by gene list(s)
if rec_gene_list_tsv!='NA':
    gene_list_uris = pd.read_csv(rec_gene_list_tsv, sep='\t', header=None).set_index(0)[1].to_dict()
    gene_lists = {gene_list_name: pd.read_csv(uri, sep='\t', header=None)[0].tolist() 
                for gene_list_name, uri in gene_list_uris.items()}

    gene_phased_tm = gene_phased_tm.annotate_rows(
        gene_lists=hl.array([hl.or_missing(hl.array(gene_list).contains(gene_phased_tm.vep.transcript_consequences.SYMBOL), gene_list_name) 
            for gene_list_name, gene_list in gene_lists.items()]).filter(hl.is_defined))

    in_rec_gene_list = (gene_phased_tm.gene_lists.size()>0)
else:
    in_rec_gene_list = False

not_in_omim = (gene_phased_tm.vep.transcript_consequences.inheritance_code=='')
# OMIM recessive code
omim_rec_code = (gene_phased_tm.vep.transcript_consequences.inheritance_code.matches('2'))
# OMIM XLR code
omim_xlr_code = (gene_phased_tm.vep.transcript_consequences.inheritance_code.matches('4'))
# NEW 3/10/2025: cohort_AC OR cohort_AF filter
passes_ac_af_rec = ((gene_phased_tm.info.cohort_AC<=ac_rec_threshold) | (gene_phased_tm.info.cohort_AF<=af_rec_threshold))
# gnomAD AF popmax filter
passes_gnomad_af_rec = ((gene_phased_tm.info.gnomad_popmax_af<=gnomad_af_rec_threshold) | (hl.is_missing(gene_phased_tm.info.gnomad_popmax_af)))
# MPC filter
passes_mpc_rec = ((gene_phased_tm.info.MPC>=mpc_rec_threshold) | (hl.is_missing(gene_phased_tm.info.MPC)))
# AlphaMissense filter
# NEW 1/7/2025 only apply on missense variants
is_missense_var = (hl.set(['missense_variant']).intersection(
            hl.set(gene_phased_tm.vep.transcript_consequences.Consequence)).size()>0)
passes_alpha_missense_score = (hl.if_else(gene_phased_tm.vep.transcript_consequences.am_pathogenicity=='', 1, 
                hl.float(gene_phased_tm.vep.transcript_consequences.am_pathogenicity))>=am_rec_threshold)
passes_alpha_missense = ((is_missense_var & passes_alpha_missense_score) | (~is_missense_var))

if include_not_omim:
    omim_rec_gene_phased_tm = gene_phased_tm.filter_rows(
        (passes_ac_af_rec) &
        (passes_alpha_missense) &
        (
            omim_rec_code |
            omim_xlr_code |
            in_rec_gene_list |
            (
                not_in_omim &
                passes_mpc_rec &
                passes_gnomad_af_rec
            )
        )
    )
else:
    omim_rec_gene_phased_tm = gene_phased_tm.filter_rows(
        (passes_ac_af_rec) &
        (passes_alpha_missense) &
            (omim_rec_code | omim_xlr_code | in_rec_gene_list)        
    )

omim_rec_gene_phased_tm = (omim_rec_gene_phased_tm.group_rows_by(omim_rec_gene_phased_tm.locus, omim_rec_gene_phased_tm.alleles)
    .aggregate_rows(vep = hl.agg.collect(omim_rec_gene_phased_tm.vep))).result()

fields = list(omim_rec_gene_phased_tm.vep.transcript_consequences[0])
new_csq = omim_rec_gene_phased_tm.vep.transcript_consequences.scan(lambda i, j: 
                                      hl.str('|').join(hl.array([i]))
                                      +','+hl.str('|').join(hl.array([j[col] if col!='Consequence' else 
                                                                  hl.str('&').join(j[col]) 
                                                                  for col in list(fields)])), '')[-1][1:]
omim_rec_gene_phased_tm = omim_rec_gene_phased_tm.annotate_rows(CSQ=new_csq)
omim_rec_mt = mt.semi_join_rows(omim_rec_gene_phased_tm.rows())
omim_rec_mt = omim_rec_mt.annotate_rows(info=omim_rec_mt.info.annotate(CSQ=omim_rec_gene_phased_tm.rows()[omim_rec_mt.row_key].CSQ))

# Output 3: OMIM Dominant
# Filter by gene list(s)
if dom_gene_list_tsv!='NA':
    gene_list_uris = pd.read_csv(dom_gene_list_tsv, sep='\t', header=None).set_index(0)[1].to_dict()
    gene_lists = {gene_list_name: pd.read_csv(uri, sep='\t', header=None)[0].tolist() 
                for gene_list_name, uri in gene_list_uris.items()}
    # overrides recessive gene_lists but already saved as intermediate
    gene_phased_tm = gene_phased_tm.annotate_rows(
        gene_lists=hl.array([hl.or_missing(hl.array(gene_list).contains(gene_phased_tm.vep.transcript_consequences.SYMBOL), gene_list_name) 
            for gene_list_name, gene_list in gene_lists.items()]).filter(hl.is_defined))

    in_dom_gene_list = (gene_phased_tm.gene_lists.size()>0)
else:
    in_dom_gene_list = False

not_in_omim = (gene_phased_tm.vep.transcript_consequences.inheritance_code=='')
# OMIM dominant code
omim_dom_code = (gene_phased_tm.vep.transcript_consequences.inheritance_code.matches('1')) 
# OMIM XLD code
omim_xld_code = (gene_phased_tm.vep.transcript_consequences.inheritance_code.matches('3'))
# NEW 3/10/2025: cohort_AC OR cohort_AF filter
passes_ac_af_dom = ((gene_phased_tm.info.cohort_AC<=ac_dom_threshold) | (gene_phased_tm.info.cohort_AF<=af_dom_threshold))
# gnomAD AF popmax filter
passes_gnomad_af_dom = ((gene_phased_tm.info.gnomad_popmax_af<=gnomad_af_dom_threshold) | (hl.is_missing(gene_phased_tm.info.gnomad_popmax_af)))
# MPC filter
passes_mpc_dom = ((gene_phased_tm.info.MPC>=mpc_dom_threshold) | (hl.is_missing(gene_phased_tm.info.MPC)))
# AlphaMissense filter
# NEW 1/7/2025 only apply on missense variants
is_missense_var = (hl.set(['missense_variant']).intersection(
            hl.set(gene_phased_tm.vep.transcript_consequences.Consequence)).size()>0)
passes_alpha_missense_score = (hl.if_else(gene_phased_tm.vep.transcript_consequences.am_pathogenicity=='', 1, 
                hl.float(gene_phased_tm.vep.transcript_consequences.am_pathogenicity))>=am_dom_threshold)
passes_alpha_missense = ((is_missense_var & passes_alpha_missense_score) | (~is_missense_var))
# LOEUF v2/v4 filters
passes_loeuf_v2 = (hl.if_else(gene_phased_tm.vep.transcript_consequences.LOEUF_v2=='', 0, 
                        hl.float(gene_phased_tm.vep.transcript_consequences.LOEUF_v2))<=loeuf_v2_threshold)

passes_loeuf_v4 = (hl.if_else(gene_phased_tm.vep.transcript_consequences.LOEUF_v4=='', 0, 
                        hl.float(gene_phased_tm.vep.transcript_consequences.LOEUF_v4))<=loeuf_v4_threshold)

if include_not_omim:
    omim_dom = gene_phased_tm.filter_rows(
        (passes_ac_af_dom) &
        (passes_gnomad_af_dom) &
        (passes_alpha_missense) &
        (
            omim_dom_code |
            omim_xld_code |
            in_dom_gene_list |
            (
                not_in_omim &
                passes_mpc_dom &
                (passes_loeuf_v2 | passes_loeuf_v4)
            )
        )
    )
else:
    omim_dom = gene_phased_tm.filter_rows(
        (passes_ac_af_dom) &
        (passes_gnomad_af_dom) &
        (passes_alpha_missense) &
        (omim_dom_code | omim_xld_code | in_dom_gene_list)
    )

omim_dom = omim_dom.filter_entries((omim_dom.proband_entry.GT.is_non_ref()) | 
                                   (omim_dom.mother_entry.GT.is_non_ref()) |
                                   (omim_dom.father_entry.GT.is_non_ref()))

# filter by AD of alternate allele in trio
omim_dom = omim_dom.filter_entries((omim_dom.proband_entry.AD[1]>=ad_alt_threshold) |
                                   (omim_dom.mother_entry.AD[1]>=ad_alt_threshold) |
                                   (omim_dom.father_entry.AD[1]>=ad_alt_threshold))

# variant must be in at least 1 trio
omim_dom = omim_dom.filter_rows(hl.agg.count_where((hl.is_defined(omim_dom.proband_entry.GT)) |
                                                    (hl.is_defined(omim_dom.mother_entry.GT)) |
                                                    (hl.is_defined(omim_dom.father_entry.GT)))>0)
omim_dom = omim_dom.annotate_rows(variant_category='dominant')

# NEW 1/17/2025: export OMIM Recessive TSV
omim_rec = gene_phased_tm.semi_join_rows(omim_rec_mt.rows())

# NEW 1/27/2025: filter by AD before outputting recessive_tsv
# filter by AD of alternate allele in trio
omim_rec = omim_rec.filter_entries(omim_rec.proband_entry.AD[1]>=ad_alt_threshold)
# variant must be in at least 1 trio
omim_rec = omim_rec.filter_rows(hl.agg.count_where(hl.is_defined(omim_rec.proband_entry.GT))>0)
omim_rec = omim_rec.annotate_rows(variant_category='recessive')

# export OMIM Recessive VCF
hl.export_vcf(omim_rec_mt, prefix+'_recessive.vcf.bgz', metadata=header, tabix=True)

# export OMIM Dominant TSV
# NEW 1/17/2025: only include fetal sample in output (mother_entry will be filled)
omim_dom.filter_cols(omim_dom.proband.s.matches('_fetal')).entries().flatten().export(prefix+'_dominant.tsv.gz', delimiter='\t')

# export OMIM Recessive TSV
# NEW 1/17/2025: only include fetal sample in output (mother_entry will be filled)
omim_rec.filter_cols(omim_rec.proband.s.matches('_fetal')).entries().flatten().export(prefix+'_recessive.tsv.gz', delimiter='\t')