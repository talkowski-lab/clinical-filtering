###
# New script based on some code from filterClinicalVariantsSV_v0.1.wdl, 
# tested in test_sv_annotations.ipynb. 
# Created 1/30/2025.

## CHANGE LOG:
'''
1/30/2025:
- added dominant_gt and recessive_gt annotations
1/31/2025:
- added remove_parent_probands_trio_matrix function --> removes redundant "trios"
2/3/2025:
- changed complete_trio annotation to trio_status to match comphet script
2/5/2025:
- moved ped format standardization to before ped_ht
2/19/2025:
- use dominant_freq and recessive_freq flag in INFO to filter for dominant and recessive outputs
- annotate categories for merged output, instead of filtering for five separate outputs
- get affected/unaffected counts at family-level
2/20/2025:
- rewrote variant_category annotation because using append was messing things up in Hail
5/29/2025:
- annotate family unaffected/affected AC
- move affected/unaffected fields out of INFO
6/2/2025:
- only consider restrictive genes for dominant and recessive inheritance_code filtering
7/14/2025:
- restrict large regions to not BND or CPX, impacts a gene from restrictive fields
- restrict large_region only outputs to private to a family
- new 'genic' category for any restrictive fields' genes that are in a gene list
'''
###

from clinical_helper_functions import remove_parent_probands_trio_matrix, get_transmission
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

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

def filter_and_annotate_tm(tm, variant_category=None, filter_type='trio'):
    '''
    Filter entries by presence in trio --> annotate by variant_category (if given)
    filter_type can be 'trio' or 'proband'
    '''
    if filter_type=='trio':
        tm = tm.filter_entries((tm.proband_entry.GT.is_non_ref()) | 
                                    (tm.mother_entry.GT.is_non_ref()) |
                                    (tm.father_entry.GT.is_non_ref()))
    if filter_type=='proband':
        tm = tm.filter_entries(tm.proband_entry.GT.is_non_ref())
    if variant_category:
        tm = tm.annotate_rows(variant_category=variant_category)
    return tm

sv_mt = hl.import_vcf(sv_vcf, force_bgz=sv_vcf.split('.')[-1] in ['gz', 'bgz'], 
    reference_genome=genome_build, array_elements_required=False, call_fields=[])
header = hl.get_vcf_metadata(sv_vcf)

# Annotate affected status/phenotype from pedigree
# NEW 2/5/2025: Moved ped format standardization to before ped_ht
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
cropped_ped_uri = f"{os.path.basename(ped_uri).split('.ped')[0]}_crop.ped"
tmp_ped.to_csv(cropped_ped_uri, sep='\t', index=False)

ped_ht = hl.import_table(cropped_ped_uri, delimiter='\t').key_by('sample_id')
sv_mt = sv_mt.annotate_cols(phenotype=ped_ht[sv_mt.s].phenotype)

# Phasing
pedigree = hl.Pedigree.read(cropped_ped_uri, delimiter='\t')

# NEW 2/19/2025: Get affected/unaffected counts at family-level
fam_sv_mt = sv_mt.annotate_cols(family_id=ped_ht[sv_mt.s].family_id)
grouped_fam_sv_mt = fam_sv_mt.group_cols_by(fam_sv_mt.family_id).aggregate(**{
    "n_family_het_unaffected": hl.agg.filter(fam_sv_mt.phenotype=='1', hl.agg.sum(fam_sv_mt.GT.is_het())),
    "n_family_hom_var_unaffected": hl.agg.filter(fam_sv_mt.phenotype=='1', hl.agg.sum(fam_sv_mt.GT.is_hom_var())),
    "n_family_het_affected": hl.agg.filter(fam_sv_mt.phenotype=='2', hl.agg.sum(fam_sv_mt.GT.is_het())),
    "n_family_hom_var_affected": hl.agg.filter(fam_sv_mt.phenotype=='2', hl.agg.sum(fam_sv_mt.GT.is_hom_var()))
})

# NEW 5/29/2025: Annotate family unaffected/affected AC  
grouped_fam_sv_mt = grouped_fam_sv_mt.annotate_entries(**{
    "family_AC_unaffected": grouped_fam_sv_mt.n_family_het_unaffected + 2*grouped_fam_sv_mt.n_family_hom_var_unaffected,
    "family_AC_affected": grouped_fam_sv_mt.n_family_het_affected + 2*grouped_fam_sv_mt.n_family_hom_var_affected
})

sv_tm = hl.trio_matrix(sv_mt, pedigree, complete_trios=False)
sv_tm = remove_parent_probands_trio_matrix(sv_tm)  # NEW 1/31/2025: Removes redundant "trios"  
phased_sv_tm = hl.experimental.phase_trio_matrix_by_transmission(sv_tm, call_field='GT', phased_call_field='PBT_GT')

# Mendel errors
all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(sv_mt['GT'], pedigree)
all_errors_sv_mt = all_errors.key_by().to_matrix_table(row_key=['locus','alleles'], col_key=['s'], row_fields=['fam_id'])
phased_sv_tm = phased_sv_tm.annotate_entries(mendel_code=all_errors_sv_mt[phased_sv_tm.row_key, phased_sv_tm.col_key].mendel_code)
phased_sv_tm = get_transmission(phased_sv_tm)

# Annotate if complete trio or not
# NEW 2/3/2025: Change complete_trio annotation to trio_status to match comphet script
complete_trio_probands = [trio.s for trio in pedigree.complete_trios()]
phased_sv_tm = phased_sv_tm.annotate_cols(trio_status=hl.if_else(phased_sv_tm.fam_id=='-9', 'not_in_pedigree', 
                                                   hl.if_else(hl.array(complete_trio_probands).contains(phased_sv_tm.id), 'trio', 'non_trio')))

# Annotate affected status/phenotype from pedigree
phased_sv_tm = phased_sv_tm.annotate_cols(
    proband=phased_sv_tm.proband.annotate(
        phenotype=ped_ht[phased_sv_tm.proband.s].phenotype),
    mother=phased_sv_tm.mother.annotate(
            phenotype=ped_ht[phased_sv_tm.mother.s].phenotype),
    father=phased_sv_tm.father.annotate(
            phenotype=ped_ht[phased_sv_tm.father.s].phenotype))

# NEW 5/29/2025: Move affected/unaffected fields out of INFO
cohort_affected_cols = [field for field in list(sv_mt.info) if '_affected' in field or '_unaffected' in field]
sv_mt = sv_mt.annotate_rows(**{field: sv_mt.info[field] for field in cohort_affected_cols})

phased_sv_tm = phased_sv_tm.annotate_rows(**{col: sv_mt.rows()[phased_sv_tm.row_key][col] 
                                             for col in cohort_affected_cols}  | 
                                            {'info': phased_sv_tm.info.drop(*cohort_affected_cols)})
family_affected_cols = [field for field in list(grouped_fam_sv_mt.entry) if '_affected' in field or '_unaffected' in field]
phased_sv_tm = phased_sv_tm.annotate_entries(**{col: grouped_fam_sv_mt[phased_sv_tm.row_key, phased_sv_tm.fam_id][col]
                                             for col in family_affected_cols})

## Annotate dominant_gt and recessive_gt
# denovo
dom_trio_criteria = ((phased_sv_tm.trio_status=='trio') &  
                     (phased_sv_tm.mendel_code==2))
# het absent in unaff
dom_non_trio_criteria = ((phased_sv_tm.trio_status!='trio') & 
                        (phased_sv_tm.n_cohort_hom_var_unaffected==0) & 
                        (phased_sv_tm.n_cohort_het_unaffected==0) & 
                        (phased_sv_tm.proband_entry.GT.is_het())
                        )

# homozygous and het parents
rec_trio_criteria = ((phased_sv_tm.trio_status=='trio') &  
                     (phased_sv_tm.proband_entry.GT.is_hom_var()) &
                     (phased_sv_tm.mother_entry.GT.is_het()) &
                     (phased_sv_tm.father_entry.GT.is_het())
                    )  
# hom and unaff are not hom
rec_non_trio_criteria = ((phased_sv_tm.trio_status!='trio') &  
                        (phased_sv_tm.n_cohort_hom_var_unaffected==0) & 
                        (phased_sv_tm.proband_entry.GT.is_hom_var())
                        )

phased_sv_tm = phased_sv_tm.annotate_entries(dominant_gt=((dom_trio_criteria) | (dom_non_trio_criteria)),
                              recessive_gt=((rec_trio_criteria) | (rec_non_trio_criteria)))

# NEW 2/19/2025: Annotate categories for merged output, instead of filtering
# NEW 2/20/2025: Rewrote variant_category annotation because using append was messing things up in Hail?
# NEW 5/29/2025: Added PREDICTED_LOF only criterion for P/LP output
size_field = [x for x in list(sv_mt.info) if 'passes_SVLEN_filter_' in x][0]
phased_sv_tm = phased_sv_tm.annotate_rows(
    variant_category = hl.array([
        # Category 1: Pathogenic only (P/LP)
        hl.if_else(
            ((hl.any(lambda x: x.matches('athogenic'), phased_sv_tm.info.clinical_interpretation)) |  # ClinVar P/LP
            (hl.is_defined(phased_sv_tm.info.dbvar_pathogenic))) &  # dbVar Pathogenic 
            (phased_sv_tm.info.PREDICTED_LOF.size()!=0),  # PREDICTED_LOF only
            'P/LP', 
            hl.missing(hl.tstr)
        ),
        
        # Category 2: GD only (GD region)
        hl.if_else(
            hl.is_defined(phased_sv_tm.info.gd_sv_name),  # GD region
            'GD', 
            hl.missing(hl.tstr)
        ),
        
        # Category 3: Large regions (passes_SVLEN_filter)
        # NEW 7/14/2025: restrict large regions to not BND or CPX, impacts a gene from restrictive fields
        hl.if_else(
            (phased_sv_tm.info[size_field]) &
            (~hl.array(['BND','CPX']).contains(phased_sv_tm.info.SVTYPE)) &
            (phased_sv_tm.info.restrictive_csq_genes.size()>0), 
            'large_region', 
            hl.missing(hl.tstr)
        ),
        
        # Category 4: OMIM AD and XLD (dominant_freq)
        hl.if_else(
            (hl.any(lambda x: x.matches('1') | x.matches('3'), phased_sv_tm.info.restrictive_inheritance_code)) & (phased_sv_tm.info.dominant_freq),
            'dominant', 
            hl.missing(hl.tstr)
        ),
        
        # Category 5: OMIM AR and XLR (recessive_freq)
        hl.if_else(
            (hl.any(lambda x: x.matches('2') | x.matches('4'), phased_sv_tm.info.restrictive_inheritance_code)) & (phased_sv_tm.info.recessive_freq),
            'recessive', 
            hl.missing(hl.tstr)
        ),

        # Category 6: Genic
        # NEW 7/14/2025: new 'genic' category for any restrictive fields' genes that are in a gene list
        hl.if_else(
            (phased_sv_tm.info.restrictive_gene_list.size()>0),
            'genic', 
            hl.missing(hl.tstr)
        )
    ]).filter(lambda x: hl.is_defined(x))  # Filter out null values
)

# NEW 7/14/2025: restrict large_region-only outputs to private to a family
phased_sv_tm = phased_sv_tm.filter_entries(
    (
        # large_region-only variants that are private to a family
        (
            (phased_sv_tm.variant_category.size() == 1) &
            (phased_sv_tm.variant_category[0] == 'large_region')
        ) &
        (
            (phased_sv_tm.n_cohort_het_unaffected == phased_sv_tm.n_family_het_unaffected) &
            (phased_sv_tm.n_cohort_hom_var_unaffected == phased_sv_tm.n_family_hom_var_unaffected) &
            (phased_sv_tm.n_cohort_het_affected == phased_sv_tm.n_family_het_affected) &
            (phased_sv_tm.n_cohort_hom_var_affected == phased_sv_tm.n_family_hom_var_affected)
        )
    ) |
    # variants that are not exclusively large_region
    (
        hl.set(['large_region']).intersection(hl.set(phased_sv_tm.variant_category)).size() <
        hl.set(phased_sv_tm.variant_category).size()
    )
)

merged_output_tm = phased_sv_tm.filter_rows(phased_sv_tm.variant_category.size()>0)
merged_output_tm = filter_and_annotate_tm(merged_output_tm)

# export merged TSV
merged_output_tm.entries().flatten().export(os.path.basename(sv_vcf).split('.vcf')[0] + '_merged_variants.tsv.gz', delimiter='\t')