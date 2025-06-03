###
# Copied over from hail_filter_comphets_xlr_hom_var_NIFS_v0.1.py on 1/27/2025
# to incorporate updates made to NIFS pipeline

## CHANGE LOG:
'''
1/27/2025:
- removed all code related to CA (NIFS-specific)
- edited get_non_trio_comphets to remove NIFS-specific criteria
- added back filter_entries before aggregation (non-NIFS-specific) in phase_by_transmission_aggregate_by_gene
- edited mat_carrier filtering to use mother_entry.GT instead of CA for mother het status
1/28/2025:
- allow for missing proband_entry.AD (e.g. for SVs)
- remove OMIM_MIM_number as SNV/Indel annotation
- variant_source annotation if not including ClinVar
- dummy variant_source annotation for SVs (TODO: edit with OMIM/ClinVar equivalent after updating SVs!)
1/30/2025:
- added comments to get_trio_comphets, added variant_type and variant_source annotations to trio_phased_tm
- allow for missing proband_entry.AD (e.g. for SVs) in get_subset_tm
- removed omim_uri, sv_gene_fields, rec_gene_list_tsv (annotation task added to filterClinicalVariantsSV_v0.1.wdl)
- flatten INFO and vep.transcript_consequences fields
- filter out rows where CHR2 is not the same chromosome (SV spans multiple chromosomes)
1/31/2025:
- added remove_parent_probands_trio_matrix function --> removes redundant "trios"
2/3/2025:
- annotate phenotype and unaffected/affected counts
- added annotate_and_filter_trio_matrix function to match SV outputs and reduce redundancy in code (run after phasing in phase_by_transmission_aggregate_by_gene)
- added get_transmission in phase_by_transmission_aggregate_by_gene
- added get_mendel_errors in annotate_and_filter_trio_matrix
- changed variant_type annotation to variant_types to retain original variant_type field for comphets
3/4/2025:
- change OMIM_recessive/OMIM_dominant to just recessive/dominant
- remove redundant gene field from output
5/14/2025: 
- don't drop original info/vep fields
- drop renamed INFO and VEP fields and retain originals (flattened later) to match other outputs
- add sex annotation to annotate_and_filter_trio_matrix
6/2/2025:
- use restrictive CSQ fields for SV gene fields
6/3/2025:
- rename all INFO struct fields to "info.{field}" and VEP struct fields to "vep.transcript_consequences.{field}"
- annotate inheritance_code field from vep.transcript_consequences
- use restrictive_csq_genes as "gene" field for comphets, restrictive_inheritance_code as "inheritance_code" field
- don't drop renamed INFO and VEP fields because already renamed above
- keep 'gene' column
'''
###

from pyspark.sql import SparkSession
from clinical_helper_functions import filter_mt, remove_parent_probands_trio_matrix, load_split_vep_consequences, require_biallelic, mendel_errors, annotate_trio_matrix
import hail as hl
import numpy as np
import pandas as pd
import sys
import ast
import os

from typing import Tuple

import hail.expr.aggregators as agg
from hail.expr import expr_call, expr_float64
from hail.genetics.pedigree import Pedigree
from hail.matrixtable import MatrixTable
from hail.table import Table
from hail.typecheck import numeric, typecheck
from hail.utils.java import Env

snv_indel_vcf = sys.argv[1]
clinvar_vcf = sys.argv[2]
sv_vcf = sys.argv[3]
ped_uri = sys.argv[4]
prefix = sys.argv[5]
build = sys.argv[6]
cores = sys.argv[7]  # string
mem = int(np.floor(float(sys.argv[8])))
ad_alt_threshold = int(sys.argv[9])
carrier_gene_list = sys.argv[10]
                               
hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": "true",
                    "spark.driver.extraJavaOptions": "-XX:ReservedCodeCacheSize=512m",
                    "spark.executor.extraJavaOptions": "-XX:ReservedCodeCacheSize=512m"
                    },
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

## STEP 1: Merge SNV/Indel VCF with SV VCF (or just one of them)
# Load SNV/Indel VCF
if snv_indel_vcf!='NA':
    locus_expr = 'locus'
    snv_mt = load_split_vep_consequences(snv_indel_vcf, build)

    # Load and merge SNV/Indel ClinVar P/LP VCF
    if clinvar_vcf!='NA':
        clinvar_mt = load_split_vep_consequences(clinvar_vcf, build) 
        # NEW 1/14/2025: added variant_source —— ClinVar_P/LP or recessive or both
        snv_mt_no_clinvar = snv_mt
        snv_mt = snv_mt.union_rows(clinvar_mt).distinct_by_row()
        snv_mt = snv_mt.annotate_rows(variant_source=hl.if_else(hl.is_defined(clinvar_mt.rows()[snv_mt.row_key]),  # if in ClinVar
                                                 hl.if_else(hl.is_defined(snv_mt_no_clinvar.rows()[snv_mt.row_key]), 'ClinVar_P/LP_recessive',  # if also in recessive_vcf
                                                            'ClinVar_P/LP'), 'recessive'))
    # NEW 1/28/2025: variant_source annotation if not including ClinVar
    else:
        snv_mt = snv_mt.annotate_rows(variant_source='recessive')

    # Explode rows by transcript
    snv_mt = snv_mt.explode_rows(snv_mt.vep.transcript_consequences)

    # NEW 1/9/2025: annotate gnomad_popmax_af after exploding by transcript
    # NEW 1/15/2025: commented out because now annotated (in INFO) in hail_filter_clinical_variants_v0.1.py :)
    # gnomad_fields = [x for x in list(snv_mt.vep.transcript_consequences) if 'gnomAD' in x]
    # snv_mt = snv_mt.annotate_rows(
    #     gnomad_popmax_af=hl.max([hl.or_missing(snv_mt.vep.transcript_consequences[gnomad_field]!='',
    #                                     hl.float(snv_mt.vep.transcript_consequences[gnomad_field])) 
    #                             for gnomad_field in gnomad_fields]))
    
    # Filter SNV/Indel MT
    snv_mt = filter_mt(snv_mt)

    # Filter out empty gene fields
    # NEW 6/3/2025: annotate inheritance_code field from vep.transcript_consequences
    snv_mt = snv_mt.annotate_rows(gene=snv_mt['vep']['transcript_consequences']['SYMBOL'],
                                 inheritance_code=snv_mt['vep']['transcript_consequences']['inheritance_code'])
    snv_mt = snv_mt.filter_rows(snv_mt.gene!='')

    snv_mt = snv_mt.annotate_rows(variant_type='SNV/Indel', 
                                  gene_source=['vep'])
    
    # NEW 1/30/2025: flatten INFO and vep.transcript_consequences fields
    snv_info_fields = list(snv_mt.info)
    vep_fields = list(snv_mt.vep.transcript_consequences)
            
    # NEW 6/3/2025: rename all INFO struct fields to "info.{field}" and VEP struct fields to "vep.transcript_consequences.{field}"
    new_info_field_map = {og_field: f"info.{og_field}" for og_field in snv_info_fields}
    new_vep_field_map = {og_field: f"vep.transcript_consequences.{og_field}" for og_field in vep_fields}

    snv_mt = snv_mt.annotate_rows(**{new_field: snv_mt.info[og_field] for og_field, new_field in new_info_field_map.items()} | 
                    {new_field: snv_mt.vep.transcript_consequences[og_field] for og_field, new_field in new_vep_field_map.items()})
    
    snv_mt = snv_mt.drop('info','vep')

# Load SV VCF
if sv_vcf!='NA':
    locus_expr = 'locus_interval'
    sv_mt = hl.import_vcf(sv_vcf, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)

    # filter out BNDs
    sv_mt = sv_mt.filter_rows(sv_mt.info.SVTYPE!='BND') 
    # NEW 1/30/2025: filter out rows where CHR2 is not the same chromosome (SV spans multiple chromosomes)
    sv_mt = sv_mt.filter_rows(sv_mt.info.CHR2==sv_mt.locus.contig)
    
    sv_mt = sv_mt.annotate_rows(variant_type='SV')

    # NEW 1/30/2025: flatten INFO fields
    sv_info_fields = list(sv_mt.info)

    # NEW 6/3/2025: rename all INFO struct fields to "info.{field}"
    new_info_field_map = {og_field: f"info.{og_field}" for og_field in sv_info_fields}

    sv_mt = sv_mt.annotate_rows(**{new_field: sv_mt.info[og_field] for og_field, new_field in new_info_field_map.items()})    
    
    # NEW 1/30/2025: combine gene-level annotations in INFO where there is a value for each gene
    # Annotate gene to match SNV/Indels (to explode on and keep original genes annotation)
    # NEW 6/3/2025: use restrictive_csq_genes as "gene" field for comphets, restrictive_inheritance_code as "inheritance_code" field
    sv_gene_fields = ['gene','inheritance_code']
    sv_mt = sv_mt.annotate_rows(**{'gene': sv_mt['info.restrictive_csq_genes'],
                               'inheritance_code': sv_mt['info.restrictive_inheritance_code']})

    sv_mt = sv_mt.annotate_rows(gene_level=hl.zip(*[sv_mt[field] for field in sv_gene_fields])\
            .map(lambda x: hl.struct(**{field: x[i] 
                                        for i, field in enumerate(sv_gene_fields)})))\
        .explode_rows('gene_level')
    sv_mt = sv_mt.annotate_rows(**{field: sv_mt.gene_level[field] 
                                   for field in sv_gene_fields}).drop('gene_level')
    
    # NEW 1/28/2025: dummy variant_source annotation for SVs
    sv_mt = sv_mt.annotate_rows(variant_source='SV')

    # OMIM recessive code
    omim_rec_code = (sv_mt['inheritance_code'].matches('2'))
    # OMIM XLR code
    omim_xlr_code = (sv_mt['inheritance_code'].matches('4'))
    sv_mt = sv_mt.filter_rows(omim_rec_code | omim_xlr_code)
    
    sv_mt = sv_mt.drop('info')

# Unify SNV/Indel MT and SV MT row and entry fields
# NEW 1/30/2025: adjusted for INFO field flattened above
if (snv_indel_vcf!='NA') and (sv_vcf!='NA'):
    sv_row_fields, sv_entry_fields = list(sv_mt.row), list(sv_mt.entry)
    snv_row_fields, snv_entry_fields = list(snv_mt.row), list(snv_mt.entry)

    sv_missing_entry_fields = {field: str(snv_mt[field].dtype) for field in np.setdiff1d(snv_entry_fields, sv_entry_fields)}
    snv_missing_entry_fields = {field: str(sv_mt[field].dtype) for field in np.setdiff1d(sv_entry_fields, snv_entry_fields)}

    sv_missing_row_fields = {field: str(snv_mt[field].dtype) for field in np.setdiff1d(snv_row_fields, sv_row_fields)}
    snv_missing_row_fields = {field: str(sv_mt[field].dtype) for field in np.setdiff1d(sv_row_fields, snv_row_fields)}

    sv_mt = sv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in sv_missing_entry_fields.items()})
    snv_mt = snv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in snv_missing_entry_fields.items()})

    sv_mt = sv_mt.select_entries(*sorted(list(sv_mt.entry)))
    snv_mt = snv_mt.select_entries(*sorted(list(snv_mt.entry)))

    sv_mt = sv_mt.annotate_rows(**{field: hl.missing(dtype) for field, dtype in sv_missing_row_fields.items()})
    snv_mt = snv_mt.annotate_rows(**{field: hl.missing(dtype) for field, dtype in snv_missing_row_fields.items()})

    sv_mt = sv_mt.key_rows_by().select_rows(*sorted(list(sv_mt.row))).key_rows_by('locus','alleles')
    snv_mt = snv_mt.key_rows_by().select_rows(*sorted(list(snv_mt.row))).key_rows_by('locus','alleles')

    # Subset shared samples 
    sv_samps = sv_mt.s.collect()
    snv_samps = snv_mt.s.collect()
    shared_samps = list(np.intersect1d(sv_samps, snv_samps))

    if len(shared_samps)==0:
        shared_samps = ['']

    # Match column order before merging
    def align_mt2_cols_to_mt1(mt1, mt2):
        mt1 = mt1.add_col_index()
        mt2 = mt2.add_col_index()
        new_col_order = mt2.index_cols(mt1.col_key).col_idx.collect()
        return mt2.choose_cols(new_col_order)
    
    sv_mt = sv_mt.filter_cols(hl.array(shared_samps).contains(sv_mt.s))
    snv_mt = align_mt2_cols_to_mt1(sv_mt, snv_mt)

    variant_types = 'SV_SNV_Indel'
    merged_mt = sv_mt.union_rows(snv_mt)
        
elif snv_indel_vcf!='NA':
    variant_types = 'SNV_Indel'
    merged_mt = snv_mt

elif sv_vcf!='NA':
    variant_types = 'SV'
    merged_mt = sv_mt

# Clean up merged SV VCF with SNV/Indel VCF
# NEW 1/30/2025: adjusted for INFO field flattened above
# NEW 6/3/2025: adjust INFO fields for 'info.' prefix
if sv_vcf!='NA':
    # Change locus to locus_interval to include END for SVs
    merged_mt = merged_mt.annotate_rows(end=hl.if_else(hl.is_defined(merged_mt['info.END2']), merged_mt['info.END2'], merged_mt['info.END']))
    # Account for INS having same END but different SVLEN
    merged_mt = merged_mt.annotate_rows(end=hl.if_else(merged_mt['info.SVTYPE']=='INS', merged_mt.end + merged_mt['info.SVLEN'], merged_mt.end))
    # Account for empty END for SNV/Indels
    merged_mt = merged_mt.annotate_rows(end=hl.if_else(merged_mt.variant_type=='SNV/Indel', merged_mt.locus.position, merged_mt.end))
    merged_mt = merged_mt.key_rows_by()
    merged_mt = merged_mt.annotate_rows(locus_interval=hl.locus_interval(contig=merged_mt.locus.contig, 
                                                                              start=merged_mt.locus.position,
                                                                              end=merged_mt.end, reference_genome=build))
    merged_mt = merged_mt.key_rows_by(locus_expr, 'alleles')

# NEW 1/14/2025: Annotate PAR status (moved up from end of script)
merged_mt = merged_mt.annotate_rows(in_non_par=~(merged_mt.locus.in_autosome_or_par()))

# NEW 2/3/2025: Annotate phenotype and unaffected/affected counts
# Annotate affected status/phenotype from pedigree
ped_ht = hl.import_table(ped_uri, delimiter='\t').key_by('sample_id')
merged_mt = merged_mt.annotate_cols(phenotype=ped_ht[merged_mt.s].phenotype)

# Get cohort unaffected/affected het and homvar counts
merged_mt = merged_mt.annotate_rows(**{
    "n_het_unaffected": hl.agg.filter(merged_mt.phenotype=='1', hl.agg.sum(merged_mt.GT.is_het())),
    "n_hom_var_unaffected": hl.agg.filter(merged_mt.phenotype=='1', hl.agg.sum(merged_mt.GT.is_hom_var())),
    "n_het_affected": hl.agg.filter(merged_mt.phenotype=='2', hl.agg.sum(merged_mt.GT.is_het())),
    "n_hom_var_affected": hl.agg.filter(merged_mt.phenotype=='2', hl.agg.sum(merged_mt.GT.is_hom_var()))
})

## STEP 2: Get CompHets
# Mendel errors
def get_mendel_errors(mt, phased_tm, pedigree):
    all_errors, per_fam, per_sample, per_variant = mendel_errors(mt['GT'], pedigree)  # edited Hail function, see above
    all_errors_mt = all_errors.key_by().to_matrix_table(row_key=[locus_expr,'alleles'], col_key=['s'], col_fields=['fam_id'])
    phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)
    return phased_tm

def get_transmission(phased_tm):
    phased_tm = phased_tm.annotate_entries(transmission=hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|0'), 'uninherited',
            hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('0|1'), 'inherited_from_mother',
                        hl.if_else(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|0'), 'inherited_from_father',
                                hl.or_missing(phased_tm.proband_entry.PBT_GT==hl.parse_call('1|1'), 'inherited_from_both'))))
    )
    return phased_tm

def phase_by_transmission_aggregate_by_gene(tm, mt, pedigree):
    # filter out calls that are hom ref in proband
    # NEw 1/27/2025: added back filter_entries before aggregation (non-NIFS-specific)
    tm = tm.filter_entries(tm.proband_entry.GT.is_non_ref())

    phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')
    
    # NEW 2/3/2025: Run annotate_and_filter_trio_matrix after phasing in phase_by_transmission_aggregate_by_gene
    phased_tm = annotate_and_filter_trio_matrix(phased_tm, mt, pedigree, ped_ht, locus_expr) 
    
    phased_tm = phased_tm.key_rows_by(locus_expr,'alleles','gene')

    # NEW 2/3/2025: get_transmission in phase_by_transmission_aggregate_by_gene
    phased_tm = get_transmission(phased_tm)

    # NEW 1/14/2025: filter_entries before aggregation
    phased_tm = phased_tm.filter_entries(hl.is_defined(phased_tm.proband_entry.GT))
    # NEW 1/21/2025: filter by proband GT before aggregating rows
    phased_tm = phased_tm.filter_rows(hl.agg.count_where(hl.is_defined(phased_tm.proband_entry.GT))>0)
    gene_agg_phased_tm = (phased_tm.group_rows_by(phased_tm.gene)
        .aggregate_rows(locus_alleles = hl.agg.collect(phased_tm.row_key),
                       variant_type = hl.agg.collect_as_set(phased_tm.variant_type),  # SET
                       variant_source = hl.agg.collect_as_set(phased_tm.variant_source)  # SET
                       )  # NEW 1/14/2025: added variant_source
        .aggregate_entries(all_locus_alleles=hl.agg.collect(phased_tm.row_key),
                          proband_PBT_GT = hl.agg.collect(phased_tm.proband_entry.PBT_GT),
                          proband_GT = hl.agg.collect(phased_tm.proband_entry.GT)
        )).result()
    
    return phased_tm, gene_agg_phased_tm

def annotate_and_filter_trio_matrix(tm, mt, pedigree, ped_ht, locus_expr):
    tm = annotate_trio_matrix(tm, mt, pedigree, ped_ht, locus_expr=locus_expr) 
    # filter by AD of alternate allele in proband
    # NEW 1/30/2025: allow for missing proband_entry.AD (e.g. for SVs)
    if 'AD' in list(mt.entry):
        tm = tm.filter_entries((tm.proband_entry.AD[1]>=ad_alt_threshold) |
                                                    (hl.is_missing(tm.proband_entry.AD)))
    return tm    

def get_subset_tm(mt, samples, pedigree, keep=True, complete_trios=False):
    subset_mt = mt.filter_cols(hl.array(samples).contains(mt.s), keep=keep)

    # remove variants missing in subset samples
    subset_mt = hl.variant_qc(subset_mt)
    subset_mt = subset_mt.filter_rows(subset_mt.variant_qc.AC[1]>0)
    subset_mt = subset_mt.drop('variant_qc')

    subset_tm = hl.trio_matrix(subset_mt, pedigree, complete_trios=complete_trios)
    subset_tm = remove_parent_probands_trio_matrix(subset_tm)  # NEW 1/31/2025: Removes redundant "trios"  
    return subset_mt, subset_tm

def get_non_trio_comphets(mt):
    non_trio_mt, non_trio_tm = get_subset_tm(mt, non_trio_samples, non_trio_pedigree)
    non_trio_phased_tm, non_trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(non_trio_tm, non_trio_mt, non_trio_pedigree)

    # Filter to genes where at least one sample has multiple variants
    potential_comp_hets_non_trios = non_trio_gene_agg_phased_tm.filter_rows(
            hl.agg.count_where(non_trio_gene_agg_phased_tm.proband_GT.size()>1)>0
    )
    # Explode by variant (locus_expr, alleles) --> key by locus_expr, alleles, gene
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.explode_rows(potential_comp_hets_non_trios.locus_alleles)
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.key_rows_by(potential_comp_hets_non_trios.locus_alleles[locus_expr], potential_comp_hets_non_trios.locus_alleles.alleles, 'gene')

    # Filter to variants that hit genes with multiple variants
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.filter_entries(potential_comp_hets_non_trios.proband_GT.size()>1)
   
    # Annotate non-gene-aggregated TM using potential comphets from gene-aggregated TM 
    non_trio_phased_tm = non_trio_phased_tm.key_rows_by(locus_expr, 'alleles', 'gene')
    non_trio_phased_tm = non_trio_phased_tm.annotate_entries(locus_alleles=  
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].all_locus_alleles,
                                                            proband_GT=
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].proband_GT,
                                                            proband_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].proband_GT),
                                                            proband_PBT_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].proband_PBT_GT),
    )

    # Filter non-gene-aggregated TM to variants that hit genes with multiple unique variants
    non_trio_phased_tm = non_trio_phased_tm.filter_entries((hl.set(non_trio_phased_tm.locus_alleles).size()>1) &
                                                                     (non_trio_phased_tm.proband_GT.size()>1))  
    
    # NEW 1/14/2025: Annotate variant_type and variant_source as comma-separated strings of unique values (basically per gene)
    # NEW 2/3/2025: Change variant_type annotation to variant_types to retain original variant_type field for comphets
    non_trio_phased_tm = non_trio_phased_tm.annotate_rows(variant_types=  
        hl.str(', ').join(hl.sorted(hl.array(potential_comp_hets_non_trios.rows()[non_trio_phased_tm.row_key].variant_type))),
                                                        variant_source= 
        hl.str(', ').join(hl.sorted(hl.array(potential_comp_hets_non_trios.rows()[non_trio_phased_tm.row_key].variant_source))),
    )
    # Grab rows (variants) from non-gene-aggregated TM, of potential comphets from gene-aggregated TM
    phased_tm_comp_hets_non_trios = non_trio_phased_tm.semi_join_rows(potential_comp_hets_non_trios.rows()).key_rows_by(locus_expr, 'alleles')
    return phased_tm_comp_hets_non_trios.drop('locus_alleles')  

def get_trio_comphets(mt):
    trio_mt, trio_tm = get_subset_tm(mt, trio_samples, trio_pedigree, keep=True, complete_trios=True)
    trio_phased_tm, trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(trio_tm, trio_mt, trio_pedigree)

    # Filter to genes where at least one sample has multiple *phased* variants
    potential_comp_hets_trios = trio_gene_agg_phased_tm.filter_rows(
        hl.agg.count_where(hl.set(trio_gene_agg_phased_tm.proband_PBT_GT).size()>1)>0
    )
    # Explode by variant (locus_expr, alleles) --> key by locus_expr, alleles, gene
    potential_comp_hets_trios = potential_comp_hets_trios.explode_rows(potential_comp_hets_trios.locus_alleles)
    potential_comp_hets_trios = potential_comp_hets_trios.key_rows_by(potential_comp_hets_trios.locus_alleles[locus_expr], potential_comp_hets_trios.locus_alleles.alleles, 'gene')

    # Filter to variants that hit genes with multiple *phased* variants
    potential_comp_hets_trios = potential_comp_hets_trios.filter_entries(hl.set(potential_comp_hets_trios.proband_PBT_GT).size()>1)

    # Annotate non-gene-aggregated TM using potential comphets from gene-aggregated TM 
    trio_phased_tm = trio_phased_tm.key_rows_by(locus_expr, 'alleles', 'gene')
    trio_phased_tm = trio_phased_tm.annotate_entries(locus_alleles=  
        potential_comp_hets_trios[trio_phased_tm.row_key, trio_phased_tm.col_key].all_locus_alleles,
                                                            proband_GT=
        potential_comp_hets_trios[trio_phased_tm.row_key, trio_phased_tm.col_key].proband_GT,
                                                            proband_GT_set=hl.set(
        potential_comp_hets_trios[trio_phased_tm.row_key, trio_phased_tm.col_key].proband_GT),
                                                            proband_PBT_GT_set=hl.set(
        potential_comp_hets_trios[trio_phased_tm.row_key, trio_phased_tm.col_key].proband_PBT_GT),
    )
    
    # Filter non-gene-aggregated TM to variants that hit genes with multiple unique variants *that are in trans*
    trio_phased_tm = trio_phased_tm.filter_entries(trio_phased_tm.proband_PBT_GT_set.size()>1)  

    # NEW 1/30/2025: (from get_trio_comphets) Annotate variant_type and variant_source as comma-separated strings of unique values (basically per gene)
    # NEW 2/3/2025: Change variant_type annotation to variant_types to retain original variant_type field for comphets
    trio_phased_tm = trio_phased_tm.annotate_rows(variant_types=  
        hl.str(', ').join(hl.sorted(hl.array(potential_comp_hets_trios.rows()[trio_phased_tm.row_key].variant_type))),
                                                        variant_source= 
        hl.str(', ').join(hl.sorted(hl.array(potential_comp_hets_trios.rows()[trio_phased_tm.row_key].variant_source))),
    )
    # Grab rows (variants) from non-gene-aggregated TM, of potential comphets from gene-aggregated TM    
    phased_tm_comp_hets_trios = trio_phased_tm.semi_join_rows(potential_comp_hets_trios.rows()).key_rows_by(locus_expr, 'alleles')
    return phased_tm_comp_hets_trios.drop('locus_alleles')  

# Subset pedigree to samples in VCF, edit parental IDs
vcf_samples = merged_mt.s.collect()
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped = tmp_ped[tmp_ped.sample_id.isin(vcf_samples)].copy()
tmp_ped['paternal_id'] = tmp_ped.paternal_id.apply(lambda id: id if id in vcf_samples else '0')
tmp_ped['maternal_id'] = tmp_ped.maternal_id.apply(lambda id: id if id in vcf_samples else '0')
tmp_ped.to_csv(f"{prefix}.ped", sep='\t', index=False)

pedigree = hl.Pedigree.read(f"{prefix}.ped", delimiter='\t')
pedigree = pedigree.filter_to(vcf_samples)

trio_samples = list(np.intersect1d(vcf_samples,
                              list(np.array([[trio.s, trio.pat_id, trio.mat_id] 
                                             for trio in pedigree.complete_trios() if trio.fam_id!='-9']).flatten())))
non_trio_samples = list(np.setdiff1d(vcf_samples, trio_samples))

trio_pedigree = pedigree.filter_to(trio_samples)
non_trio_pedigree = pedigree.filter_to(non_trio_samples)

## Get CompHets 
# Filter to only in autosomes or PAR
comphet_mt = merged_mt.filter_rows(merged_mt.locus.in_autosome_or_par())

if len(trio_samples)>0:
    merged_trio_comphets = get_trio_comphets(comphet_mt)
    merged_comphets = merged_trio_comphets.entries()

if len(non_trio_samples)>0:
    merged_non_trio_comphets = get_non_trio_comphets(comphet_mt)
    merged_comphets = merged_non_trio_comphets.entries()

if (len(trio_samples)>0) and (len(non_trio_samples)>0):
    merged_comphets = merged_trio_comphets.entries().union(merged_non_trio_comphets.entries())

if len(trio_samples)==0:
    trio_samples = ['']

# Trio matrix
merged_tm = hl.trio_matrix(merged_mt, pedigree, complete_trios=False)
merged_tm = remove_parent_probands_trio_matrix(merged_tm)  # NEW 1/31/2025: Removes redundant "trios"  

gene_phased_tm, gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(merged_tm, merged_mt, pedigree)

# NEW 1/13/2025: maternal carrier variants
# NEW 1/30/2025: edited gene_phased_tm.vep.transcript_consequences.SYMBOL --> gene_phased_tm.gene,
# in carrier gene list and mother is het
carrier_genes = pd.read_csv(carrier_gene_list, sep='\t', header=None)[0].tolist()
mat_carrier = gene_phased_tm.filter_rows(hl.array(carrier_genes).contains(gene_phased_tm.gene))
mat_carrier = mat_carrier.filter_entries(mat_carrier.mother_entry.GT.is_het()).key_rows_by(locus_expr, 'alleles').entries()

# XLR only
xlr_phased_tm = gene_phased_tm.filter_rows(gene_phased_tm['inheritance_code'].matches('4'))   # OMIM XLR
xlr_phased = xlr_phased_tm.filter_entries((xlr_phased_tm.proband_entry.GT.is_non_ref()) &
                            (~xlr_phased_tm.is_female)).key_rows_by(locus_expr, 'alleles').entries()

# HomVar in proband only
phased_hom_var = gene_phased_tm.filter_entries(gene_phased_tm.proband_entry.GT.is_hom_var())
phased_hom_var = phased_hom_var.filter_entries((phased_hom_var.locus.in_x_nonpar()) &
                            (~phased_hom_var.is_female), keep=False)  # filter out non-PAR chrX in males
phased_hom_var = phased_hom_var.filter_rows(hl.agg.count_where(
    hl.is_defined(phased_hom_var.proband_entry.GT))>0).key_rows_by(locus_expr, 'alleles').entries()

xlr_phased = xlr_phased.annotate(variant_category='XLR')
phased_hom_var = phased_hom_var.annotate(variant_category='hom_var')
merged_comphets = merged_comphets.annotate(variant_category='comphet')
mat_carrier = mat_carrier.annotate(variant_category='maternal_carrier')

# NEW 5/14/2025: drop renamed INFO and VEP fields and retain originals (flattened later) to match other outputs
# NEW 6/3/2025: don't drop renamed INFO and VEP fields because already renamed above
# merged_comphets = merged_comphets.drop(*(list(new_vep_field_map.values()) + list(new_info_field_map.values())))
# xlr_phased = xlr_phased.drop(*(list(new_vep_field_map.values()) + list(new_info_field_map.values())))
# phased_hom_var = phased_hom_var.drop(*(list(new_vep_field_map.values()) + list(new_info_field_map.values())))
# mat_carrier = mat_carrier.drop(*(list(new_vep_field_map.values()) + list(new_info_field_map.values())))

# NEW 1/14/2025: use to_pandas() to bypass ClassTooLargeException in Hail tables union
# NEW 1/22/2025: use export() and then load in pandas instead of to_pandas() to match formatting with other outputs
# NEW 3/4/2025: remove redundant gene field from output
# NEW 6/3/2025: keep 'gene' column
merged_comphets.drop('proband_GT','proband_GT_set','proband_PBT_GT_set').flatten().export('comphets.tsv.gz')
xlr_phased.flatten().export('xlr.tsv.gz')
phased_hom_var.flatten().export('hom_var.tsv.gz')
mat_carrier.flatten().export('mat_carrier.tsv.gz')

merged_comphets_df = pd.read_csv('comphets.tsv.gz', sep='\t')
xlr_phased_df = pd.read_csv('xlr.tsv.gz', sep='\t')
phased_hom_var_df = pd.read_csv('hom_var.tsv.gz', sep='\t')
mat_carrier_df = pd.read_csv('mat_carrier.tsv.gz', sep='\t')

merged_comphets_xlr_hom_var_mat_carrier_df = pd.concat([merged_comphets_df, xlr_phased_df, phased_hom_var_df, mat_carrier_df])

output_filename = f"{prefix}_{variant_types}_comp_hets_xlr_hom_var_mat_carrier.tsv.gz"
if len(output_filename)>(os.pathconf('/', 'PC_NAME_MAX')-len('/cromwell_root/.')):  # if filename too long
    output_filename = f"{variant_types}_comp_hets_xlr_hom_var_mat_carrier.tsv.gz"

merged_comphets_xlr_hom_var_mat_carrier_df.to_csv(output_filename, sep='\t', index=False)