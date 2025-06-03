###
# New script based on some code from hail_filter_comphets_xlr_hom_var_v0.1.py
# and filterClinicalVariantsSV_v0.1.wdl, tested in test_sv_annotations.ipynb. 
# Created 1/29/2025.

## CHANGE LOG:
'''
1/30/2025:
- added gnomAD_popmax_AF to INFO and gnomad_popmax_freq flag
- changed gene_lists annotation to gene_list in INFO to match SNV/Indels
2/19/2025:
- remove sv_gene_fields input and change Python variable to be a union of restrictive_csq_fields and permissive_csq_fields
2/20/2025:
- allow for missing gnomAD AFs
5/29/2025:
- more AC and AF cutoffs, PED inputs
- annotate cohort unaffected/affected counts and AC (moved from hail_filter_clinical_sv_v0.1.py to hail_annotate_sv_gene_level_v0.1.py)
6/2/2025:
- rename permissive_csq, restrictive_csq to permissive_csq_fields, restrictive_csq_fields
- add permissive_gene_source, restrictive_gene_source, permissive_inheritance_code, restrictive_inheritance_code
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
parser.add_argument('-o', dest='output_filename', help='Output filename')
parser.add_argument('-l', dest='gene_list_tsv', help='Gene list tsv for annotations')
parser.add_argument('-s', dest='size_threshold', help='Size threshold in BP for SVLEN flag')
parser.add_argument('--inheritance', dest='inheritance_uri', help='inheritance file for annotations')
parser.add_argument('--cores', dest='cores', help='CPU cores')
parser.add_argument('--mem', dest='mem', help='Memory')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--permissive-csq-fields', dest='permissive_csq_fields', help='PREDICTED_* fields to consider for permissive_csq_genes field in INFO')
parser.add_argument('--restrictive-csq-fields', dest='restrictive_csq_fields', help='PREDICTED_* fields to consider for restrictive_csq_genes field in INFO')
parser.add_argument('--constrained-uri', dest='constrained_uri', help='File for constrained genes')
parser.add_argument('--prec-uri', dest='prec_uri', help='File for pRec genes')
parser.add_argument('--hi-uri', dest='hi_uri', help='File for HI genes')
parser.add_argument('--ts-uri', dest='ts_uri', help='File for TS genes')
parser.add_argument('--ped', dest='ped_uri', help='Input ped file')
parser.add_argument('--dom-af', dest='dom_af_threshold', help='Cohort AF threshold for dominants')
parser.add_argument('--rec-af', dest='rec_af_threshold', help='Cohort AF threshold for recessives')
parser.add_argument('--gnomad-dom-af', dest='gnomad_af_dom_threshold', help='gnomAD AF threshold for dominants')
parser.add_argument('--gnomad-rec-af', dest='gnomad_af_rec_threshold', help='gnomAD AF threshold for recessives')
parser.add_argument('--gnomad-af-field', dest='gnomad_af_field', help='Field for gnomAD AFs in INFO')
parser.add_argument('--gnomad-popmax-af', dest='gnomad_popmax_af_threshold', help='gnomAD popmax AF threshold')
parser.add_argument('--rec-n-hom-var', dest='rec_n_cohort_hom_var_threshold', help='Number of hom var individuals threshold for recessives')
parser.add_argument('--dom-ac', dest='dom_ac_threshold', help='Cohort AC threshold for dominants')
parser.add_argument('--dom-ac-unaffected', dest='dom_ac_unaffected_threshold', help='AC of unaffected individuals threshold for dominants')

args = parser.parse_args()

sv_vcf = args.vcf_file
output_filename = args.output_filename
cores = args.cores  # string
mem = int(np.floor(float(args.mem)))
genome_build = args.build
gene_list_tsv = args.gene_list_tsv
inheritance_uri = args.inheritance_uri
size_threshold = int(args.size_threshold)
restrictive_csq_fields = (args.restrictive_csq_fields).split(',')
permissive_csq_fields = (args.permissive_csq_fields).split(',')
constrained_uri = args.constrained_uri
ped_uri = args.ped_uri
prec_uri = args.prec_uri
hi_uri = args.hi_uri
ts_uri = args.ts_uri
dom_af_threshold = float(args.dom_af_threshold)
rec_af_threshold = float(args.rec_af_threshold)
gnomad_af_dom_threshold = float(args.gnomad_af_dom_threshold)
gnomad_af_rec_threshold = float(args.gnomad_af_rec_threshold)
gnomad_af_field = args.gnomad_af_field
gnomad_popmax_af_threshold = float(args.gnomad_popmax_af_threshold)
# NEW 5/29/2025: More AC and AF cutoffs, PED inputs
rec_n_cohort_hom_var_threshold = int(args.rec_n_cohort_hom_var_threshold)
dom_ac_threshold = int(args.dom_ac_threshold)
dom_ac_unaffected_threshold =int(args.dom_ac_unaffected_threshold)

# NEW 2/19/2025: Remove sv_gene_fields input and change Python variable to be a union of restrictive_csq_fields and permissive_csq_fields
sv_gene_fields = list(np.union1d(permissive_csq_fields, restrictive_csq_fields))

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

# Load SV VCF
locus_expr = 'locus_interval'
sv_mt = hl.import_vcf(sv_vcf, reference_genome=genome_build, force_bgz=True, call_fields=[], array_elements_required=False)
header = hl.get_vcf_metadata(sv_vcf)

# NEW 5/29/2025: Moved from hail_filter_clinical_sv_v0.1.py to hail_annotate_sv_gene_level_v0.1.py
# Annotate affected status/phenotype from pedigree
# NEW 2/5/2025: Moved ped format standardization to before ped_ht
tmp_ped = pd.read_csv(ped_uri, sep='\t').iloc[:,:6]
tmp_ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
cropped_ped_uri = f"{os.path.basename(ped_uri).split('.ped')[0]}_crop.ped"
tmp_ped.to_csv(cropped_ped_uri, sep='\t', index=False)

ped_ht = hl.import_table(cropped_ped_uri, delimiter='\t').key_by('sample_id')
sv_mt = sv_mt.annotate_cols(phenotype=ped_ht[sv_mt.s].phenotype)

# Get cohort unaffected/affected het and homvar counts
sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(**{
    "n_cohort_het_unaffected": hl.agg.filter(sv_mt.phenotype=='1', hl.agg.sum(sv_mt.GT.is_het())),
    "n_cohort_hom_var_unaffected": hl.agg.filter(sv_mt.phenotype=='1', hl.agg.sum(sv_mt.GT.is_hom_var())),
    "n_cohort_het_affected": hl.agg.filter(sv_mt.phenotype=='2', hl.agg.sum(sv_mt.GT.is_het())),
    "n_cohort_hom_var_affected": hl.agg.filter(sv_mt.phenotype=='2', hl.agg.sum(sv_mt.GT.is_hom_var()))
    })
)

# NEW 5/29/2025: Annotate cohort unaffected/affected AC  
sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(**{
    "cohort_AC_unaffected": sv_mt.info.n_cohort_het_unaffected + 2*sv_mt.info.n_cohort_hom_var_unaffected,
    "cohort_AC_affected": sv_mt.info.n_cohort_het_affected + 2*sv_mt.info.n_cohort_hom_var_affected
    })
)

# Frequency flags
# Use "max" AC and AF because Hail parses as array with 1 element
# NEW 2/20/2025: Allow for missing gnomAD AFs
# NEW 5/29/2025: Incorporate cohort unaffected/affected counts and AC
sv_mt = sv_mt.annotate_rows(
    info=sv_mt.info.annotate(
        dominant_freq=(
            (
                (hl.max(sv_mt.info.AC) <= dom_ac_threshold) |
                (hl.max(sv_mt.info.AF) <= dom_af_threshold)
            ) &
            (
                (sv_mt.info[gnomad_af_field] <= gnomad_af_dom_threshold) |
                hl.is_missing(sv_mt.info[gnomad_af_field])
            ) &
            (sv_mt.info.cohort_AC_unaffected <= dom_ac_unaffected_threshold)
        ),
        recessive_freq=(
            (
                (
                    (sv_mt.info.n_cohort_hom_var_unaffected + sv_mt.info.n_cohort_hom_var_affected)
                    <= rec_n_cohort_hom_var_threshold
                ) |
                (hl.max(sv_mt.info.AF) <= rec_af_threshold)
            ) &
            (
                (sv_mt.info[gnomad_af_field] <= gnomad_af_rec_threshold) |
                hl.is_missing(sv_mt.info[gnomad_af_field])
            )
        )
    )
)

n_tot_samples = sv_mt.count_cols()

## ANNOTATIONS
# Annotate gene_source
def get_predicted_sources_expr(mt, sv_gene_fields, gene_expr):
    return hl.array(
        [hl.or_missing(hl.array(mt.info[col]).contains(gene_expr), col) for col in sv_gene_fields]
    ).filter(hl.is_defined)

# Function to get gene-level annotations (including gene_list)
def get_gene_level_annotations(mt, gene_field, gene_field_list, inheritance_ht, prefix, gene_lists=None):
    """
    Generalized function to annotate gene_source, inheritance_code, and gene_list.

    Args:
        mt: Hail MatrixTable
        gene_field: str - Field in mt.info to explode by (e.g. 'genes')
        gene_field_list: list - INFO fields to infer gene sources
        inheritance_ht: Hail Table - keyed by gene symbol
        prefix: str - Prefix for output field names (e.g. 'restrictive_', 'permissive_', '')
        gene_lists: dict (optional) - {gene_list_name: list of gene symbols}

    Returns:
        MatrixTable aggregated by rsid with prefixed fields.
    """
    gene_expr = getattr(mt.info, gene_field)
    gene_mt = mt.explode_rows(gene_expr)  # Exploding rows to work with individual genes

    # Now gene_expr is a scalar, so we can safely index the inheritance_ht table
    gene_mt = gene_mt.annotate_rows(
        **{f"{prefix}gene_source": get_predicted_sources_expr(gene_mt, gene_field_list, gene_mt.info[gene_field])},
        **{f"{prefix}inheritance_code": hl.or_missing(
            hl.is_defined(inheritance_ht[gene_mt.info[gene_field]]),  # Using exploded gene field
            inheritance_ht[gene_mt.info[gene_field]].inheritance_code)}
    )

    if gene_lists:
        gene_mt = gene_mt.annotate_rows(
            **{f"{prefix}gene_list": hl.array([
                hl.or_missing(hl.array(gene_set).contains(gene_mt.info[gene_field]), list_name)
                for list_name, gene_set in gene_lists.items()
            ]).filter(hl.is_defined)}
        )

        # Convert gene_list to &-joined string if defined
        gene_mt = gene_mt.annotate_rows(
            **{f"{prefix}gene_list": hl.or_missing(
                hl.len(gene_mt[f"{prefix}gene_list"]) > 0,
                hl.str("&").join(gene_mt[f"{prefix}gene_list"]))}
        )

    # Convert gene_source to &-joined string
    gene_mt = gene_mt.annotate_rows(
        **{f"{prefix}gene_source": hl.or_missing(
            hl.len(gene_mt[f"{prefix}gene_source"]) > 0,
            hl.str("&").join(gene_mt[f"{prefix}gene_source"]))}
    )

    # Define aggregation fields
    agg_fields = {
        f"{prefix}gene_source": hl.agg.collect(gene_mt[f"{prefix}gene_source"]),
        f"{prefix}inheritance_code": hl.agg.collect(gene_mt[f"{prefix}inheritance_code"])
    }
    if gene_lists:
        agg_fields[f"{prefix}gene_list"] = hl.agg.collect(gene_mt[f"{prefix}gene_list"])

    return gene_mt.group_rows_by(gene_mt.rsid).aggregate_rows(**agg_fields).result()

# Load gene list file(s) if available
if gene_list_tsv != 'NA':
    gene_list_uris = pd.read_csv(gene_list_tsv, sep='\t', header=None).set_index(0)[1].to_dict()
    gene_list_dict = {
        name: hl.literal(set(pd.read_csv(uri, header=None)[0].tolist()))
        for name, uri in gene_list_uris.items()
    }
else:
    gene_list_dict = None

# Load inheritance table
inheritance_ht = hl.import_table(inheritance_uri).key_by('approvedGeneSymbol')

# Annotate genes in INFO for sv_mt
sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(
    genes=hl.array(hl.set(hl.flatmap(lambda x: x, [sv_mt.info[field] for field in sv_gene_fields]))),
    restrictive_csq_genes=hl.array(hl.set(hl.flatmap(lambda x: x, [sv_mt.info[field] for field in restrictive_csq_fields]))),
    permissive_csq_genes=hl.array(hl.set(hl.flatmap(lambda x: x, [sv_mt.info[field] for field in permissive_csq_fields]))))
)

# Annotate gene-level annotations for all gene fields (including restrictive and permissive)
sv_gene_agg_mt = get_gene_level_annotations(
    sv_mt, 'genes', sv_gene_fields, inheritance_ht, '', gene_lists=gene_list_dict)

sv_restrictive_gene_agg_mt = get_gene_level_annotations(
    sv_mt, 'restrictive_csq_genes', restrictive_csq_fields, inheritance_ht, 'restrictive_', gene_lists=gene_list_dict)

sv_permissive_gene_agg_mt = get_gene_level_annotations(
    sv_mt, 'permissive_csq_genes', permissive_csq_fields, inheritance_ht, 'permissive_', gene_lists=gene_list_dict)

# Annotate all gene-level results back onto the original sv_mt
sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(
    **{field: sv_gene_agg_mt.rows()[sv_mt.rsid][field]
       for field in ['gene_source', 'inheritance_code', 'gene_list']},
    **{field: sv_restrictive_gene_agg_mt.rows()[sv_mt.rsid][field]
       for field in ['restrictive_gene_source', 'restrictive_inheritance_code', 'restrictive_gene_list'] if field in sv_restrictive_gene_agg_mt.row},
    **{field: sv_permissive_gene_agg_mt.rows()[sv_mt.rsid][field]
       for field in ['permissive_gene_source', 'permissive_inheritance_code', 'permissive_gene_list'] if field in sv_permissive_gene_agg_mt.row}
))

## FLAGS
constrained_gene_list =pd.read_csv(constrained_uri, sep='\t', header=None)[0].tolist()
prec_gene_list =pd.read_csv(prec_uri, sep='\t', header=None)[0].tolist()
hi_gene_list =pd.read_csv(hi_uri, sep='\t', header=None)[0].tolist()
ts_gene_list =pd.read_csv(ts_uri, sep='\t', header=None)[0].tolist()

def get_gene_list_overlap(gene_list, mt):
    '''
    Gets intersection of gene_list and mt.info.genes and returns as a Hail array
    '''
    return hl.array(hl.set(gene_list).intersection(hl.set(mt.info.genes)))

# Annotate constrained, pRec, HI, TS genes
gene_list_map = {'constrained_genes': constrained_gene_list,
    'prec_genes': prec_gene_list,
    'hi_genes': hi_gene_list,
    'ts_genes': ts_gene_list}

sv_mt = sv_mt.annotate_rows(
info=sv_mt.info.annotate(
    **{name: get_gene_list_overlap(gene_list, sv_mt) for name, gene_list in gene_list_map.items()})
)

# Add flags for constrained, pRec, HI, TS as any_{category}
sv_mt = sv_mt.annotate_rows(
info=sv_mt.info.annotate(
    **{f"any_{name.split('_genes')[0]}": sv_mt.info[name].size()>0 for name in gene_list_map.keys()})
)

# Add flag for any genes in gene lists (from gene_list_tsv) and any genes in inheritance_uri
sv_mt = sv_mt.annotate_rows(
    info=sv_mt.info.annotate(
        any_genelist=sv_mt.info.gene_list.filter(hl.is_defined).size()>0,
        any_inheritance=sv_mt.info.inheritance_code.filter(hl.is_defined).size()>0
    )
)

# Frequency flags
# NEW 2/20/2025: Allow for missing gnomAD AFs
sv_mt = sv_mt.annotate_rows(
    info=sv_mt.info.annotate(
        dominant_freq=((hl.max(sv_mt.info.AF)<=dom_af_threshold) & 
            ((sv_mt.info[gnomad_af_field]<=gnomad_af_dom_threshold) | 
            (hl.is_missing(sv_mt.info[gnomad_af_field])))
        ),
        recessive_freq=((hl.max(sv_mt.info.AF)<=rec_af_threshold) & 
            ((sv_mt.info[gnomad_af_field]<=gnomad_af_rec_threshold) |
            (hl.is_missing(sv_mt.info[gnomad_af_field])))
        )
    )
)

# Annotate gnomAD_popmax_AF and gnomad_popmax_freq flag
gnomad_fields = [x for x in list(sv_mt.info) if 'gnomad' in x.lower() 
                and 'ac' not in x.lower() and 'an' not in x.lower() 
                and 'af' in x.lower()]
sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(gnomAD_popmax_AF=hl.max([sv_mt.info[field] for field in gnomad_fields])))
sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(gnomad_popmax_freq=(sv_mt.info.gnomAD_popmax_AF<=gnomad_popmax_af_threshold) |
                                                                        (hl.is_missing(sv_mt.info.gnomAD_popmax_AF))))

# SVLEN flag
suffixes = ['bp', 'kbp', 'mbp', 'gbp', 'tbp', 'pbp']
def humansize(bps):
    i = 0
    while bps >= 1000 and i < len(suffixes)-1:
        bps /= 1000.
        i += 1
    f = ('%.2f' % bps).rstrip('0').rstrip('.')
    return '%s%s' % (f, suffixes[i])

size_threshold_field = f"passes_SVLEN_filter_{humansize(size_threshold)}"
# Flag size threshold
sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(**{size_threshold_field: (sv_mt.info.SVLEN>=size_threshold)}))

# Update header with all new annotations and flags
# Annotations
header['info']['genes'] = {'Description': f"All genes from (union of restrictive_csq_fields and permissive_csq_fields) {', '.join(sv_gene_fields)}.", 'Number': '.', 'Type': 'String'}
header['info']['restrictive_csq_genes'] = {'Description': f"All genes from {', '.join(restrictive_csq_fields)}.", 'Number': '.', 'Type': 'String'}
header['info']['permissive_csq_genes'] = {'Description': f"All genes from {', '.join(permissive_csq_fields)}.", 'Number': '.', 'Type': 'String'}
header['info']['gene_source'] = {'Description': f"Sources for genes in genes field, considered fields: {', '.join(sv_gene_fields)}.", 'Number': '.', 'Type': 'String'}
header['info']['restrictive_gene_source'] = {'Description': f"Sources for genes in genes field, considered fields: {', '.join(restrictive_csq_fields)}.", 'Number': '.', 'Type': 'String'}
header['info']['permissive_gene_source'] = {'Description': f"Sources for genes in genes field, considered fields: {', '.join(permissive_csq_fields)}.", 'Number': '.', 'Type': 'String'}
header['info']['inheritance_code'] = {'Description': f"Inheritance codes from {os.path.basename(inheritance_uri)} for all genes in genes field.", 'Number': '.', 'Type': 'String'}
header['info']['restrictive_inheritance_code'] = {'Description': f"Inheritance codes from {os.path.basename(inheritance_uri)} for genes in restrictive_csq_genes field.", 'Number': '.', 'Type': 'String'}
header['info']['permissive_inheritance_code'] = {'Description': f"Inheritance codes from {os.path.basename(inheritance_uri)} for genes in permissive_csq_genes field.", 'Number': '.', 'Type': 'String'}
header['info']['gene_list'] = {'Description': f"Gene lists for each gene in genes field (&-delimited for multiple gene lists) from {os.path.basename(gene_list_tsv)}.", 'Number': '.', 'Type': 'String'}
header['info']['restrictive_gene_list'] = {'Description': f"Gene lists for each gene in restrictive_csq_genes field (&-delimited for multiple gene lists) from {os.path.basename(gene_list_tsv)}.", 'Number': '.', 'Type': 'String'}
header['info']['permissive_gene_list'] = {'Description': f"Gene lists for each gene in permissive_csq_genes field (&-delimited for multiple gene lists) from {os.path.basename(gene_list_tsv)}.", 'Number': '.', 'Type': 'String'}
header['info']['constrained_genes'] = {'Description': f"All genes in genes field that are in {os.path.basename(constrained_uri)}.", 'Number': '.', 'Type': 'String'}
header['info']['prec_genes'] = {'Description': f"All genes in genes field that are in {os.path.basename(prec_uri)}.", 'Number': '.', 'Type': 'String'}
header['info']['hi_genes'] = {'Description': f"All genes in genes field that are in {os.path.basename(hi_uri)}.", 'Number': '.', 'Type': 'String'}
header['info']['ts_genes'] = {'Description': f"All genes in genes field that are in {os.path.basename(ts_uri)}.", 'Number': '.', 'Type': 'String'}
header['info']['gnomAD_popmax_AF'] = {'Description': f"gnomAD popmax AF taken from fields: {', '.join(gnomad_fields)}.", 'Number': '1', 'Type': 'Float'}
# NEW 5/29/2025: Cohort affected/unaffected annotations
header['info']['n_cohort_het_unaffected'] = {'Description': f"Number of het unaffected individuals (out of {n_tot_samples} total individuals).", 'Number': '1', 'Type': 'Int'}
header['info']['n_cohort_hom_var_unaffected'] = {'Description': f"Number of hom var unaffected individuals (out of {n_tot_samples} total individuals).", 'Number': '1', 'Type': 'Int'}
header['info']['n_cohort_het_affected'] = {'Description': f"Number of het affected individuals (out of {n_tot_samples} total individuals).", 'Number': '1', 'Type': 'Int'}
header['info']['n_cohort_hom_var_affected'] = {'Description': f"Number of hom var affected individuals (out of {n_tot_samples} total individuals).", 'Number': '1', 'Type': 'Int'}
header['info']['cohort_AC_unaffected'] = {'Description': f"AC from unaffected individuals (out of {n_tot_samples} total individuals).", 'Number': '1', 'Type': 'Int'}
header['info']['cohort_AC_affected'] = {'Description': f"AC from affected individuals (out of {n_tot_samples} total individuals).", 'Number': '1', 'Type': 'Int'}
# Flags
header['info']['any_constrained'] = {'Description': f"Any gene in genes field is in {os.path.basename(constrained_uri)}.", 'Number': '0', 'Type': 'Flag'}
header['info']['any_prec'] = {'Description': f"Any gene in genes field is in {os.path.basename(prec_uri)}.", 'Number': '0', 'Type': 'Flag'}
header['info']['any_hi'] = {'Description': f"Any gene in genes field is in {os.path.basename(hi_uri)}.", 'Number': '0', 'Type': 'Flag'}
header['info']['any_ts'] = {'Description': f"Any gene in genes field is in {os.path.basename(ts_uri)}.", 'Number': '0', 'Type': 'Flag'}
header['info']['any_genelist'] = {'Description': f"Any gene in genes field is in gene lists from {os.path.basename(gene_list_tsv)}.", 'Number': '0', 'Type': 'Flag'}
header['info']['any_inheritance'] = {'Description': f"Any gene in genes field is in {os.path.basename(inheritance_uri)}.", 'Number': '0', 'Type': 'Flag'}
header['info']['dominant_freq'] = {'Description': f"Passes cohort AF <= {dom_af_threshold} AND gnomAD AF <= {gnomad_af_dom_threshold}.", 'Number': '0', 'Type': 'Flag'}
header['info']['recessive_freq'] = {'Description': f"Passes cohort AF <= {rec_af_threshold} AND gnomAD AF <= {gnomad_af_rec_threshold}.", 'Number': '0', 'Type': 'Flag'}
header['info']['gnomad_popmax_freq'] = {'Description': f"Passes gnomAD popmax AF <= {gnomad_popmax_af_threshold}.", 'Number': '0', 'Type': 'Flag'}
header['info'][size_threshold_field] = {'Description': f"Passes SVLEN size filter of {humansize(size_threshold).replace('_', ' ')}.", 'Number': '0', 'Type': 'Flag'}

hl.export_vcf(sv_mt, output_filename, metadata=header, tabix=True)