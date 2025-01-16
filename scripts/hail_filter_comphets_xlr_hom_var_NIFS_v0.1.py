###
# Finally forked from hail_filter_comphets_xlr_hom_var_v0.1.py on 1/13/2025
# to include maternal carrier variants in output and use cluster info 
# (CA in FORMAT field) for comphets. 
# SVs ignored for NIFS for now.

## CHANGE LOG:
'''
1/14/2025:
- added variant_source column
- added workaround for empty CA field (usually cluster 5)
- removed get_transmission function (irrelevant for NIFS)
- moved non-PAR annotations up
- changed export() to to_pandas() to bypass ClassTooLargeException in Hail tables union
1/15/2025:
- changed comphets maternally-inherited clusters from clusters 2-5 to 1-4 (cluster 5 workaround no longer necessary)
- commented out all_csqs and gnomad_popmax_af annotations because now annotated (in INFO) in hail_filter_clinical_variants_NIFS_v0.1.py :)
1/16/2025:
- filter out cluster 5 rows for comphets that were "coming along for the ride"
'''
###

from pyspark.sql import SparkSession
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
omim_uri = sys.argv[6]
sv_gene_fields = sys.argv[7].split(',') 
build = sys.argv[8]
cores = sys.argv[9]  # string
mem = int(np.floor(float(sys.argv[10])))
ad_alt_threshold = int(sys.argv[11])
rec_gene_list_tsv = sys.argv[12]
carrier_gene_list = sys.argv[13]

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

def filter_mt(mt):
    '''
    mt: can be trio matrix (tm) or matrix table (mt) but must be transcript-level, not variant-level
    '''
    # filter by Consequence —— exclude rows where Consequence contains only terms in this list
    exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                    'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']
    mt = mt.filter_rows(hl.set(exclude_csqs).intersection(
        hl.set(mt.vep.transcript_consequences.Consequence)).size()!=hl.set(mt.vep.transcript_consequences.Consequence).size())

    # filter only canonical transcript or MANE PLUS CLINICAL
    mt = mt.filter_rows((mt.vep.transcript_consequences.CANONICAL=='YES') | 
                        (mt.vep.transcript_consequences.MANE_PLUS_CLINICAL!=''))

    # filter by HIGH/MODERATE Impact OR Consequence contains at least one splice/noncoding consequence
    splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
    keep_vars = ['non_coding_transcript_exon_variant']
    mt = mt.filter_rows(
        (hl.set(splice_vars + keep_vars).intersection(
            hl.set(mt.vep.transcript_consequences.Consequence)).size()>0) |
        (hl.array(['HIGH', 'MODERATE']).contains(
        mt.vep.transcript_consequences.IMPACT))
        )
    return mt 

def load_split_vep_consequences(vcf_uri):
    mt = hl.import_vcf(vcf_uri, reference_genome=build, find_replace=('null', ''), force_bgz=True, call_fields=[], array_elements_required=False)
    csq_columns = hl.get_vcf_metadata(vcf_uri)['info']['CSQ']['Description'].split('Format: ')[1].split('|')

    mt = mt.annotate_rows(vep=mt.info)
    transcript_consequences = mt.vep.CSQ.map(lambda x: x.split('\|'))

    transcript_consequences_strs = transcript_consequences.map(lambda x: hl.if_else(hl.len(x)>1, hl.struct(**
                                                        {col: x[i] if col!='Consequence' else x[i].split('&')  
                                                            for i, col in enumerate(csq_columns)}), 
                                                            hl.struct(**{col: hl.missing('str') if col!='Consequence' else hl.array([hl.missing('str')])  
                                                            for i, col in enumerate(csq_columns)})))

    mt = mt.annotate_rows(vep=mt.vep.annotate(transcript_consequences=transcript_consequences_strs))
    mt = mt.annotate_rows(vep=mt.vep.select('transcript_consequences'))
    # NEW 1/9/2025: annotate all_csqs (TODO: many transcripts already filtered out in upstream steps?)
    # NEW 1/15/2025: commented out because now annotated (in INFO) in hail_filter_clinical_variants_NIFS_v0.1.py :)
    # mt = mt.annotate_rows(all_csqs=hl.set(hl.flatmap(lambda x: x, mt.vep.transcript_consequences.Consequence)))
    return mt

## STEP 1: Merge SNV/Indel VCF with SV VCF (or just one of them)
# Load SNV/Indel VCF
if snv_indel_vcf!='NA':
    locus_expr = 'locus'
    snv_mt = load_split_vep_consequences(snv_indel_vcf)

    # Load and merge SNV/Indel ClinVar P/LP VCF
    if clinvar_vcf!='NA':
        clinvar_mt = load_split_vep_consequences(clinvar_vcf) 
        # NEW 1/14/2025: added variant_source —— ClinVar_P/LP or OMIM_recessive or both
        snv_mt_no_clinvar = snv_mt
        snv_mt = snv_mt.union_rows(clinvar_mt).distinct_by_row()
        snv_mt = snv_mt.annotate_rows(variant_source=hl.if_else(hl.is_defined(clinvar_mt.rows()[snv_mt.row_key]),  # if in ClinVar
                                                 hl.if_else(hl.is_defined(snv_mt_no_clinvar.rows()[snv_mt.row_key]), 'ClinVar_P/LP_OMIM_recessive',  # if also in omim_recessive_vcf
                                                            'ClinVar_P/LP'), 'OMIM_recessive'))

    # Explode rows by transcript
    snv_mt = snv_mt.explode_rows(snv_mt.vep.transcript_consequences)

    # NEW 1/9/2025: annotate gnomad_popmax_af after exploding by transcript
    # NEW 1/15/2025: commented out because now annotated (in INFO) in hail_filter_clinical_variants_NIFS_v0.1.py :)
    # gnomad_fields = [x for x in list(snv_mt.vep.transcript_consequences) if 'gnomAD' in x]
    # snv_mt = snv_mt.annotate_rows(
    #     gnomad_popmax_af=hl.max([hl.or_missing(snv_mt.vep.transcript_consequences[gnomad_field]!='',
    #                                     hl.float(snv_mt.vep.transcript_consequences[gnomad_field])) 
    #                             for gnomad_field in gnomad_fields]))
    
    # Filter SNV/Indel MT
    snv_mt = filter_mt(snv_mt)

    # Filter out empty gene fields
    snv_mt = snv_mt.annotate_rows(gene=snv_mt['vep']['transcript_consequences']['SYMBOL'])
    snv_mt = snv_mt.filter_rows(snv_mt.gene!='')

    snv_mt = snv_mt.annotate_rows(variant_type='SNV/Indel', 
                                  gene_source=['vep'])
    
    # NEW 1/7/2025
    # Annotate by gene list(s)
    if rec_gene_list_tsv!='NA':
        gene_list_uris = pd.read_csv(rec_gene_list_tsv, sep='\t', header=None).set_index(0)[1].to_dict()
        gene_lists = {gene_list_name: pd.read_csv(uri, sep='\t', header=None)[0].tolist() 
                    for gene_list_name, uri in gene_list_uris.items()}

        snv_mt = snv_mt.annotate_rows(
            gene_lists=hl.array([hl.or_missing(hl.array(gene_list).contains(snv_mt.gene), gene_list_name) 
                for gene_list_name, gene_list in gene_lists.items()]).filter(hl.is_defined))
    
# Load SV VCF
if sv_vcf!='NA':
    locus_expr = 'locus_interval'
    sv_mt = hl.import_vcf(sv_vcf, reference_genome=build, force_bgz=True, call_fields=[], array_elements_required=False)

    # filter out BNDs
    sv_mt = sv_mt.filter_rows(sv_mt.info.SVTYPE!='BND') 

    # Annotate genes —— ignore genes for CPX SVs
    sv_mt = sv_mt.annotate_rows(gene=hl.or_missing(sv_mt.info.SVTYPE!='CPX', 
                                                   hl.array(hl.set(hl.flatmap(lambda x: x, [sv_mt.info[field] for field in sv_gene_fields])))))
    sv_mt = sv_mt.annotate_rows(variant_type='SV')

    sv_mt = sv_mt.explode_rows(sv_mt.gene)

    # get gene source
    def get_predicted_sources_expr(row_expr, sv_gene_fields):
        return hl.array(
            [hl.or_missing(hl.array(row_expr.info[col]).contains(row_expr.gene), col) for col in sv_gene_fields]
        ).filter(hl.is_defined)

    sv_mt = sv_mt.annotate_rows(gene_source=get_predicted_sources_expr(sv_mt, sv_gene_fields))

    # Make "fake" VEP annotations for SV MT formatting
    if (snv_indel_vcf!='NA'):
        snv_vep_fields = {field: str(snv_mt.vep.transcript_consequences[field].dtype) for field in list(snv_mt.row.vep.transcript_consequences)}
    else:
        snv_vep_fields = {'OMIM_MIM_number': 'array<str>', 'OMIM_inheritance_code': 'str'}
    sv_mt = sv_mt.annotate_rows(vep=hl.struct(transcript_consequences=
            {field: hl.missing(dtype) for field, dtype in snv_vep_fields.items()}))

    # Annotate OMIM in SVs
    omim = hl.import_table(omim_uri).key_by('approvedGeneSymbol')
    sv_mt = sv_mt.key_rows_by('gene')
    sv_mt = sv_mt.annotate_rows(vep=sv_mt.vep.annotate(
        transcript_consequences=sv_mt.vep.transcript_consequences.annotate(
        OMIM_MIM_number=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].mimNumber, ''),
        OMIM_inheritance_code=hl.if_else(hl.is_defined(omim[sv_mt.row_key]), omim[sv_mt.row_key].inheritance_code, ''))))
    sv_mt = sv_mt.key_rows_by('locus', 'alleles')

    # NEW 1/7/2025 
    # Annotate and filter by gene list(s)
    if rec_gene_list_tsv!='NA':
        gene_list_uris = pd.read_csv(rec_gene_list_tsv, sep='\t', header=None).set_index(0)[1].to_dict()
        gene_lists = {gene_list_name: pd.read_csv(uri, sep='\t', header=None)[0].tolist() 
                    for gene_list_name, uri in gene_list_uris.items()}

        sv_mt = sv_mt.annotate_rows(
            gene_lists=hl.array([hl.or_missing(hl.array(gene_list).contains(sv_mt.gene), gene_list_name) 
                for gene_list_name, gene_list in gene_lists.items()]).filter(hl.is_defined))

        # OMIM recessive code
        omim_rec_code = (sv_mt.vep.transcript_consequences.OMIM_inheritance_code.matches('2'))
        # OMIM XLR code
        omim_xlr_code = (sv_mt.vep.transcript_consequences.OMIM_inheritance_code.matches('4'))
        in_gene_list = (sv_mt.gene_lists.size()>0)
        sv_mt = sv_mt.filter_rows(omim_rec_code | omim_xlr_code | in_gene_list)

# Unify SNV/Indel MT and SV MT INFO (row) and entry fields
if (snv_indel_vcf!='NA') and (sv_vcf!='NA'):
    sv_info_fields, sv_entry_fields = list(sv_mt.row.info), list(sv_mt.entry)
    snv_info_fields, snv_entry_fields = list(snv_mt.row.info), list(snv_mt.entry)

    sv_missing_entry_fields = {field: str(snv_mt[field].dtype) for field in np.setdiff1d(snv_entry_fields, sv_entry_fields)}
    snv_missing_entry_fields = {field: str(sv_mt[field].dtype) for field in np.setdiff1d(sv_entry_fields, snv_entry_fields)}

    sv_missing_info_fields = {field: str(snv_mt.info[field].dtype) for field in np.setdiff1d(snv_info_fields, sv_info_fields)}
    snv_missing_info_fields = {field: str(sv_mt.info[field].dtype) for field in np.setdiff1d(sv_info_fields, snv_info_fields)}

    sv_mt = sv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in sv_missing_entry_fields.items()})
    snv_mt = snv_mt.annotate_entries(**{field: hl.missing(dtype) for field, dtype in snv_missing_entry_fields.items()})

    sv_mt = sv_mt.select_entries(*sorted(list(sv_mt.entry)))
    snv_mt = snv_mt.select_entries(*sorted(list(snv_mt.entry)))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in sv_missing_info_fields.items()}))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.annotate(**{field: hl.missing(dtype) for field, dtype in snv_missing_info_fields.items()}))

    sv_mt = sv_mt.annotate_rows(info=sv_mt.info.select(*sorted(list(sv_mt.info))))
    snv_mt = snv_mt.annotate_rows(info=snv_mt.info.select(*sorted(list(snv_mt.info))))

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
if sv_vcf!='NA':
    # Change locus to locus_interval to include END for SVs
    merged_mt = merged_mt.annotate_rows(end=hl.if_else(hl.is_defined(merged_mt.info.END2), merged_mt.info.END2, merged_mt.info.END))
    # Account for INS having same END but different SVLEN
    merged_mt = merged_mt.annotate_rows(end=hl.if_else(merged_mt.info.SVTYPE=='INS', merged_mt.end + merged_mt.info.SVLEN, merged_mt.end))
    # Account for empty END for SNV/Indels
    merged_mt = merged_mt.annotate_rows(end=hl.if_else(merged_mt.variant_type=='SNV/Indel', merged_mt.locus.position, merged_mt.end))
    merged_mt = merged_mt.key_rows_by()
    merged_mt = merged_mt.annotate_rows(locus_interval=hl.locus_interval(contig=merged_mt.locus.contig, 
                                                                              start=merged_mt.locus.position,
                                                                              end=merged_mt.end, reference_genome=build))
    merged_mt = merged_mt.key_rows_by(locus_expr, 'alleles')

# NEW 1/14/2025: Annotate PAR status (moved up from end of script)
merged_mt = merged_mt.annotate_rows(in_non_par=~(merged_mt.locus.in_autosome_or_par()))

## EDITED HAIL FUNCTIONS
# EDITED: don't check locus struct
@typecheck(dataset=MatrixTable, method=str, tolerate_generic_locus=bool)
def require_biallelic(dataset, method, tolerate_generic_locus: bool = False) -> MatrixTable:
    return dataset._select_rows(
        method,
        hl.case()
        .when(dataset.alleles.length() == 2, dataset._rvrow)
        .or_error(
            f"'{method}' expects biallelic variants ('alleles' field of length 2), found "
            + hl.str(dataset.locus)
            + ", "
            + hl.str(dataset.alleles)
        ),
    )

# EDITED: custom require_biallelic function
@typecheck(call=expr_call, pedigree=Pedigree)
def mendel_errors(call, pedigree) -> Tuple[Table, Table, Table, Table]:
    r"""Find Mendel errors; count per variant, individual and nuclear family.

    .. include:: ../_templates/req_tstring.rst

    .. include:: ../_templates/req_tvariant.rst

    .. include:: ../_templates/req_biallelic.rst

    Examples
    --------

    Find all violations of Mendelian inheritance in each (dad, mom, kid) trio in
    a pedigree and return four tables (all errors, errors by family, errors by
    individual, errors by variant):

    >>> ped = hl.Pedigree.read('data/trios.fam')
    >>> all_errors, per_fam, per_sample, per_variant = hl.mendel_errors(dataset['GT'], ped)

    Export all mendel errors to a text file:

    >>> all_errors.export('output/all_mendel_errors.tsv')

    Annotate columns with the number of Mendel errors:

    >>> annotated_samples = dataset.annotate_cols(mendel=per_sample[dataset.s])

    Annotate rows with the number of Mendel errors:

    >>> annotated_variants = dataset.annotate_rows(mendel=per_variant[dataset.locus, dataset.alleles])

    Notes
    -----

    The example above returns four tables, which contain Mendelian violations
    grouped in various ways. These tables are modeled after the `PLINK mendel
    formats <https://www.cog-genomics.org/plink2/formats#mendel>`_, resembling
    the ``.mendel``, ``.fmendel``, ``.imendel``, and ``.lmendel`` formats,
    respectively.

    **First table:** all Mendel errors. This table contains one row per Mendel
    error, keyed by the variant and proband id.

        - `locus` (:class:`.tlocus`) -- Variant locus, key field.
        - `alleles` (:class:`.tarray` of :py:data:`.tstr`) -- Variant alleles, key field.
        - (column key of `dataset`) (:py:data:`.tstr`) -- Proband ID, key field.
        - `fam_id` (:py:data:`.tstr`) -- Family ID.
        - `mendel_code` (:py:data:`.tint32`) -- Mendel error code, see below.

    **Second table:** errors per nuclear family. This table contains one row
    per nuclear family, keyed by the parents.

        - `pat_id` (:py:data:`.tstr`) -- Paternal ID. (key field)
        - `mat_id` (:py:data:`.tstr`) -- Maternal ID. (key field)
        - `fam_id` (:py:data:`.tstr`) -- Family ID.
        - `children` (:py:data:`.tint32`) -- Number of children in this nuclear family.
        - `errors` (:py:data:`.tint64`) -- Number of Mendel errors in this nuclear family.
        - `snp_errors` (:py:data:`.tint64`) -- Number of Mendel errors at SNPs in this
          nuclear family.

    **Third table:** errors per individual. This table contains one row per
    individual. Each error is counted toward the proband, father, and mother
    according to the `Implicated` in the table below.

        - (column key of `dataset`) (:py:data:`.tstr`) -- Sample ID (key field).
        - `fam_id` (:py:data:`.tstr`) -- Family ID.
        - `errors` (:py:data:`.tint64`) -- Number of Mendel errors involving this
          individual.
        - `snp_errors` (:py:data:`.tint64`) -- Number of Mendel errors involving this
          individual at SNPs.

    **Fourth table:** errors per variant.

        - `locus` (:class:`.tlocus`) -- Variant locus, key field.
        - `alleles` (:class:`.tarray` of :py:data:`.tstr`) -- Variant alleles, key field.
        - `errors` (:py:data:`.tint64`) -- Number of Mendel errors in this variant.

    This method only considers complete trios (two parents and proband with
    defined sex). The code of each Mendel error is determined by the table
    below, extending the
    `Plink classification <https://www.cog-genomics.org/plink2/basic_stats#mendel>`__.

    In the table, the copy state of a locus with respect to a trio is defined
    as follows, where PAR is the `pseudoautosomal region
    <https://en.wikipedia.org/wiki/Pseudoautosomal_region>`__ (PAR) of X and Y
    defined by the reference genome and the autosome is defined by
    :meth:`~.LocusExpression.in_autosome`.

    - Auto -- in autosome or in PAR or female child
    - HemiX -- in non-PAR of X and male child
    - HemiY -- in non-PAR of Y and male child

    `Any` refers to the set \{ HomRef, Het, HomVar, NoCall \} and `~`
    denotes complement in this set.

    +------+---------+---------+--------+----------------------------+
    | Code | Dad     | Mom     | Kid    | Copy State | Implicated    |
    +======+=========+=========+========+============+===============+
    |    1 | HomVar  | HomVar  | Het    | Auto       | Dad, Mom, Kid |
    +------+---------+---------+--------+------------+---------------+
    |    2 | HomRef  | HomRef  | Het    | Auto       | Dad, Mom, Kid |
    +------+---------+---------+--------+------------+---------------+
    |    3 | HomRef  | ~HomRef | HomVar | Auto       | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    4 | ~HomRef | HomRef  | HomVar | Auto       | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    5 | HomRef  | HomRef  | HomVar | Auto       | Kid           |
    +------+---------+---------+--------+------------+---------------+
    |    6 | HomVar  | ~HomVar | HomRef | Auto       | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    7 | ~HomVar | HomVar  | HomRef | Auto       | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |    8 | HomVar  | HomVar  | HomRef | Auto       | Kid           |
    +------+---------+---------+--------+------------+---------------+
    |    9 | Any     | HomVar  | HomRef | HemiX      | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |   10 | Any     | HomRef  | HomVar | HemiX      | Mom, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |   11 | HomVar  | Any     | HomRef | HemiY      | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+
    |   12 | HomRef  | Any     | HomVar | HemiY      | Dad, Kid      |
    +------+---------+---------+--------+------------+---------------+

    See Also
    --------
    :func:`.mendel_error_code`

    Parameters
    ----------
    dataset : :class:`.MatrixTable`
    pedigree : :class:`.Pedigree`

    Returns
    -------
    (:class:`.Table`, :class:`.Table`, :class:`.Table`, :class:`.Table`)
    """
    source = call._indices.source
    if not isinstance(source, MatrixTable):
        raise ValueError(
            "'mendel_errors': expected 'call' to be an expression of 'MatrixTable', found {}".format(
                "expression of '{}'".format(source.__class__) if source is not None else 'scalar expression'
            )
        )

    source = source.select_entries(__GT=call)
    dataset = require_biallelic(source, 'mendel_errors', tolerate_generic_locus=True)
    tm = hl.trio_matrix(dataset, pedigree, complete_trios=True)
    tm = tm.select_entries(
        mendel_code=hl.mendel_error_code(
            tm.locus, tm.is_female, tm.father_entry['__GT'], tm.mother_entry['__GT'], tm.proband_entry['__GT']
        )
    )
    ck_name = next(iter(source.col_key))
    tm = tm.filter_entries(hl.is_defined(tm.mendel_code))
    tm = tm.rename({'id': ck_name})

    entries = tm.entries()

    table1 = entries.select('fam_id', 'mendel_code')

    t2 = tm.annotate_cols(errors=hl.agg.count(), snp_errors=hl.agg.count_where(hl.is_snp(tm.alleles[0], tm.alleles[1])))
    table2 = t2.key_cols_by().cols()
    table2 = table2.select(
        pat_id=table2.father[ck_name],
        mat_id=table2.mother[ck_name],
        fam_id=table2.fam_id,
        errors=table2.errors,
        snp_errors=table2.snp_errors,
    )
    table2 = table2.group_by('pat_id', 'mat_id').aggregate(
        fam_id=hl.agg.take(table2.fam_id, 1)[0],
        children=hl.int32(hl.agg.count()),
        errors=hl.agg.sum(table2.errors),
        snp_errors=hl.agg.sum(table2.snp_errors),
    )
    table2 = table2.annotate(
        errors=hl.or_else(table2.errors, hl.int64(0)), snp_errors=hl.or_else(table2.snp_errors, hl.int64(0))
    )

    # in implicated, idx 0 is dad, idx 1 is mom, idx 2 is child
    implicated = hl.literal(
        [
            [0, 0, 0],  # dummy
            [1, 1, 1],
            [1, 1, 1],
            [1, 0, 1],
            [0, 1, 1],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [0, 0, 1],
            [0, 1, 1],
            [0, 1, 1],
            [1, 0, 1],
            [1, 0, 1],
        ],
        dtype=hl.tarray(hl.tarray(hl.tint64)),
    )

    table3 = (
        tm.annotate_cols(
            all_errors=hl.or_else(hl.agg.array_sum(implicated[tm.mendel_code]), [0, 0, 0]),
            snp_errors=hl.or_else(
                hl.agg.filter(hl.is_snp(tm.alleles[0], tm.alleles[1]), hl.agg.array_sum(implicated[tm.mendel_code])),
                [0, 0, 0],
            ),
        )
        .key_cols_by()
        .cols()
    )

    table3 = table3.select(
        xs=[
            hl.struct(**{
                ck_name: table3.father[ck_name],
                'fam_id': table3.fam_id,
                'errors': table3.all_errors[0],
                'snp_errors': table3.snp_errors[0],
            }),
            hl.struct(**{
                ck_name: table3.mother[ck_name],
                'fam_id': table3.fam_id,
                'errors': table3.all_errors[1],
                'snp_errors': table3.snp_errors[1],
            }),
            hl.struct(**{
                ck_name: table3.proband[ck_name],
                'fam_id': table3.fam_id,
                'errors': table3.all_errors[2],
                'snp_errors': table3.snp_errors[2],
            }),
        ]
    )
    table3 = table3.explode('xs')
    table3 = table3.select(**table3.xs)
    table3 = (
        table3.group_by(ck_name, 'fam_id')
        .aggregate(errors=hl.agg.sum(table3.errors), snp_errors=hl.agg.sum(table3.snp_errors))
        .key_by(ck_name)
    )

    table4 = tm.select_rows(errors=hl.agg.count_where(hl.is_defined(tm.mendel_code))).rows()

    return table1, table2, table3, table4

## STEP 2: Get CompHets
# Mendel errors
def get_mendel_errors(mt, phased_tm, pedigree):
    all_errors, per_fam, per_sample, per_variant = mendel_errors(mt['GT'], pedigree)  # edited Hail function, see above
    all_errors_mt = all_errors.key_by().to_matrix_table(row_key=[locus_expr,'alleles'], col_key=['s'], col_fields=['fam_id'])
    phased_tm = phased_tm.annotate_entries(mendel_code=all_errors_mt[phased_tm.row_key, phased_tm.col_key].mendel_code)
    return phased_tm

def phase_by_transmission_aggregate_by_gene(tm, mt, pedigree):
    # filter out calls that are hom ref in proband
    tm = tm.filter_entries(tm.proband_entry.GT.is_non_ref())

    phased_tm = hl.experimental.phase_trio_matrix_by_transmission(tm, call_field='GT', phased_call_field='PBT_GT')
    phased_tm = get_mendel_errors(mt, phased_tm, pedigree)
    phased_tm = phased_tm.key_rows_by(locus_expr,'alleles','gene')

    # NEW 1/14/2025: filter_entries before aggregation
    phased_tm = phased_tm.filter_entries(hl.is_defined(phased_tm.proband_entry.GT))
    gene_agg_phased_tm = (phased_tm.group_rows_by(phased_tm.gene)
        .aggregate_rows(locus_alleles = hl.agg.collect(phased_tm.row_key),
                       variant_type = hl.agg.collect_as_set(phased_tm.variant_type),  # SET
                       variant_source = hl.agg.collect_as_set(phased_tm.variant_source)  # SET
                       )  # NEW 1/14/2025: added variant_source
        .aggregate_entries(all_locus_alleles=hl.agg.collect(phased_tm.row_key),
                          proband_PBT_GT = hl.agg.collect(phased_tm.proband_entry.PBT_GT),
                          proband_GT = hl.agg.collect(phased_tm.proband_entry.GT),
                          proband_CA = hl.agg.collect(phased_tm.proband_entry.CA))  # NEW FOR NIFS CLUSTER
        ).result()
    
    return phased_tm, gene_agg_phased_tm

def get_subset_tm(mt, samples, pedigree, keep=True, complete_trios=False):
    subset_mt = mt.filter_cols(hl.array(samples).contains(mt.s), keep=keep)

    # remove variants missing in subset samples
    subset_mt = hl.variant_qc(subset_mt)
    subset_mt = subset_mt.filter_rows(subset_mt.variant_qc.AC[1]>0)
    subset_mt = subset_mt.drop('variant_qc')

    subset_tm = hl.trio_matrix(subset_mt, pedigree, complete_trios=complete_trios)
    # filter by AD of alternate allele in proband
    if 'AD' in list(subset_mt.entry):
        subset_tm = subset_tm.filter_entries(subset_tm.proband_entry.AD[1]>=ad_alt_threshold)
    return subset_mt, subset_tm

def get_non_trio_comphets(mt):  # EDITED FOR NIFS
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

    # NIFS-specific: comphets must have at least one variant in cluster 0 (fetal 0/1, maternal 0/0) 
    # and one maternally-inherited variant (clusters 1-4)
    in_cluster0 = (hl.set([0]).intersection(hl.set(potential_comp_hets_non_trios.proband_CA)).size()>0)
    in_maternally_inherited_cluster = (hl.set([1, 2, 3, 4]).intersection(hl.set(potential_comp_hets_non_trios.proband_CA)).size()>0)  
    potential_comp_hets_non_trios = potential_comp_hets_non_trios.filter_entries(
        in_cluster0 & in_maternally_inherited_cluster  # NEW FOR NIFS     
    )

    # Annotate non-gene-aggregated TM using potential comphets from gene-aggregated TM 
    non_trio_phased_tm = non_trio_phased_tm.key_rows_by(locus_expr, 'alleles', 'gene')
    non_trio_phased_tm = non_trio_phased_tm.annotate_entries(locus_alleles=  # EDITED
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].all_locus_alleles,
                                                            proband_GT=
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].proband_GT,
                                                            proband_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].proband_GT),
                                                            proband_PBT_GT_set=hl.set(
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].proband_PBT_GT),
                                                            proband_CA=  # NEW FOR NIFS
        potential_comp_hets_non_trios[non_trio_phased_tm.row_key, non_trio_phased_tm.col_key].proband_CA)
    
    # NEW 1/14/2025: Annotate variant_type and variant_source as comma-separated strings of unique values (basically per gene)
    non_trio_phased_tm = non_trio_phased_tm.annotate_rows(variant_type=  
        hl.str(', ').join(hl.sorted(hl.array(potential_comp_hets_non_trios.rows()[non_trio_phased_tm.row_key].variant_type))),
                                                        variant_source= 
        hl.str(', ').join(hl.sorted(hl.array(potential_comp_hets_non_trios.rows()[non_trio_phased_tm.row_key].variant_source)))
    )  

    # Apply same NIFS-specific comphet filter on entries to non-gene-aggregated TM
    in_cluster0 = (hl.set([0]).intersection(hl.set(non_trio_phased_tm.proband_CA)).size()>0)
    in_maternally_inherited_cluster = (hl.set([1, 2, 3, 4]).intersection(hl.set(non_trio_phased_tm.proband_CA)).size()>0) 
    non_trio_phased_tm = non_trio_phased_tm.filter_entries((hl.set(non_trio_phased_tm.locus_alleles).size()>1) &
                                                    in_cluster0 & in_maternally_inherited_cluster  # NEW FOR NIFS     
                                                    )  
    # NEW 1/16/2025: filter out cluster 5 rows for comphets that were "coming along for the ride"
    non_trio_phased_tm = non_trio_phased_tm.filter_rows(hl.array([0, 1, 2, 3, 4]).contains(non_trio_phased_tm.info.CA))

    # Grab rows (variants) from non-gene-aggregated TM, of potential comphets from gene-aggregated TM
    phased_tm_comp_hets_non_trios = non_trio_phased_tm.semi_join_rows(potential_comp_hets_non_trios.rows()).key_rows_by(locus_expr, 'alleles')
    return phased_tm_comp_hets_non_trios.drop('locus_alleles')

def get_trio_comphets(mt):  # IRRELEVANT FOR NIFS —— IGNORE
    trio_mt, trio_tm = get_subset_tm(mt, trio_samples, trio_pedigree, keep=True, complete_trios=True)
    trio_phased_tm, trio_gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(trio_tm, trio_mt, trio_pedigree)

    # different criteria for trios (requires phasing)
    potential_comp_hets_trios = trio_gene_agg_phased_tm.filter_rows(
        hl.agg.count_where(hl.set(trio_gene_agg_phased_tm.proband_PBT_GT).size()>1)>0
    )
    potential_comp_hets_trios = potential_comp_hets_trios.explode_rows(potential_comp_hets_trios.locus_alleles)
    potential_comp_hets_trios = potential_comp_hets_trios.key_rows_by(potential_comp_hets_trios.locus_alleles[locus_expr], potential_comp_hets_trios.locus_alleles.alleles, 'gene')

    potential_comp_hets_trios = potential_comp_hets_trios.filter_entries(hl.set(potential_comp_hets_trios.proband_PBT_GT).size()>1)

    trio_phased_tm = trio_phased_tm.key_rows_by(locus_expr, 'alleles', 'gene')
    trio_phased_tm = trio_phased_tm.annotate_entries(proband_GT=
        potential_comp_hets_trios[trio_phased_tm.row_key, trio_phased_tm.col_key].proband_GT,
                                                               proband_GT_set=hl.set(
        potential_comp_hets_trios[trio_phased_tm.row_key, trio_phased_tm.col_key].proband_GT),
                                                               proband_PBT_GT_set=hl.set(
        potential_comp_hets_trios[trio_phased_tm.row_key, trio_phased_tm.col_key].proband_PBT_GT))

    trio_phased_tm = trio_phased_tm.filter_entries(trio_phased_tm.proband_PBT_GT_set.size()>1)  
    phased_tm_comp_hets_trios = trio_phased_tm.semi_join_rows(potential_comp_hets_trios.rows()).key_rows_by(locus_expr, 'alleles')
    return phased_tm_comp_hets_trios

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
    merged_trio_comphets = merged_trio_comphets.annotate_cols(trio_status='trio')
    merged_comphets = merged_trio_comphets.entries()

if len(non_trio_samples)>0:
    merged_non_trio_comphets = get_non_trio_comphets(comphet_mt)
    merged_non_trio_comphets = merged_non_trio_comphets.annotate_cols(trio_status=
                                                              hl.if_else(merged_non_trio_comphets.fam_id=='-9', 
                                                              'not_in_pedigree', 'non_trio'))
    merged_comphets = merged_non_trio_comphets.entries()

if (len(trio_samples)>0) and (len(non_trio_samples)>0):
    merged_comphets = merged_trio_comphets.entries().union(merged_non_trio_comphets.entries())

if len(trio_samples)==0:
    trio_samples = ['']

# Trio matrix
merged_tm = hl.trio_matrix(merged_mt, pedigree, complete_trios=False)

# filter by AD of alternate allele in proband
# TODO: might filter out all if AD is empty (e.g. all SVs after merging)
if 'AD' in list(merged_mt.entry):
    merged_tm = merged_tm.filter_entries(merged_tm.proband_entry.AD[1]>=ad_alt_threshold)

gene_phased_tm, gene_agg_phased_tm = phase_by_transmission_aggregate_by_gene(merged_tm, merged_mt, pedigree)
gene_phased_tm = gene_phased_tm.annotate_cols(trio_status=hl.if_else(gene_phased_tm.fam_id=='-9', 'not_in_pedigree', 
                                                   hl.if_else(hl.array(trio_samples).contains(gene_phased_tm.id), 'trio', 'non_trio')))

# NEW 1/13/2025: maternal carrier variants
# in carrier gene list and in cluster where mother is het (clusters 1-3)
carrier_genes = pd.read_csv(carrier_gene_list, sep='\t', header=None)[0].tolist()
mat_carrier = gene_phased_tm.filter_entries((hl.array(carrier_genes).contains(gene_phased_tm.vep.transcript_consequences.SYMBOL)) &
                           (hl.array([1, 2, 3]).contains(gene_phased_tm.proband_entry.CA))).key_rows_by(locus_expr, 'alleles').entries()

# XLR only
xlr_phased_tm = gene_phased_tm.filter_rows(gene_phased_tm.vep.transcript_consequences.OMIM_inheritance_code.matches('4'))   # OMIM XLR
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

# NEW 1/14/2025: use to_pandas() to bypass ClassTooLargeException in Hail tables union
merged_comphets_df = merged_comphets.drop('proband_GT','proband_GT_set','proband_PBT_GT_set','proband_CA').flatten().to_pandas()
xlr_phased_df = xlr_phased.flatten().to_pandas()
phased_hom_var_df = phased_hom_var.flatten().to_pandas()
mat_carrier_df = mat_carrier.flatten().to_pandas()

merged_comphets_xlr_hom_var_mat_carrier_df = pd.concat([merged_comphets_df, xlr_phased_df, phased_hom_var_df, mat_carrier_df])

output_filename = f"{prefix}_{variant_types}_comp_hets_xlr_hom_var_mat_carrier.tsv.gz"
if len(output_filename)>(os.pathconf('/', 'PC_NAME_MAX')-len('/cromwell_root/.')):  # if filename too long
    output_filename = f"{variant_types}_comp_hets_xlr_hom_var_mat_carrier.tsv.gz"

merged_comphets_xlr_hom_var_mat_carrier_df.to_csv(output_filename, sep='\t', index=False)