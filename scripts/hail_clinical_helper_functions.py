###
# Created 3/10/2025 for better organization for clinical filtering pipeline(s).

## CHANGE LOG:
'''

'''
###
from typing import Tuple

import pandas as pd
import numpy as np
import hail as hl
import hail.expr.aggregators as agg
from hail.expr import expr_call, expr_float64
from hail.genetics.pedigree import Pedigree
from hail.matrixtable import MatrixTable
from hail.table import Table
from hail.typecheck import numeric, typecheck
from hail.utils.java import Env

def filter_mt(mt, filter_csq=True, filter_impact=True, filter_by_in_gene_list=True):
    '''
    mt: can be trio matrix (tm) or matrix table (mt) but must be transcript-level, not variant-level
    '''
    # filter by Consequence
    if filter_csq:
        exclude_csqs = ['intergenic_variant', 'upstream_gene_variant', 'downstream_gene_variant',
                        'synonymous_variant', 'coding_sequence_variant', 'sequence_variant']
        mt = mt.filter_rows(hl.set(exclude_csqs).intersection(
            hl.set(mt.vep.transcript_consequences.Consequence)).size()!=hl.set(mt.vep.transcript_consequences.Consequence).size())

    # filter only canonical transcript or MANE PLUS CLINICAL
    mt = mt.filter_rows((mt.vep.transcript_consequences.CANONICAL=='YES') | 
                        (mt.vep.transcript_consequences.MANE_PLUS_CLINICAL!=''))
    
    # NEW 3/10/2025: Filter by in gene list
    # NEW 4/18/2025: Make filter_by_in_gene_list optional (default True)
    if filter_by_in_gene_list:
        mt = mt.filter_rows(mt.vep.transcript_consequences.gene_list!='')

    # filter by Impact and splice/noncoding consequence
    if filter_impact:
        splice_vars = ['splice_donor_5th_base_variant', 'splice_region_variant', 'splice_donor_region_variant']
        keep_vars = ['non_coding_transcript_exon_variant']
        mt = mt.filter_rows(
            (hl.set(splice_vars + keep_vars).intersection(
                hl.set(mt.vep.transcript_consequences.Consequence)).size()>0) |
            (hl.array(['HIGH', 'MODERATE']).contains(
            mt.vep.transcript_consequences.IMPACT))
            )
    return mt 

def remove_parent_probands_trio_matrix(tm):
    '''
    Function to bypass peculiarity of Hail's trio_matrix() function when complete_trios=False
    removes "trios" where the "proband" is a parent --> only leaves trios/duos/singletons as entries
    '''
    fathers = tm.father.s.collect()
    mothers = tm.mother.s.collect()
    return tm.filter_cols(hl.array(fathers + mothers).contains(tm.proband.s), keep=False)

def annotate_trio_matrix(phased_tm, mt, pedigree, ped_ht):
    # Annotate trio_status
    complete_trio_probands = [trio.s for trio in pedigree.complete_trios()]
    if len(complete_trio_probands)==0:
        complete_trio_probands = ['']
    phased_tm = phased_tm.annotate_cols(trio_status=hl.if_else(phased_tm.fam_id=='-9', 'not_in_pedigree', 
                                                       hl.if_else(hl.array(complete_trio_probands).contains(phased_tm.id), 'trio', 'non_trio')))
    # Annotate phenotype in MT
    mt = mt.annotate_cols(phenotype=ped_ht[mt.s].phenotype)

    # Get cohort unaffected/affected het and homvar counts
    mt = mt.annotate_rows(**{
        "n_het_unaffected": hl.agg.filter(mt.phenotype=='1', hl.agg.sum(mt.GT.is_het())),
        "n_hom_var_unaffected": hl.agg.filter(mt.phenotype=='1', hl.agg.sum(mt.GT.is_hom_var())),
        "n_het_affected": hl.agg.filter(mt.phenotype=='2', hl.agg.sum(mt.GT.is_het())),
        "n_hom_var_affected": hl.agg.filter(mt.phenotype=='2', hl.agg.sum(mt.GT.is_hom_var()))
    })

    # Get Mendel code/errors, get transmission
    phased_tm = get_mendel_errors(mt, phased_tm, pedigree)
    phased_tm = get_transmission(phased_tm)
    
    # Annotate sex in TM
    phased_tm = phased_tm.annotate_cols(sex=ped_ht[phased_tm.id].sex)
    
    # Annotate affected status/phenotype from pedigree
    phased_tm = phased_tm.annotate_cols(
        proband=phased_tm.proband.annotate(
            phenotype=ped_ht[phased_tm.proband.s].phenotype),
    mother=phased_tm.mother.annotate(
            phenotype=ped_ht[phased_tm.mother.s].phenotype),
    father=phased_tm.father.annotate(
            phenotype=ped_ht[phased_tm.father.s].phenotype))

    affected_cols = ['n_het_unaffected', 'n_hom_var_unaffected', 'n_het_affected', 'n_hom_var_affected']
    phased_tm = phased_tm.annotate_rows(**{col: mt.rows()[phased_tm.row_key][col] 
                                                 for col in affected_cols})

    ## Annotate dominant_gt and recessive_gt
    # denovo
    dom_trio_criteria = ((phased_tm.trio_status=='trio') &  
                         (phased_tm.mendel_code==2))
    # het absent in unaff
    dom_non_trio_criteria = ((phased_tm.trio_status!='trio') & 
                            (phased_tm.n_hom_var_unaffected==0) &
                            (phased_tm.n_het_unaffected==0) & 
                            (phased_tm.proband_entry.GT.is_het())
                            )

    # homozygous and het parents
    rec_trio_criteria = ((phased_tm.trio_status=='trio') &  
                         (phased_tm.proband_entry.GT.is_hom_var()) &
                         (phased_tm.mother_entry.GT.is_het()) &
                         (phased_tm.father_entry.GT.is_het())
                        )  
    # hom and unaff are not hom
    rec_non_trio_criteria = ((phased_tm.trio_status!='trio') &  
                            (phased_tm.n_hom_var_unaffected==0) &
                            (phased_tm.proband_entry.GT.is_hom_var())
                            )

    phased_tm = phased_tm.annotate_entries(dominant_gt=((dom_trio_criteria) | (dom_non_trio_criteria)),
                            recessive_gt=((rec_trio_criteria) | (rec_non_trio_criteria)))
    return phased_tm    

def load_split_vep_consequences(vcf_uri, build):
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
    return mt

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

def get_mendel_errors(mt, phased_tm, pedigree, locus_expr='locus'):
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

def sort_final_merged_output_by_tiers(merged_df):
    # Sort by tier (lower = higher priority)
    def get_len_of_top_numeric_tier(row):
        top_numeric_tier = row.top_numeric_tier
        numeric_tiers = row.numeric_tiers_list
        top_tier = row.tiers_list[numeric_tiers.index(top_numeric_tier)]
        return len(top_tier)

    merged_df['tiers_list'] = merged_df.Tier.str.split(',')  # with * and flags
    merged_df['numeric_tiers_list'] = merged_df.tiers_list.apply(lambda lst: [int(x[0])  if x!='' else 6 for x in lst])  # Assign missing tier to tier 6 for sorting
    merged_df['top_numeric_tier'] = merged_df.numeric_tiers_list.apply(min)
    merged_df['top_tier_len'] = merged_df.apply(get_len_of_top_numeric_tier, axis=1)  # Longer = worse tier!
    merged_df['all_tiers_len'] = merged_df['Tier'].apply(len)  # Longer = worse (more "bad" tiers included)
    # Get best tier for comphets
    merged_df['top_numeric_tier_comphet'] = merged_df['top_numeric_tier']
    merged_df.loc[merged_df.variant_category.str.contains('comphet'), 'top_numeric_tier_comphet'] = merged_df.comphet_ID.map(
                    merged_df[merged_df.variant_category.str.contains('comphet')].groupby('comphet_ID')['top_numeric_tier'].min())

    # Pull flags from Tier column and move to filters column
    def get_flags_from_tier_list(tier_list):
        flags = [';'.join(x.split(';')[1:]) for x in tier_list]
        flags = [x for x in flags if x!='']
        unique_flags = list(set(flags))
        if len(unique_flags)==0:
            return np.nan
        return unique_flags[0]

    def update_filters_with_flags(row):
        if pd.isna(row.tier_filter):  # no new flags to add
            return row.filters
        if pd.isna(row.filters):  # no existing filters
            return row.tier_filter
        return row.filters + ',' + row.tier_filter

    merged_df['tier_filter'] = merged_df.Tier.str.split(',').apply(get_flags_from_tier_list)
    merged_df['filters'] = merged_df.apply(update_filters_with_flags, axis=1)

    # Update Tier column to just have Tier (including *, but NO flags)
    merged_df['Tier'] = merged_df.Tier.str.split(',').apply(lambda lst: [x.split(';')[0] for x in lst]).apply(','.join)

    tmp_tier_cols = ['tiers_list','numeric_tiers_list','top_numeric_tier_comphet','top_numeric_tier','top_tier_len','all_tiers_len','tier_filter']
    # Sort by sample ID first!
    merged_df = merged_df.sort_values(['id','top_numeric_tier_comphet','comphet_ID','top_numeric_tier','top_tier_len','all_tiers_len']).drop(tmp_tier_cols, axis=1)

    # Strip out leading commas
    for col in merged_df.columns:
        if merged_df[col].dtype=='object':
            merged_df[col] = merged_df[col].replace({np.nan: ''}).astype(str).str.lstrip(',').replace({'': np.nan})
    
    return merged_df