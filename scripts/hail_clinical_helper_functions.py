###
# Created 3/10/2025 for better organization for clinical filtering pipeline(s).

## CHANGE LOG:
'''

'''
###
import hail as hl

def filter_mt(mt, filter_csq=True, filter_impact=True):
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