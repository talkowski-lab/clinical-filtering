import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", dest="input_uri", help="Input TSV")
parser.add_argument("-o", "--output", dest="output_uri", help="Output filename")
parser.add_argument("--ECNT-threshold", dest="ECNT_threshold", help="Upper bound for ECNT (for Tiers 1 and 2)")
parser.add_argument("--NCount-over-proband-DP-threshold", dest="ncount_over_proband_DP_threshold", help="Upper bound for NCount/proband DP (for Tiers 1 and 2)")
parser.add_argument("--GQ-threshold", dest="GQ_threshold", help="Lower bound for proband GQ (for Tier 1)")
parser.add_argument("-t", "--type", dest="type", help="Type (dominant, recessive, other)")

args = parser.parse_args()

# INPUTS
uri = args.input_uri
output_uri = args.output_uri
ECNT_threshold = int(args.ECNT_threshold)
ncount_over_proband_DP_threshold = float(args.ncount_over_proband_DP_threshold)
GQ_threshold = int(args.GQ_threshold)
inheritance_type = args.type

clnrevstat_one_star_plus = ["criteria_provided,_multiple_submitters,_no_conflicts", "criteria_provided,_single_submitter", "practice_guideline", "reviewed_by_expert_panel"]
PTVs = ['frameshift_variant', 'stop_gained', 'splice_donor_variant', 'splice_acceptor_variant', 'transcript_ablation']

df = pd.read_csv(uri, sep='\t')
# Strip quotes etc. from every column
for col in df.columns:
    if df[col].dtype=='object':
        df[col] = df[col].str.strip('\n').str.replace('\"','').str.replace('[','').str.replace(']','')

# Tier 5: default/lowest tier
df['Tier'] = 5

# Tier 4: Only native NIFS filters (except ABBINOM)
passes_filters_except_abbinom = ((df.filters.isna()) | 
                                 (df.filters.str.contains('ABBINOM')))
df.loc[passes_filters_except_abbinom, 'Tier'] = 4

# Tier 3: Include VUS or Conflicting in ClinVar
vus_or_conflicting_in_clinvar = (df['info.CLNSIG'].str.contains('Uncertain') | df['info.CLNSIG'].str.contains('Conflicting'))
df.loc[passes_filters_except_abbinom & vus_or_conflicting_in_clinvar, 'Tier'] = 3
# CRITERIA FOR BOTH TIER 1 AND TIER 2
is_clinvar_P_LP = ((df['info.CLNSIG'].astype(str).str.contains('athogenic')) 
                & (~df['info.CLNSIG'].astype(str).str.contains('Conflicting')))
is_clnrevstat_one_star_plus = (df['info.CLNREVSTAT'].isin(clnrevstat_one_star_plus))
is_clinvar_P_LP_one_star_plus = is_clinvar_P_LP & is_clnrevstat_one_star_plus
not_in_clinvar = df['info.CLNSIG'].isna()
is_clinvar_P_LP_one_star_plus_or_not_in_clinvar = is_clinvar_P_LP_one_star_plus | not_in_clinvar

in_any_gene_list = ~df['vep.transcript_consequences.gene_list'].isna()

not_in_segdup_str_simplerep = (~df[['info.SEGDUP', 'info.STR', 'info.SIMPLEREP']].any(axis=1))

ncount_over_proband_DP = df['info.NCount'] / df['proband_entry.DP']
passes_ncount_over_proband_DP = (ncount_over_proband_DP < ncount_over_proband_DP_threshold)

passes_ECNT = (df['info.ECNT'] < ECNT_threshold)

high_or_moderate_impact = df['vep.transcript_consequences.IMPACT'].isin(['HIGH','MODERATE'])

not_inframe_indel = (~df['vep.transcript_consequences.Consequence'].str.contains('inframe'))
# GT CRITERIA
if inheritance_type=='recessive':
    # Treat XLR like AD
    xlr_proband_GT = ((df.locus.str.contains('X')) & 
                      (df['proband_entry.GT'].isin(['0/1', '1/1'])))
    # Only HomVar for Tier 1 (except XLR)
    tier_1_proband_GT = (df['proband_entry.GT']=='1/1') | (xlr_proband_GT)
    # Allow Het for Tier 2 if mother has variant
    tier_2_proband_GT = ((tier_1_proband_GT) | 
                         ((df['proband_entry.GT']=='0/1') & (df['mother_entry.GT']!='0/0')))

if inheritance_type=='dominant':
    tier_1_proband_GT = (df['proband_entry.GT'].isin(['0/1', '1/1']))
    tier_2_proband_GT = tier_1_proband_GT

passes_tier_1_and_2 = (is_clinvar_P_LP_one_star_plus_or_not_in_clinvar &
                        in_any_gene_list &
                        not_in_segdup_str_simplerep &
                        passes_ncount_over_proband_DP &
                        passes_ECNT &
                        high_or_moderate_impact &
                        not_inframe_indel &
                        passes_filters_except_abbinom)

# Tier 2: ClinVar P/LP 1*+, in gene list, not in SEGDUP/STR/SIMPLEREP, etc.
df.loc[passes_tier_1_and_2 & tier_2_proband_GT, 'Tier'] = 2
# CRITERIA FOR JUST TIER 1
passes_GQ = (df['proband_entry.GQ']>GQ_threshold)
is_PTV = df['vep.transcript_consequences.Consequence'].apply(lambda csq_str:
                                                             any([csq in csq_str for csq in PTVs]))

# Tier 1: ClinVar P/LP 1*+, in gene list, not in SEGDUP/STR/SIMPLEREP, etc. + GQ filter + PTVs
df.loc[passes_tier_1_and_2 & 
       passes_GQ &
       is_PTV &
       tier_1_proband_GT, 'Tier'] = 1

# Add flags for LQ and VUS
df['Tier'] = df['Tier'].astype(str)
df.loc[~passes_GQ, 'Tier'] = df.loc[~passes_GQ, 'Tier'] + ';LQ'
df.loc[vus_or_conflicting_in_clinvar, 'Tier'] = df.loc[vus_or_conflicting_in_clinvar, 'Tier'] + ';VUS'

df.to_csv(output_uri, sep='\t', index=False)