###
# Created 4/19/2025 to adjust NIFS tiering to small variant pipeline.

## CHANGE LOG:
'''
4/19/2025:
- change is_female field to sex to account for missing sex information
'''
###

import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", dest="input_uri", help="Input TSV")
parser.add_argument("-p", "--prefix", dest="prefix", help="Prefix for output filename")
parser.add_argument("-t", "--type", dest="type", help="Type (dominant, recessive, other)")

args = parser.parse_args()

# INPUTS
uri = args.input_uri
prefix = args.prefix
inheritance_type = args.type

clnrevstat_one_star_plus = ['practice_guideline', 'reviewed_by_expert_panel', 'criteria_provided,_multiple_submitters,_no_conflicts', 'criteria_provided,_conflicting_classifications', 'criteria_provided,_single_submitter']

df = pd.read_csv(uri, sep='\t')
# Strip quotes etc. from every column
for col in df.columns:
    if df[col].dtype=='object':
        df[col] = df[col].replace({np.nan: ''}).astype(str).str.strip('\n').str.replace('\"','').str.replace('[','').str.replace(']','')

# Tier 7: default/lowest tier
df['Tier'] = 7

# Tier 6: Only native NIFS filters
passes_filters = (df.filters=='')
df.loc[passes_filters, 'Tier'] = 6

# ClinVar criteria
is_clinvar_P_LP = ((df['info.CLNSIG'].astype(str).str.contains('athogenic')) 
                & (~df['info.CLNSIG'].astype(str).str.contains('Conflicting')))
is_clnrevstat_one_star_plus = (df['info.CLNREVSTAT'].isin(clnrevstat_one_star_plus))
is_clinvar_P_LP_one_star_plus = is_clinvar_P_LP & is_clnrevstat_one_star_plus
is_not_clinvar_B_LB = (~df['info.CLNSIG'].astype(str).str.contains('enign'))
# NEW 7/11/2025: Require ClinVar gene match for Tiers 1-3
# NEW 7/28/2025: Require exact match for ClinVar gene match
if df.empty:
    clinvar_gene_matches = False
else:
    clinvar_gene_matches = (
        df.apply(
            lambda row: row['vep.transcript_consequences.SYMBOL'] in [
                gene.split(':')[0] for gene in row['info.GENEINFO'].split('|')
            ],
            axis=1
        )
        | (df['info.GENEINFO'] == '')
    )

# CRITERIA FOR TIERS 3-5
vus_or_conflicting_in_clinvar = (df['info.CLNSIG'].str.contains('Uncertain') | df['info.CLNSIG'].str.contains('Conflicting'))

# Tier 5: Include VUS or Conflicting in ClinVar
df.loc[passes_filters &
        (vus_or_conflicting_in_clinvar | is_clinvar_P_LP_one_star_plus), 'Tier'] = 5

# CRITERIA FOR TIERS 1-4: STRONG/DEFINITIVE
has_strong_definitive_evidence = (df['vep.transcript_consequences.genCC_classification'].astype(str).str.match('Strong|Definitive'))

# Tier 4: Include VUS or Conflicting in ClinVar AND Strong/Definitive
df.loc[passes_filters & has_strong_definitive_evidence &
       (vus_or_conflicting_in_clinvar | is_clinvar_P_LP_one_star_plus), 'Tier'] = 4
# NEW 4/11/2025: Tier 3: Only include Conflicting with at least one P/LP in CLNSIGCONF
conflicting_P_LP = ((df['info.CLNSIGCONF'].astype(str).str.contains('athogenic')) |  # if Conflicting, must have at least one P/LP
                    (df['info.CLNSIGCONF'].isna()))  # if not Conflicting, CLNSIGCONF is empty
df.loc[passes_filters & has_strong_definitive_evidence & 
        (((vus_or_conflicting_in_clinvar & conflicting_P_LP) | is_clinvar_P_LP_one_star_plus) & 
                clinvar_gene_matches), 'Tier'] = 3

# CRITERIA FOR BOTH TIER 1 AND TIER 2
high_impact = (df['vep.transcript_consequences.IMPACT']=='HIGH')  # Tier 1
high_or_moderate_impact = df['vep.transcript_consequences.IMPACT'].isin(['HIGH','MODERATE'])  # Tier 2

# GT CRITERIA
if inheritance_type=='recessive':
    ar_proband_GT = df['proband_entry.GT'].isin(['1/1','1|1'])
    # Males: treat XLR like AD (proband has alt allele) except for PARs (treat like AR)
    xlr_proband_cond_male_par = ((df.sex=='1') & 
                                (df.locus.str.contains('X')) &
                                (~df.in_non_par) & 
                                (ar_proband_GT))
    xlr_proband_cond_male_non_par = ((df.sex=='1') & 
                                (df.locus.str.contains('X')) &
                                (df.in_non_par) & 
                                (df['proband_entry.GT'].str.contains('1')))
    # Females: treat XLR like AR
    xlr_proband_cond_female = ((df.sex=='2') & 
                               (df.locus.str.contains('X')) &
                               (ar_proband_GT))    
    xlr_proband_GT = ((xlr_proband_cond_male_par | xlr_proband_cond_male_non_par) |
                      (xlr_proband_cond_female))
    # Tier 1: Only HomVar (except XLR)
    tier_1_proband_GT = (ar_proband_GT) | (xlr_proband_GT)
    # Tier 2: Same as Tier 1 but allow 0/1 for comphets
    comphet_proband_het = (df['variant_category']=='comphet') & (df['proband_entry.GT'].isin(['0/1','0|1']))
    tier_2_proband_GT = (tier_1_proband_GT) | (comphet_proband_het) 

elif inheritance_type=='dominant':
    # Tier 1: Proband has alt allele
    tier_1_proband_GT = (df['proband_entry.GT'].str.contains('1'))
    # Tier 2: Same as Tier 1
    tier_2_proband_GT = tier_1_proband_GT

else:  # 'other' inheritance_type --> force into Tiers 3 and above
    tier_1_proband_GT = (df['proband_entry.GT'].isna())  # Assumes proband_entry.GT is defined for all
    tier_2_proband_GT = tier_1_proband_GT
    
passes_tier_1_and_2 = (is_not_clinvar_B_LB &
                        ~vus_or_conflicting_in_clinvar &
                        passes_filters &
                        has_strong_definitive_evidence)

# Tier 2: ClinVar P/LP 1*+ OR HIGH/MODERATE IMPACT, not in SEGDUP, etc.
df.loc[((is_clinvar_P_LP_one_star_plus & clinvar_gene_matches) | high_or_moderate_impact) &
        passes_tier_1_and_2 & 
        tier_2_proband_GT, 'Tier'] = 2

# Tier 1: ClinVar P/LP 1*+ OR HIGH IMPACT, not in SEGDUP/STR/SIMPLEREP, etc.
df.loc[((is_clinvar_P_LP_one_star_plus & clinvar_gene_matches) | high_impact) &
        passes_tier_1_and_2 & 
        tier_1_proband_GT, 'Tier'] = 1

# Add flags for LQ and VUS
df['Tier'] = df['Tier'].astype(str)
# df.loc[(~passes_ECNT | ~passes_ncount_over_proband_DP), 'Tier'] = df.loc[(~passes_ECNT | ~passes_ncount_over_proband_DP), 'Tier'] + ';LQ'
df.loc[vus_or_conflicting_in_clinvar, 'Tier'] = df.loc[vus_or_conflicting_in_clinvar, 'Tier'] + ';VUS'

# Add maternal_carrier column: same as Tier 2 except tier_2_proband_GT (alt allele) filter
is_maternal_variant = df['mother_entry.GT'].str.contains('1')
df['maternal_carrier'] = ((is_clinvar_P_LP_one_star_plus | high_or_moderate_impact) &
        passes_tier_1_and_2 & 
        is_maternal_variant)

output_filename = f"{prefix}_tiers.tsv.gz"
df.to_csv(output_filename, sep='\t', index=False)
