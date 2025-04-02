import datetime
import pandas as pd
import numpy as np
import sys
import ast
import os
import argparse

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='input_uris', help='Comma-separated list of all input TSVs')
parser.add_argument('-p', dest='prefix', help='Prefix for output filename')
parser.add_argument('-g', dest='gene_phenotype_map', help='TSV with gene_symbol, disease_title_recessive, disease_title_dominant columns')
parser.add_argument('-s', dest='sample_id', help='Sample ID')
parser.add_argument('--exclude-cols', dest='exclude_cols', help='[DEPRECATED AS OF 4/1/2025, TODO: REMOVE] Columns to exclude when calculating duplicate rows to drop')
parser.add_argument('--cols-for-varkey', dest='cols_for_varkey', help='Columns to use to create unique string for each row')
parser.add_argument('--float-cols', dest='float_cols', help='[DEPRECATED AS OF 4/1/2025, TODO: REMOVE] Columns to convert from float to int to str for uniform formatting across inputs')
parser.add_argument('--priority-cols', dest='priority_cols', help='Columns to prioritize/put at front of output')
parser.add_argument('--cols-to-rename', dest='cols_to_rename', help='TSV with columns to rename after removing vep.transcript_consequences. and info. prefixes')
parser.add_argument('--ff-estimate', dest='xgenotyping_nomat_fetal_fraction_estimate', help='Fetal fraction estimate')
parser.add_argument('--hpo-uri', dest='sample_hpo_uri', help='Path to file with HPO terms for each sample')
parser.add_argument('--hpo-col', dest='hpo_col', help='Column in HPO file to annotate with')

args = parser.parse_args()
input_uris = args.input_uris.split(',')
exclude_cols = args.exclude_cols.split(',')    
cols_for_varkey = args.cols_for_varkey.split(',')
float_cols = args.float_cols.split(',')
priority_cols = args.priority_cols.split(';')
prefix = args.prefix
pheno_uri = args.gene_phenotype_map
cols_to_rename = pd.read_csv(args.cols_to_rename, sep='\t', header=None, names=['old_name', 'new_name'])\
    .set_index('old_name').new_name.to_dict()  # dict mapping old name to new name
xgenotyping_nomat_fetal_fraction_estimate = float(args.xgenotyping_nomat_fetal_fraction_estimate)
sample_id = args.sample_id
sample_hpo_uri = args.sample_hpo_uri
hpo_col = args.hpo_col

# Fix float formatting before merging variant_category column
def convert_to_uniform_format(num):
    '''
    To convert e.g. 1384.0 to 1384 while keeping in mind formatting like '1384-1403'
    '''
    if pd.isna(num):
        return num
    try:
        as_float = float(num)
        as_int = int(as_float)
        # only convert to int and then to str if it matches the original float
        if as_int==as_float:
            return str(as_int)  
        else:
            return as_float
    except Exception as e:
        return num

merged_df = pd.DataFrame()
all_cols = []
n_inputs_to_merge = 0

for i, uri in enumerate(input_uris):
    df = pd.concat(pd.read_csv(uri, sep='\t', chunksize=100_000))
    # Skip empty df
    if df.empty:
        continue
    # Strip quotes etc. from every column
    for col in df.columns:
        if df[col].dtype=='object':
            df[col] = df[col].str.strip('\n').str.replace('\"','').str.replace('[','').str.replace(']','').replace({'': np.nan})           
            try:  # convert float column
                df[col] = df[col].astype(float)
            except:
                pass
    
    # Make unique VarKey
    df['VarKey'] = df[cols_for_varkey].astype(str).apply(':'.join, axis=1)
    # NEW 3/31/2025: run across all columns, not just select float_cols
    for col in df.columns:
        # NEW 4/2/2025: skip boolean columns
        if df[col].dtype!='bool':
            df[col] = df[col].apply(convert_to_uniform_format)
    # NEW 3/12/2025: output_category for getting unique tiers below
    df['output_category'] = df.variant_category
    # Check if variant_category already has multiple values (e.g. CompHet/XLR/hom_var/mat_carrier output)
    n_variant_categories = df['variant_category'].value_counts().index.size
    if n_variant_categories>1:
        # Merge variant_category as comma separated string
        df['variant_category'] = df.VarKey.map(df.groupby('VarKey').variant_category.unique().apply(sorted).apply(','.join).to_dict())
        # Drop duplicate rows using all columns except variant_source (can be different because of comphets) 
        all_cols_minus_variant_source = [col for col in df.columns if col!='variant_source'] 
        df = df.drop_duplicates(all_cols_minus_variant_source)

    # NEW 3/31/2025: convert Tier to string
    df['Tier'] = df.Tier.astype(str)

    all_cols += df.columns.tolist()
    n_inputs_to_merge += 1
    merged_df = pd.concat([merged_df, df])

# Merge variant_category as comma separated string for various outputs
merged_df['variant_category'] = merged_df.VarKey.map(merged_df.groupby('VarKey').variant_category.apply(','.join).to_dict())

# Temporary lists for condensing tiers
merged_df['output_category_list'] = merged_df.VarKey.map(merged_df.groupby('VarKey').output_category.apply(list).to_dict())
merged_df['Tier_List'] = merged_df.VarKey.map(merged_df.groupby('VarKey').Tier.apply(lambda lst: pd.Series(lst).dropna().astype(str).tolist()).to_dict())

recessive_substrings = ['recessive', 'XLR', 'maternal_carrier', 'hom_var']

# NEW 3/31/2025: add 'other' category
def condense_output_category_and_tier(row):
    tier_dict = {'comphet': [], 'recessive': [], 'dominant': [], 'other': []}
    
    # Iterate through the output categories and their corresponding tiers
    for output_category, tier in zip(row.output_category_list, row.Tier_List):
        if 'comphet' in output_category:
            tier_dict['comphet'].append(tier)
        if any(substring in output_category for substring in recessive_substrings):
            tier_dict['recessive'].append(tier)
        if 'dominant' in output_category:
            tier_dict['dominant'].append(tier)
        if 'other' in output_category:
            tier_dict['other'].append(tier)

    output_categories_to_return, tiers_to_return = [], []

    # Sanity check that there is only one unique tier (for recessives and dominants)
    for output_category, tiers in tier_dict.items():
        if len(set(tiers)) > 1:
            raise Exception(f"{output_category} has multiple tiers! {set(tiers)}")
        if len(tiers) > 0:
            output_categories_to_return.append(output_category)
            tiers_to_return.append(tiers[0])
    
    return ','.join(output_categories_to_return), ','.join(tiers_to_return)

# Make corresponding inheritance_mode and Tier columns as comma separated strings
merged_df[['inheritance_mode', 'Tier']] = merged_df.apply(condense_output_category_and_tier, axis=1, result_type='expand')
merged_df = merged_df.drop(['output_category','output_category_list', 'Tier_List'], axis=1).copy()

# Prioritize CompHet/XLR/hom_var/mat_carrier output because extra columns
# NEW 3/28/2025: Use n_inputs_to_merge to account for empty inputs
col_counts = pd.Series(all_cols).value_counts()
extra_cols = col_counts[col_counts<n_inputs_to_merge].index.tolist()
cols_for_duplicate = list(np.setdiff1d(merged_df.columns, extra_cols+exclude_cols))
# NEW 4/1/2025: Simplify dropping duplicates by just using select columns
cols_for_duplicate = ['Tier', 'inheritance_mode', 'VarKey']
merged_df = merged_df.drop_duplicates(cols_for_duplicate)

# Drop duplicate columns from tiering script
merged_df = merged_df.loc[:,~merged_df.columns.str.contains('\.1')]

# Remove 'info.' and 'vep.transcript_consequences.' prefixes from column names
merged_df.columns = merged_df.columns.str.replace('info.','').str.replace('vep.transcript_consequences.','')

# Drop duplicate columns after renaming
merged_df = merged_df.loc[:,~merged_df.columns.duplicated()]

# NEW 3/12/2025: Drop columns where all values are empty
merged_df = merged_df.dropna(axis=1, how='all').copy()

# NEW 4/2/2025: Remove all mother_entry columns that are identical to proband_entry
mother_entry_cols = merged_df.columns[merged_df.columns.str.contains('mother_entry')]
original_format_cols = mother_entry_cols.str.split('.').str[1]
redundant_mother_entry_cols = []
for format_col in original_format_cols:
    if (merged_df[f"mother_entry.{format_col}"]==merged_df[f"proband_entry.{format_col}"]).all():
        redundant_mother_entry_cols.append(f"mother_entry.{format_col}")
merged_df = merged_df.drop(redundant_mother_entry_cols, axis=1)

# NEW 4/2/2025: Rename columns based on cols_to_rename input
merged_df = merged_df.rename(cols_to_rename, axis=1)

# NEW 3/12/2025: Split HGVSc and HGVSp
merged_df[['HGVSc_ENST', 'HGVSc']] = merged_df['HGVSc'].str.split(':', expand=True)
merged_df[['HGVSp_ENSP', 'HGVSp']] = merged_df['HGVSp'].str.split(':', expand=True)

# Drop VarKey column before export
merged_df = merged_df.drop('VarKey', axis=1).copy()
remaining_cols = list(np.setdiff1d(merged_df.columns, priority_cols))

# Add comphet_ID column to keep comphets together when sorting
merged_df['comphet_ID'] = ''
comphet_IDs = merged_df.loc[merged_df.variant_category.str.contains('comphet'), ['id','SYMBOL']].apply(':'.join, axis=1)
merged_df.loc[merged_df.variant_category.str.contains('comphet'), 'comphet_ID'] = comphet_IDs

# Map phenotypes
pheno_df = pd.read_csv(pheno_uri, sep='\t')
merged_df['disease_title_recessive'] = merged_df.SYMBOL.map(pheno_df.set_index('gene_symbol').disease_title_recessive.to_dict())
merged_df['disease_title_dominant'] = merged_df.SYMBOL.map(pheno_df.set_index('gene_symbol').disease_title_dominant.to_dict())

# NEW 4/2/2025: Add sample fetal fraction
merged_df['Fetal_Fraction'] = xgenotyping_nomat_fetal_fraction_estimate

# NEW 4/2/2025: Add sample HPO terms
hpo_df = pd.read_csv(sample_hpo_uri, sep='\t', dtype='str').set_index('Participant')
merged_df['Case_Pheno'] = hpo_df.loc[sample_id, hpo_col]

# Add 2 empty columns as spacers after priority columns
merged_df = merged_df[priority_cols + remaining_cols].copy()
for i in range(2):
    merged_df.insert(len(priority_cols)+i, f"SPACE_{i}", np.nan)

output_filename = f"{prefix}.merged.clinical.variants.tsv"
merged_df.to_csv(output_filename, sep='\t', index=False)
