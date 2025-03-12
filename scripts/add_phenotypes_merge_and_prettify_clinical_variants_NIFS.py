import datetime
import pandas as pd
import numpy as np
import sys
import ast
import os
import argparse

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='input_uris', help='Comma-separated list of all input TSVs')
parser.add_argument('-o', dest='output_filename', help='Output filename')
parser.add_argument('-p', dest='gene_phenotype_map', help='TSV with gene_symbol, disease_title_recessive, disease_title_dominant columns')
parser.add_argument('--exclude-cols', dest='exclude_cols', help='Columns to exclude when calculating duplicate rows to drop')
parser.add_argument('--cols-for-varkey', dest='cols_for_varkey', help='Columns to use to create unique string for each row')
parser.add_argument('--float-cols', dest='float_cols', help='Columns to convert from float to int to str for uniform formatting across inputs')
parser.add_argument('--priority-cols', dest='priority_cols', help='Columns to prioritize/put at front of output')

args = parser.parse_args()
input_uris = args.input_uris.split(',')
exclude_cols = args.exclude_cols.split(',')    
cols_for_varkey = args.cols_for_varkey.split(',')
float_cols = args.float_cols.split(',')
priority_cols = args.priority_cols.split(',')
output_filename = args.output_filename
pheno_uri = args.gene_phenotype_map

# Fix float formatting before merging variant_category column
def convert_to_uniform_format(num):
    '''
    To convert e.g. 1384.0 to 1384 while keeping in mind formatting like '1384-1403'
    '''
    if pd.isna(num):
        return num
    try:
        return str(int(float(num)))
    except Exception as e:
        return str(num)

merged_df = pd.DataFrame()
all_cols = []

for i, uri in enumerate(input_uris):
    df = pd.concat(pd.read_csv(uri, sep='\t', chunksize=100_000))
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
    for col in float_cols:
        df[col] = df[col].apply(convert_to_uniform_format)
    # Check if variant_category already has multiple values (e.g. CompHet/XLR/hom_var/mat_carrier output)
    n_variant_categories = df['variant_category'].value_counts().index.size
    if n_variant_categories>1:
        # Merge variant_category as comma separated string
        df['variant_category'] = df.VarKey.map(df.groupby('VarKey').variant_category.unique().apply(sorted).apply(','.join).to_dict())
        # Drop duplicate rows using all columns except variant_source (can be different because of comphets) 
        all_cols_minus_variant_source = [col for col in df.columns if col!='variant_source'] 
        df = df.drop_duplicates(all_cols_minus_variant_source)
    # NEW 3/12/2025: output_category for getting unique tiers below
    df['output_category'] = df.variant_category.str.replace(',', '_')

    all_cols += df.columns.tolist()
    merged_df = pd.concat([merged_df, df])

# Merge variant_category as comma separated string for various outputs
merged_df['variant_category'] = merged_df.VarKey.map(merged_df.groupby('VarKey').variant_category.apply(','.join).to_dict())

# Temporary lists for condensing tiers
merged_df['output_category_list'] = merged_df.VarKey.map(merged_df.groupby('VarKey').output_category.apply(list).to_dict())
merged_df['Tier_List'] = merged_df.VarKey.map(merged_df.groupby('VarKey').Tier.apply(lambda lst: pd.Series(lst).dropna().tolist()).to_dict())

recessive_substrings = ['recessive', 'XLR', 'maternal_carrier', 'hom_var']

def condense_output_category_and_tier(row):
    tier_dict = {'comphet': [], 'recessive': [], 'dominant': []}
    
    # Iterate through the output categories and their corresponding tiers
    for output_category, tier in zip(row.output_category_list, row.Tier_List):
        if 'comphet' in output_category:
            tier_dict['comphet'].append(tier)
        if any(substring in output_category for substring in recessive_substrings):
            tier_dict['recessive'].append(tier)
        if 'dominant' in output_category:
            tier_dict['dominant'].append(tier)
    
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
col_counts = pd.Series(all_cols).value_counts()
extra_cols = col_counts[col_counts<len(input_uris)].index.tolist()
cols_for_duplicate = list(np.setdiff1d(merged_df.columns, extra_cols+exclude_cols))
merged_df = merged_df.drop_duplicates(cols_for_duplicate)

# Drop duplicate columns from tiering script
merged_df = merged_df.loc[:,~merged_df.columns.str.contains('\.1')]

# Remove 'info.' and 'vep.transcript_consequences.' prefixes from column names
merged_df.columns = merged_df.columns.str.replace('info.','').str.replace('vep.transcript_consequences.','')

# Drop duplicate columns after renaming
merged_df = merged_df.loc[:,~merged_df.columns.duplicated()]

# Drop VarKey column before export
merged_df = merged_df.drop('VarKey', axis=1).copy()
remaining_cols = list(np.setdiff1d(merged_df.columns, priority_cols))

# Add comphet_ID column to keep comphets together when sorting
merged_df['comphet_ID'] = ''
comphet_IDs = merged_df.loc[merged_df.variant_category.str.contains('comphet'), ['id','SYMBOL']].apply(':'.join, axis=1)
merged_df.loc[merged_df.variant_category.str.contains('comphet'), 'comphet_ID'] = comphet_IDs

# Map phenotypes
pheno_df = pd.read_csv(pheno_uri, sep='\t')
merged_df['disease_title_recessive'] = merged_df.HGVSc_symbol.map(pheno_df.set_index('gene_symbol').disease_title_recessive.to_dict())
merged_df['disease_title_dominant'] = merged_df.HGVSc_symbol.map(pheno_df.set_index('gene_symbol').disease_title_dominant.to_dict())

# Add new phenotype columns to priority columns, before HGVSc_symbol
priority_cols = priority_cols[:priority_cols.index('HGVSc_symbol')] + ['disease_title_recessive', 'disease_title_dominant'] + priority_cols[priority_cols.index('HGVSc_symbol'):]

# Add 2 empty columns as spacers after priority columns
merged_df = merged_df[priority_cols + remaining_cols].copy()
for i in range(2):
    merged_df.insert(len(priority_cols)+i, f"SPACE_{i}", np.nan)

merged_df.to_csv(output_filename, sep='\t', index=False)
