import datetime
import pandas as pd
import numpy as np
import hail as hl
import sys
import ast
import os
import argparse
from clinical_helper_functions import sort_final_merged_output_by_tiers

parser = argparse.ArgumentParser(description='Parse arguments')
parser.add_argument('-i', dest='input_tsv', help='Input TSV to annotate with phenotypes')
parser.add_argument('-c', dest='confirmation_vcf_uri', help='confirmation_vcf')
parser.add_argument('-m', dest='maternal_vcf_uri', help='maternal_vcf')
parser.add_argument('-p', dest='prefix', help='Prefix for output filename')
parser.add_argument('--build', dest='build', help='Genome build')
parser.add_argument('--conf-id', dest='confirmation_sample_id', help='confirmation_sample_id')
parser.add_argument('--mat-id', dest='maternal_sample_id', help='maternal_sample_id')

args = parser.parse_args()
input_uri = args.input_tsv
confirmation_vcf_uri = args.confirmation_vcf_uri
maternal_vcf_uri = args.maternal_vcf_uri
prefix = args.prefix
build = args.build
confirmation_sample_id = args.confirmation_sample_id
maternal_sample_id = args.maternal_sample_id

hl.init(default_reference=build)

merged_ht = hl.import_table(input_uri)
# Annotate with temporary Hail-friendly locus/alleles fields
merged_ht = merged_ht.annotate(hail_locus=hl.parse_locus(merged_ht.locus),
                hail_alleles=hl.array(merged_ht.alleles.split(','))).key_by('hail_locus','hail_alleles')

# confirmation_vcf
if confirmation_vcf_uri!='NA' and confirmation_sample_id!='NA':
    conf_mt = hl.import_vcf(confirmation_vcf_uri, force_bgz=True, array_elements_required=False)
    # Annotate with temporary confirmation_sample_id
    merged_ht = merged_ht.annotate(confirmation_sample_id=confirmation_sample_id)
    # Annotate GT and filters from confirmation_vcf
    merged_ht = merged_ht.annotate(confirmation_GT=hl.str(conf_mt[merged_ht.key, merged_ht.confirmation_sample_id].GT),
                                   confirmation_filter=conf_mt.rows()[merged_ht.key].filters)
    # Flag if GT matches 
    merged_ht = merged_ht.annotate(
        GT_matches_confirmation_vcf = merged_ht.confirmation_GT.replace('0\|1', '0/1').replace('1\|0', '0/1')
         == merged_ht['proband_entry.GT']
    )

# maternal_vcf
if maternal_vcf_uri!='NA' and maternal_sample_id!='NA':
    mat_mt = hl.import_vcf(maternal_vcf_uri, force_bgz=True, array_elements_required=False)
    # Annotate with temporary maternal_sample_id
    merged_ht = merged_ht.annotate(maternal_sample_id=maternal_sample_id)
    # Annotate GT and filters from maternal_vcf
    merged_ht = merged_ht.annotate(maternal_GT=hl.str(mat_mt[merged_ht.key, merged_ht.maternal_sample_id].GT),
                                    maternal_filter=mat_mt.rows()[merged_ht.key].filters)
    # Flag if GT matches 
    merged_ht = merged_ht.annotate(
        GT_matches_maternal_vcf = 
        merged_ht.maternal_GT.replace('0\|1', '0/1').replace('1\|0', '0/1')
         == merged_ht['mother_entry.GT']
    )
    
# Drop temporary fields before export
tmp_fields = ['confirmation_sample_id','maternal_sample_id','hail_locus','hail_alleles']
fields_to_drop = np.intersect1d(tmp_fields, list(merged_ht.row)).tolist()
merged_df = merged_ht.key_by().drop(*fields_to_drop).to_pandas()   

# Sort by tier (in helper functions script)
merged_df = sort_final_merged_output_by_tiers(merged_df)

# Export to Excel, replace SPACE_{i} columns with empty column names (added in addPhenotypesMergeAndPrettifyOutputs task)
output_filename = f"{prefix}.conf.mat.flag.xlsx"
space_cols = merged_df.columns[merged_df.columns.str.contains('SPACE_')].tolist()
merged_df.rename({col: '' for col in space_cols}, axis=1).to_excel(output_filename, index=False)