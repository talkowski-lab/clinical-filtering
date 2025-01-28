###
# Copied from filterClinicalVariantsSV_v0.1.wdl on 1/28/2025.

## CHANGE LOG:
'''
1/28/2025:
- added more comments to script
- added description of new fields from BED to VCF header
- added renaming f3 to rsid to match SV VCF
- added using rsid as key for annotating original VCF
- adjust number_of_overlapping_BPs by subtracting 1 to match VCF indexing
'''
###

import datetime
import pandas as pd
import hail as hl
import numpy as np
import sys
import os

vcf_file = sys.argv[1]
intersect_bed = sys.argv[2]
ref_bed_with_header_uri = sys.argv[3]
genome_build = sys.argv[4]
annot_name = sys.argv[5]
cores = sys.argv[6]
mem = int(np.floor(float(sys.argv[7])))

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'], 
    reference_genome=genome_build, array_elements_required=False, call_fields=[])
header = hl.get_vcf_metadata(vcf_file)

# BED columns: 
# 0-4 from SV VCF --> vcf2bed: chrom, start, end, name, svtype 
# 5-(?) from reference bed file: chrom, start, end, ..., number_of_overlapping_BPs 
overlap_bed = hl.import_table(intersect_bed, force_bgz=True, no_header=True, 
        types={f"f{i}": 'int' for i in [1,2,6,7]})  # cast start/end columns to int
overlap_bed = overlap_bed.annotate(f1=overlap_bed.f1 + 1,  # adjust for bed 0-based coordinates
                                f6=overlap_bed.f6 + 1)

fields = list(overlap_bed.row)
overlap_field = fields[-1]  # last column is always number_of_overlapping_BPs (from -wao flag in bedtools intersect)

# calculate proportion of overlap
# NEW 1/28/2025: adjust number_of_overlapping_BPs by subtracting 1 to match VCF indexing
overlap_bed = overlap_bed.annotate(sv_len=overlap_bed.f2-overlap_bed.f1, 
                    ref_len=overlap_bed.f7-overlap_bed.f6,
                    number_of_overlapping_BPs=hl.if_else(hl.int(overlap_bed[overlap_field])>0,
                                                         hl.int(overlap_bed[overlap_field])-1,
                                                         hl.int(overlap_bed[overlap_field])))
overlap_bed = overlap_bed.annotate(sv_prop=overlap_bed.number_of_overlapping_BPs / overlap_bed.sv_len, 
                    ref_prop=overlap_bed.number_of_overlapping_BPs / overlap_bed.ref_len)

# use ref_bed_with_header for annotation column/field names
# NEW 1/28/2025: added renaming f3 to rsid to match SV VCF
ref_bed_with_header = hl.import_table(ref_bed_with_header_uri)                                    
ref_bed_with_header_idx = range(5, len(fields)-1)
ref_bed_with_header_mapping = {f"f{ref_bed_with_header_idx[i]}": list(ref_bed_with_header.row)[i].lower().replace(' ', '_') 
                                for i in range(len(ref_bed_with_header_idx))} | {'sv_prop': f"{annot_name}_overlap"}  
overlap_bed = overlap_bed.rename(ref_bed_with_header_mapping | {'f3': 'rsid'})

# annotate alleles from SV VCF
# NEW 1/28/2025: added using rsid as key for annotating original VCF
overlap_bed = overlap_bed.key_by('rsid')
mt = mt.key_rows_by('rsid')

# annotate original VCF
annot_fields = list(ref_bed_with_header_mapping.values())[3:]
mt = mt.annotate_rows(info=mt.info.annotate(
    **{field: overlap_bed[mt.row_key][field] for field in annot_fields}))

# NEW 1/28/2025: added description of new fields from BED to VCF header
for field in annot_fields:
    header['info'][field] = {'Description': f"{field} from {annot_name}.", 'Number': '.', 'Type': 'String'}

# revert to locus, alleles row key for export_vcf
mt = mt.key_rows_by('locus', 'alleles')

# export annotated VCF
hl.export_vcf(mt, os.path.basename(intersect_bed).split('.bed')[0] + '.vcf.bgz', metadata=header, tabix=True)