import datetime
import pandas as pd
import hail as hl
import numpy as np
import sys
import os
import gzip

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

# --- Helper: check if intersect_bed is empty ---
def is_empty_bed(path):
    """Returns True if the (optionally gzipped) BED file has no data lines."""
    try:
        if path.endswith('.gz') or path.endswith('.bgz'):
            with gzip.open(path, 'rt') as f:
                for line in f:
                    if line.strip():
                        return False
                return True
        else:
            return os.path.getsize(path) == 0
    except Exception:
        return True  # treat unreadable as empty

# use ref_bed_with_header for annotation column/field names
ref_bed_with_header = hl.import_table(ref_bed_with_header_uri)
ref_fields = list(ref_bed_with_header.row)

# Determine annot_fields from ref_bed_with_header (columns 5 onward, excluding last overlap col)
# Mirrors the logic below but derived independently so it works in the empty-file branch too
# ref_bed_with_header has N columns; overlap_bed has N+5 (5 VCF cols prepended) + 1 overlap col
# So ref_bed_with_header_idx = range(5, 5 + len(ref_fields) - 1)  [same as original]
ref_annot_fields_raw = ref_fields  # all columns in ref BED
ref_annot_field_names = [f.lower().replace(' ', '_') for f in ref_annot_fields_raw]

if is_empty_bed(intersect_bed):
    print("intersect_bed is empty — skipping overlap annotation, filling with null arrays.")

    # annot_fields mirrors what the non-empty branch would produce:
    # columns 5+ of ref_bed renamed + the overlap field
    annot_fields = ref_annot_field_names + [f"{annot_name}_overlap"]

    # Annotate mt rows with missing arrays for each field
    mt = mt.annotate_rows(info=mt.info.annotate(
        **{field: hl.missing(hl.tarray(hl.tstr)) for field in annot_fields}))

    for field in annot_fields:
        header['info'][field] = {'Description': f"{field} from {annot_name}.", 'Number': '.', 'Type': 'String'}

else:
    # --- Original logic (unchanged) ---
    overlap_bed = hl.import_table(intersect_bed, force_bgz=True, no_header=True, 
            types={f"f{i}": 'int' for i in [1,2,6,7]})
    overlap_bed = overlap_bed.annotate(f1=overlap_bed.f1 + 1,
                                    f6=overlap_bed.f6 + 1)

    fields = list(overlap_bed.row)
    overlap_field = fields[-1]

    overlap_bed = overlap_bed.filter(hl.int(overlap_bed[overlap_field])>0)

    overlap_bed = overlap_bed.annotate(sv_len=overlap_bed.f2-overlap_bed.f1, 
                        ref_len=overlap_bed.f7-overlap_bed.f6,
                        number_of_overlapping_BPs=hl.int(overlap_bed[overlap_field])-1)
    overlap_bed = overlap_bed.annotate(sv_prop=overlap_bed.number_of_overlapping_BPs / overlap_bed.sv_len, 
                        ref_prop=overlap_bed.number_of_overlapping_BPs / overlap_bed.ref_len)

    ref_bed_with_header_idx = range(5, len(fields)-1)
    ref_bed_with_header_mapping = {f"f{ref_bed_with_header_idx[i]}": list(ref_bed_with_header.row)[i].lower().replace(' ', '_') 
                                    for i in range(len(ref_bed_with_header_idx))} | {'sv_prop': f"{annot_name}_overlap"}  
    overlap_bed = overlap_bed.rename(ref_bed_with_header_mapping | {'f3': 'rsid'})

    overlap_bed = overlap_bed.key_by('rsid')
    mt = mt.key_rows_by('rsid')

    annot_fields = list(ref_bed_with_header_mapping.values())[3:]
    overlap_bed = (overlap_bed.group_by(overlap_bed.rsid)
            .aggregate(**{field: hl.agg.collect(overlap_bed[field]) for field in annot_fields}))

    mt = mt.annotate_rows(info=mt.info.annotate(
        **{field: overlap_bed[mt.row_key][field] for field in annot_fields}))

    for field in annot_fields:
        header['info'][field] = {'Description': f"{field} from {annot_name}.", 'Number': '.', 'Type': 'String'}

# revert to locus, alleles row key for export_vcf
mt = mt.key_rows_by('locus', 'alleles')

hl.export_vcf(mt, os.path.basename(intersect_bed).split('.bed')[0] + '.vcf.bgz', metadata=header, tabix=True)