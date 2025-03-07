version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "mergeVCFSamples.wdl" as mergeVCFSamples
import "helpers.wdl" as helpers
import "filterClinicalVariantsTasks_v0.1.wdl" as filterClinicalVariants

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

# run on sample level
# expect only SNV/Indel VCF, no SV VCF for NIFS
workflow filterClinicalVariants {
    input {
        File vcf_file
        File empty_file  # for if include_all_maternal_carrier_variants=false

        String? confirmation_vcf  # NIFS-specific
        String? maternal_vcf  # NIFS-specific

        String sample_id  # NIFS-specific
        String? confirmation_sample_id  # NIFS-specific
        String? maternal_sample_id  # NIFS-specific

        String predicted_sex_chrom_ploidy  # XX or XY, NIFS-specific

        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_NIFS_v0.1.py"
        String filter_clinical_variants_snv_indel_omim_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_omim_NIFS_v0.1.py"
        String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_comphets_xlr_hom_var_NIFS_v0.1.py"
        String filter_final_tiers_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/tier_clinical_variants_NIFS.py"
        
        String hail_docker
        String sv_base_mini_docker

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0.8
        Float af_threshold=0.1
        Int ac_threshold=3  # NIFS-specific
        Float gnomad_af_threshold=0.05
        Float am_rec_threshold=0.56
        Float am_dom_threshold=0.56
        Float mpc_rec_threshold=2
        Float mpc_dom_threshold=2
        Float gnomad_af_rec_threshold=0.001
        Float gnomad_af_dom_threshold=0.001
        Float loeuf_v2_threshold=0.35
        Float loeuf_v4_threshold=0.6

        Int ECNT_threshold=6
        Float ncount_over_proband_DP_threshold=0.05
        Int GQ_threshold=20

        String genome_build='GRCh38'

        Boolean pass_filter=false
        Boolean include_not_omim=false  # NIFS-specific
        Boolean include_all_maternal_carrier_variants=true  # NIFS-specific

        File gene_phenotype_map  # NIFS-specific for now (2/27/2025)
        File carrier_gene_list  # NIFS-specific, TODO: not actually NIFS-specific anymore?
        String rec_gene_list_tsv='NA'  # for filtering by gene list(s), tab-separated "gene_list_name"\t"gene_list_uri"
        String dom_gene_list_tsv='NA'

        # ALL NIFS-specific, for addPhenotypesMergeAndPrettifyOutputs task
        Array[String] dup_exclude_cols=['info.CSQ','Tier','variant_source']
        Array[String] cols_for_varkey=['locus','alleles','id','vep.transcript_consequences.SYMBOL','vep.transcript_consequences.Feature','vep.transcript_consequences.Consequence','vep.transcript_consequences.HGVSc']
        Array[String] float_cols=['vep.transcript_consequences.cDNA_position', 'vep.transcript_consequences.CDS_position', 'vep.transcript_consequences.Protein_position']
        Array[String] priority_cols=['id', 'is_female', 'fam_id',
                        'variant_category','CLNSIG','Tier','locus', 'alleles',  # disease_title, classification_title inserted here
                        'HGVSc_symbol', 'HGVSc', 'HGVSp', 'Consequence', 'filters', 'gene_list',
                        'proband_entry.GT', 'proband_entry.AD', 'mother_entry.GT', 'mother_entry.AD', 
                        'am_pathogenicity', 'spliceAI_score', 'gnomad_popmax_af', 'cohort_AC', 'cohort_AF', 'comphet_ID']
        
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_filter_omim
        RuntimeAttr? runtime_attr_filter_comphets
        RuntimeAttr? runtime_attr_filter_tiers
        RuntimeAttr? runtime_attr_merge_prettify
        RuntimeAttr? runtime_attr_flag_conf_mat_gt
    }

    call makeDummyPed {
        input:
        predicted_sex_chrom_ploidy=predicted_sex_chrom_ploidy,
        vcf_file=vcf_file,
        hail_docker=hail_docker,
        genome_build=genome_build
    }

    call filterClinicalVariants.runClinicalFiltering as runClinicalFiltering {
        input:
        vcf_file=vcf_file,
        ped_uri=makeDummyPed.ped_uri,
        filter_clinical_variants_snv_indel_script=filter_clinical_variants_snv_indel_script,
        hail_docker=hail_docker,
        af_threshold=af_threshold,
        ac_threshold=ac_threshold,
        gnomad_af_threshold=gnomad_af_threshold,
        genome_build=genome_build,
        pass_filter=pass_filter,
        include_all_maternal_carrier_variants=include_all_maternal_carrier_variants,
        empty_file=empty_file,
        runtime_attr_override=runtime_attr_filter
    }

    call filterClinicalVariants.runClinicalFilteringOMIM as runClinicalFilteringOMIM {
        input:
        vcf_file=runClinicalFiltering.filtered_vcf,
        ped_uri=makeDummyPed.ped_uri,
        filter_clinical_variants_snv_indel_omim_script=filter_clinical_variants_snv_indel_omim_script,
        hail_docker=hail_docker,
        spliceAI_threshold=spliceAI_threshold,
        ad_alt_threshold=ad_alt_threshold,
        am_rec_threshold=am_rec_threshold,
        am_dom_threshold=am_dom_threshold,
        mpc_rec_threshold=mpc_rec_threshold,
        mpc_dom_threshold=mpc_dom_threshold,
        gnomad_af_rec_threshold=gnomad_af_rec_threshold,
        gnomad_af_dom_threshold=gnomad_af_dom_threshold,
        loeuf_v2_threshold=loeuf_v2_threshold,
        loeuf_v4_threshold=loeuf_v4_threshold,
        genome_build=genome_build,
        include_not_omim=include_not_omim,
        rec_gene_list_tsv=rec_gene_list_tsv,
        dom_gene_list_tsv=dom_gene_list_tsv,
        runtime_attr_override=runtime_attr_filter_omim
    }

    call filterClinicalVariants.filterCompHetsXLRHomVar as filterCompHetsXLRHomVar {
        input:
            snv_indel_vcf=runClinicalFilteringOMIM.recessive_vcf,
            clinvar_vcf=runClinicalFiltering.clinvar_vcf,
            sv_vcf='NA',
            ped_uri=makeDummyPed.ped_uri,
            filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
            genome_build=genome_build,
            hail_docker=hail_docker,
            ad_alt_threshold=ad_alt_threshold,
            carrier_gene_list=carrier_gene_list,
            runtime_attr_override=runtime_attr_filter_comphets
    }

    call finalFilteringTiers as finalFilteringTiersDominant {
        input:
            input_tsv=runClinicalFilteringOMIM.dominant,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='dominant',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersRecessive {
        input:
            input_tsv=runClinicalFilteringOMIM.recessive_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersCompHet {
        input:
            input_tsv=filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call addPhenotypesMergeAndPrettifyOutputs {
        input:
            input_uris=[finalFilteringTiersCompHet.filtered_tsv,  # ORDER MATTERS (CompHet output first, exclude mat_carrier_tsv)
                finalFilteringTiersRecessive.filtered_tsv,
                finalFilteringTiersDominant.filtered_tsv,
                runClinicalFiltering.clinvar],
            gene_phenotype_map=gene_phenotype_map,
            dup_exclude_cols=dup_exclude_cols,
            cols_for_varkey=cols_for_varkey,
            float_cols=float_cols,
            priority_cols=priority_cols,
            output_filename=sample_id + '.merged.clinical.variants.tsv',
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_merge_prettify
    }

    call flagFromConfirmationMaternalVCF {
        input:
            input_tsv=addPhenotypesMergeAndPrettifyOutputs.merged_output,
            confirmation_vcf=select_first([confirmation_vcf, 'NA']),
            maternal_vcf=select_first([maternal_vcf, 'NA']),
            confirmation_sample_id=select_first([confirmation_sample_id, 'NA']),
            maternal_sample_id=select_first([maternal_sample_id, 'NA']),
            hail_docker=hail_docker,
            genome_build=genome_build,
            runtime_attr_override=runtime_attr_flag_conf_mat_gt
    }

    output {
        File mat_carrier_tsv = runClinicalFiltering.mat_carrier_tsv
        File clinvar_tsv = runClinicalFiltering.clinvar
        File clinvar_vcf = runClinicalFiltering.clinvar_vcf
        File clinvar_vcf_idx = runClinicalFiltering.clinvar_vcf_idx
        File recessive_vcf = runClinicalFilteringOMIM.recessive_vcf
        File recessive_vcf_idx = runClinicalFilteringOMIM.recessive_vcf_idx
        File recessive_tsv = runClinicalFilteringOMIM.recessive_tsv  # NEW 1/17/2025
        File dominant_tsv = runClinicalFilteringOMIM.dominant
        File comphet_xlr_hom_var_mat_carrier_tsv = filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv

        File final_recessive_tsv = finalFilteringTiersRecessive.filtered_tsv  # NEW 2/25/2025
        File final_dominant_tsv = finalFilteringTiersDominant.filtered_tsv  # NEW 2/10/2025
        File final_comphet_xlr_hom_var_mat_carrier_tsv = finalFilteringTiersCompHet.filtered_tsv  # NEW 2/10/2025

        File final_merged_clinical_excel = flagFromConfirmationMaternalVCF.flagged_excel
    }
}

task makeDummyPed {
    input {
        String predicted_sex_chrom_ploidy
        File vcf_file
        String hail_docker
        String genome_build

        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String out_ped = basename(vcf_file, file_ext) + '.ped'

    command <<<
    set -eou pipefail
    cat <<EOF > make_ped.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    out_ped = sys.argv[2]
    genome_build = sys.argv[3]
    cores = sys.argv[4]
    mem = int(np.floor(float(sys.argv[5])))
    predicted_sex_chrom_ploidy = sys.argv[6]

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
    samples = mt.s.collect()
    probands = [s for s in samples if '_fetal' in s]
    mothers = [s for s in samples if '_maternal' in s]
    fam_ids = [s.split('_fetal')[0] if s in probands else s.split('_maternal')[0] for s in samples]

    # check for one proband and one mother
    if len(probands) != 1:
        raise Exception(f"You have {len(probands)} proband samples!")
    if len(mothers) != 1:
        raise Exception(f"You have {len(mothers)} maternal samples!")
    proband = probands[0]
    mother = mothers[0]

    ped = pd.DataFrame({
        'family_id': fam_ids,
        'sample_id': samples,
        'paternal_id': [0 for _ in range(len(samples))],
        'maternal_id': [0 for _ in range(len(samples))],
        'sex': [0 for _ in range(len(samples))],
        'phenotype': [0 for _ in range(len(samples))],
    })
    ped.index = ped.sample_id

    # set mother's sex to 2
    ped.loc[mother, 'sex'] = 2

    # use predicted_sex_chrom_ploidy for proband sex
    if predicted_sex_chrom_ploidy == 'XX':
        ped.loc[proband, 'sex'] = 2
    elif predicted_sex_chrom_ploidy == 'XY':
        ped.loc[proband, 'sex'] = 1
    # NEW 1/17/2025: set maternal_id to mother for proband (treat as duo instead of two singletons)
    ped.loc[proband, 'maternal_id'] = mother

    ped.to_csv(out_ped, sep='\t', index=False)
    EOF

    python3 make_ped.py ~{vcf_file} ~{out_ped} ~{genome_build} ~{cpu_cores} ~{memory} ~{predicted_sex_chrom_ploidy} 
    >>>

    output {
        File ped_uri = out_ped
    }
}

task finalFilteringTiers {
    input {
        File input_tsv
        Int ECNT_threshold
        Float ncount_over_proband_DP_threshold
        Int GQ_threshold

        String inheritance_type
        String hail_docker
        String filter_final_tiers_script
        
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(input_tsv, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String file_ext = if sub(basename(input_tsv), '.tsv.gz', '')!=basename(input_tsv) then '.tsv.gz' else '.tsv'
    String output_filename = basename(input_tsv, file_ext) + '.tiers.tsv'

    command <<<
    set -eou pipefail
    curl ~{filter_final_tiers_script} > tier.py
    python3 tier.py -i ~{input_tsv} -o ~{output_filename} \
        --ECNT-threshold ~{ECNT_threshold} --NCount-over-proband-DP-threshold ~{ncount_over_proband_DP_threshold} \
        --GQ-threshold ~{GQ_threshold} -t ~{inheritance_type}
    >>>

    output {
        File filtered_tsv = output_filename
    }
}

task addPhenotypesMergeAndPrettifyOutputs {
    input {
        Array[File] input_uris
        File gene_phenotype_map  # Expects TSV with gene_symbol, disease_title, classification_title columns

        Array[String] dup_exclude_cols  # Columns to exclude when calculating duplicate rows to drop
        Array[String] cols_for_varkey  # Columns to use to create unique string for each row
        Array[String] float_cols  # Columns to convert from float to int to str for uniform formatting across inputs
        Array[String] priority_cols  # Columns to prioritize/put at front of output
        
        String output_filename
        String hail_docker

        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(input_uris, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    set -eou pipefail
    cat <<EOF > add_phenotypes_merge_and_prettify.py
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
    parser.add_argument('-p', dest='gene_phenotype_map', help='TSV with gene_symbol, disease_title, classification_title columns')
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

        all_cols += df.columns.tolist()
        merged_df = pd.concat([merged_df, df])

    # Merge variant_category, Tier as comma separated string for various outputs
    merged_df['variant_category'] = merged_df.VarKey.map(merged_df.groupby('VarKey').variant_category.apply(','.join).to_dict())
    merged_df['Tier'] = merged_df.VarKey.map(merged_df.groupby('VarKey').Tier.apply(lambda lst: ','.join(pd.Series(lst).dropna().tolist())).to_dict())

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
    merged_df['disease_title'] = merged_df.HGVSc_symbol.map(pheno_df.set_index('gene_symbol').disease_title.to_dict())
    merged_df['classification_title'] = merged_df.HGVSc_symbol.map(pheno_df.set_index('gene_symbol').classification_title.to_dict())

    # Add new phenotype columns to priority columns, before HGVSc_symbol
    priority_cols = priority_cols[:priority_cols.index('HGVSc_symbol')] + ['disease_title', 'classification_title'] + priority_cols[priority_cols.index('HGVSc_symbol'):]

    # Add 2 empty columns as spacers after priority columns
    merged_df = merged_df[priority_cols + remaining_cols].copy()
    for i in range(2):
        merged_df.insert(len(priority_cols)+i, f"SPACE_{i}", np.nan)

    merged_df.to_csv(output_filename, sep='\t', index=False)
    EOF

    python3 add_phenotypes_merge_and_prettify.py -i ~{sep="," input_uris} -o ~{output_filename} -p ~{gene_phenotype_map} \
        --exclude-cols "~{sep=',' dup_exclude_cols}" --cols-for-varkey "~{sep=',' cols_for_varkey}" \
        --float-cols "~{sep=',' float_cols}" --priority-cols "~{sep=',' priority_cols}"
    >>>

    output {
        File merged_output = output_filename
    }
}

task flagFromConfirmationMaternalVCF {
    input {
        File input_tsv
        String confirmation_vcf
        String maternal_vcf

        String confirmation_sample_id
        String maternal_sample_id
        String hail_docker
        String genome_build

        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(input_tsv, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String output_filename = basename(input_tsv, '.tsv') + '.conf.mat.flag.xlsx'  # Excel file
    command <<<
    set -eou pipefail
    cat <<EOF > add_GT_flags.py
    import datetime
    import pandas as pd
    import numpy as np
    import hail as hl
    import sys
    import ast
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='input_tsv', help='Input TSV to annotate with phenotypes')
    parser.add_argument('-c', dest='confirmation_vcf_uri', help='confirmation_vcf')
    parser.add_argument('-m', dest='maternal_vcf_uri', help='maternal_vcf')
    parser.add_argument('-o', dest='output_filename', help='Output filename')
    parser.add_argument('--build', dest='build', help='Genome build')
    parser.add_argument('--conf-id', dest='confirmation_sample_id', help='confirmation_sample_id')
    parser.add_argument('--mat-id', dest='maternal_sample_id', help='maternal_sample_id')

    args = parser.parse_args()
    input_uri = args.input_tsv
    confirmation_vcf_uri = args.confirmation_vcf_uri
    maternal_vcf_uri = args.maternal_vcf_uri
    output_filename = args.output_filename
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
        # Annotate GT from confirmation_vcf
        merged_ht = merged_ht.annotate(confirmation_GT=conf_mt[merged_ht.key, merged_ht.confirmation_sample_id].GT)
        # Flag if GT matches 
        merged_ht = merged_ht.annotate(GT_matches_confirmation_vcf=merged_ht.confirmation_GT==hl.parse_call(merged_ht['proband_entry.GT']))

    # maternal_vcf
    if maternal_vcf_uri!='NA' and maternal_sample_id!='NA':
        conf_mt = hl.import_vcf(maternal_vcf_uri, force_bgz=True, array_elements_required=False)
        # Annotate with temporary maternal_sample_id
        merged_ht = merged_ht.annotate(maternal_sample_id=maternal_sample_id)
        # Annotate GT from maternal_vcf
        merged_ht = merged_ht.annotate(maternal_GT=conf_mt[merged_ht.key, merged_ht.maternal_sample_id].GT)
        # Flag if GT matches 
        merged_ht = merged_ht.annotate(GT_matches_maternal_vcf=merged_ht.maternal_GT==hl.parse_call(merged_ht['mother_entry.GT']))

    # Drop temporary fields before export
    tmp_fields = ['confirmation_sample_id','maternal_sample_id','hail_locus','hail_alleles']
    fields_to_drop = np.intersect1d(tmp_fields, list(merged_ht.row)).tolist()
    merged_df = merged_ht.key_by().drop(*fields_to_drop).to_pandas()   

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

    tmp_tier_cols = ['tiers_list','numeric_tiers_list','top_numeric_tier_comphet','top_numeric_tier','top_tier_len','tier_filter']
    merged_df = merged_df.sort_values(['top_numeric_tier_comphet','comphet_ID','top_numeric_tier','top_tier_len']).drop(tmp_tier_cols, axis=1)

    # Strip out leading commas
    for col in merged_df.columns:
        if merged_df[col].dtype=='object':
            merged_df[col] = merged_df[col].str.lstrip(',')

    # Export to Excel, replace SPACE_{i} columns with empty column names (added in addPhenotypesMergeAndPrettifyOutputs task)
    space_cols = merged_df.columns[merged_df.columns.str.contains('SPACE_')].tolist()
    merged_df.rename({col: '' for col in space_cols}, axis=1).to_excel(output_filename, index=False)
    EOF
    
    python3 add_GT_flags.py -i ~{input_tsv} -c ~{confirmation_vcf} -m ~{maternal_vcf} -o ~{output_filename} \
        --build ~{genome_build} --conf-id ~{confirmation_sample_id} --mat-id ~{maternal_sample_id}
    >>>

    output {
        File flagged_excel = output_filename
    }

}