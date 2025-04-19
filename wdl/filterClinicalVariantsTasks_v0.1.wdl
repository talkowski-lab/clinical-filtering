version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

task runClinicalFiltering {
    input {
        File vcf_file
        File ped_uri
        File empty_file

        String prefix='NA'  # optional, for if vcf_file has a very long filename (e.g. NIFS)
        String helper_functions_script
        String filter_clinical_variants_snv_indel_script
        String hail_docker
        String genome_build

        Int ac_threshold
        Float af_threshold
        Float gnomad_af_threshold
        Boolean pass_filter
        Boolean include_all_maternal_carrier_variants

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
    String new_prefix = if prefix!='NA' then prefix else basename(vcf_file, file_ext)

    command {
        curl ~{helper_functions_script} > clinical_helper_functions.py
        curl ~{filter_clinical_variants_snv_indel_script} > filter_vcf.py
        python3 filter_vcf.py --vcf_file ~{vcf_file} --prefix ~{new_prefix} --cores ~{cpu_cores} --mem ~{memory} \
            --ped_uri ~{ped_uri} --af_threshold ~{af_threshold} --ac_threshold ~{ac_threshold} \
            --gnomad_af_threshold ~{gnomad_af_threshold} --build ~{genome_build} --pass_filter ~{pass_filter} \
            --include_all_maternal_carrier_variants ~{include_all_maternal_carrier_variants}
    }

    output {
        File mat_carrier_tsv = if include_all_maternal_carrier_variants then new_prefix + '_mat_carrier_variants.tsv.gz' else empty_file
        File clinvar_tsv = new_prefix + '_clinvar_variants.tsv.gz'
        File clinvar_vcf = new_prefix + '_clinvar_variants.vcf.bgz'
        File clinvar_vcf_idx = new_prefix + '_clinvar_variants.vcf.bgz.tbi'
        File filtered_vcf = new_prefix + '_clinical.vcf.bgz'
    }
}

task runClinicalFilteringInheritance {
    input {
        File vcf_file
        File ped_uri

        String prefix='NA'  # optional, for if vcf_file has a very long filename (e.g. NIFS)
        String helper_functions_script
        String filter_clinical_variants_snv_indel_inheritance_script
        String hail_docker
        String genome_build
        
        Int ac_rec_threshold
        Int ac_dom_threshold
        Float af_rec_threshold
        Float af_dom_threshold
        Int ad_alt_threshold
        Float spliceAI_threshold
        Float am_rec_threshold
        Float am_dom_threshold
        Float mpc_rec_threshold
        Float mpc_dom_threshold
        Float gnomad_af_rec_threshold
        Float gnomad_af_dom_threshold
        Float loeuf_v2_threshold
        Float loeuf_v4_threshold

        Boolean include_not_genCC_OMIM
        String rec_gene_list_tsv
        String dom_gene_list_tsv
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
    String new_prefix = if prefix!='NA' then prefix else basename(vcf_file, file_ext)

    command {
        curl ~{helper_functions_script} > clinical_helper_functions.py
        curl ~{filter_clinical_variants_snv_indel_inheritance_script} > filter_vcf.py
        python3 filter_vcf.py --vcf_file ~{vcf_file} --prefix ~{new_prefix} --cores ~{cpu_cores} --mem ~{memory} --ped_uri ~{ped_uri} \
            --ac_rec_threshold ~{ac_rec_threshold} --af_rec_threshold ~{af_rec_threshold} --ac_dom_threshold ~{ac_dom_threshold} --af_dom_threshold ~{af_dom_threshold} \
            --am_rec_threshold ~{am_rec_threshold} --am_dom_threshold ~{am_dom_threshold} --mpc_rec_threshold ~{mpc_rec_threshold} --mpc_dom_threshold ~{mpc_dom_threshold} \
            --gnomad_af_rec_threshold ~{gnomad_af_rec_threshold} --gnomad_af_dom_threshold ~{gnomad_af_dom_threshold} --loeuf_v2_threshold ~{loeuf_v2_threshold} --loeuf_v4_threshold ~{loeuf_v4_threshold} \
            --build ~{genome_build} --ad_alt_threshold ~{ad_alt_threshold} --include_not_genCC_OMIM ~{include_not_genCC_OMIM} --spliceAI_threshold ~{spliceAI_threshold} --rec_gene_list_tsv ~{rec_gene_list_tsv} --dom_gene_list_tsv ~{dom_gene_list_tsv}
    }

    output {
        File recessive_vcf = new_prefix + '_recessive.vcf.bgz'
        File recessive_vcf_idx = new_prefix + '_recessive.vcf.bgz.tbi'
        File dominant_tsv = new_prefix + '_dominant.tsv.gz'
        File recessive_tsv = new_prefix + '_recessive.tsv.gz'  # NEW 1/17/2025
        File inheritance_other_tsv = new_prefix + '_inheritance_other_variants.tsv.gz'
    }
}

task filterCompHetsXLRHomVar {
    input {
        String snv_indel_vcf
        String clinvar_vcf
        String sv_vcf
        File ped_uri
        File carrier_gene_list

        Int ad_alt_threshold

        String? prefix  # optional, for if vcf_file has a very long filename (e.g. NIFS)
        String genome_build

        String helper_functions_script
        String filter_comphets_xlr_hom_var_script
        String hail_docker
        
        RuntimeAttr? runtime_attr_override
    }
    String variant_types_ = if (snv_indel_vcf!='NA') then 'SV_SNV_Indel' else 'SV'
    String variant_types = if (sv_vcf!='NA') then variant_types_ else 'SNV_Indel'
    Map[String, Array[String]] vcf_files = {'SV_SNV_Indel': [snv_indel_vcf, clinvar_vcf, sv_vcf], 'SV': [sv_vcf], 'SNV_Indel': [snv_indel_vcf, clinvar_vcf]}

    Float input_size = size(vcf_files[variant_types], 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0

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

    String vcf_file = if (variant_types=='SV') then sv_vcf else snv_indel_vcf
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String prefix = select_first([prefix, basename(vcf_file, file_ext)])

    command {
        curl ~{helper_functions_script} > clinical_helper_functions.py
        curl ~{filter_comphets_xlr_hom_var_script} > filter_vcf.py
        python3 filter_vcf.py ~{snv_indel_vcf} ~{clinvar_vcf} ~{sv_vcf} ~{ped_uri} ~{prefix} ~{genome_build} \
        ~{cpu_cores} ~{memory} ~{ad_alt_threshold} ~{carrier_gene_list}
    }

    output {
        File comphet_xlr_hom_var_mat_carrier_tsv = glob('*_comp_hets_xlr_hom_var_mat_carrier.tsv.gz')[0]
    }
}

task splitByInheritance {
    input {
        File input_tsv

        String inheritance_code_col
        String hail_docker
        
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

    command <<<
    set -eou pipefail
    cat <<EOF > split_by_inheritance.py
    import datetime
    import pandas as pd
    import numpy as np
    import sys
    import ast
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='input_uri', help='Input TSV')
    parser.add_argument('-c', dest='inheritance_code_col', help='Column containing the (numeric) inheritance code (e.g. vep.transcript_consequences.inheritance_code)')
    parser.add_argument('--file-ext', dest='file_ext', help='File extension (.tsv or .tsv.gz)')

    args = parser.parse_args()
    input_uri = args.input_uri
    inheritance_code_col = args.inheritance_code_col
    file_ext = args.file_ext

    df = pd.read_csv(input_uri, sep='\t')

    rec_df = df[(df[inheritance_code_col].astype(str).str.contains('2')) | 
                (df[inheritance_code_col].astype(str).str.contains('4'))]
    dom_df = df[(df[inheritance_code_col].astype(str).str.contains('1')) | 
                (df[inheritance_code_col].astype(str).str.contains('3'))]
    other_df = df[(df[inheritance_code_col].astype(str).str.contains('5')) | 
                (df[inheritance_code_col].astype(str).str.contains('6'))]

    rec_df.loc[:, 'variant_category'] = rec_df['variant_category'] + '_recessive'
    dom_df.loc[:, 'variant_category'] = dom_df['variant_category'] + '_dominant'
    other_df.loc[:, 'variant_category'] = other_df['variant_category'] + '_other'

    rec_df.to_csv(os.path.basename(input_uri).split(file_ext)[0] + '.recessive.tsv', sep='\t', index=False)
    dom_df.to_csv(os.path.basename(input_uri).split(file_ext)[0] + '.dominant.tsv', sep='\t', index=False)
    other_df.to_csv(os.path.basename(input_uri).split(file_ext)[0] + '.other.tsv', sep='\t', index=False)
    EOF
    
    python3 split_by_inheritance.py -i ~{input_tsv} -c ~{inheritance_code_col} --file-ext ~{file_ext}
    >>>

    output {
        File recessive_tsv = basename(input_tsv, file_ext) + '.recessive.tsv'
        File dominant_tsv = basename(input_tsv, file_ext) + '.dominant.tsv'
        File other_tsv = basename(input_tsv, file_ext) + '.other.tsv'
    }
}

task addPhenotypesMergeAndPrettifyOutputs {
    input {
        Array[File] input_uris
        File gene_phenotype_map  # From GenCC, expects TSV with gene_symbol, disease_title_recessive, disease_title_dominant columns

        Array[String] cols_for_varkey  # Columns to use to create unique string for each row
        Array[String] priority_cols  # Columns to prioritize/put at front of output
        Map[String, String] cols_to_rename  # Columns to rename after removing 'vep.transcript_consequences.' and 'info.' prefixes
        String prefix

        # ALL OPTIONAL INPUTS FOR NIFS ONLY
        String sample_hpo_uri='NA'  # Maps samples to HPO IDs and phenotypes
        String gene_hpo_uri='NA'  # Maps genes to HPO IDs
        String hpo_id_to_name_uri='NA'  # Maps HPO IDs to HPO names
        String sample_id='NA'
        String hpo_id_col='NA'
        String phenotype_col='NA'
        Float xgenotyping_nomat_fetal_fraction_estimate=-1

        String add_phenotypes_merge_and_prettify_script
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
    curl ~{add_phenotypes_merge_and_prettify_script} > add_phenotypes_merge_and_prettify.py

    python3 add_phenotypes_merge_and_prettify.py -i ~{sep="," input_uris} -p ~{prefix} -g ~{gene_phenotype_map} \
        -s ~{sample_id} --ff-estimate ~{xgenotyping_nomat_fetal_fraction_estimate} \
        --sample-hpo-uri ~{sample_hpo_uri} --gene-hpo-uri ~{gene_hpo_uri} --hpo-id-to-name-uri ~{hpo_id_to_name_uri} \
        --hpo-id-col "~{hpo_id_col}" --phenotype-col "~{phenotype_col}" \
        --cols-for-varkey "~{sep=',' cols_for_varkey}" --priority-cols "~{sep=';' priority_cols}" \
        --cols-to-rename ~{write_map(cols_to_rename)}
    >>>

    output {
        File merged_output = prefix + '.merged.clinical.variants.tsv'
    }
}