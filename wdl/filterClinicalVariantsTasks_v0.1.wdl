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

        String? prefix  # optional, for if vcf_file has a very long filename (e.g. NIFS)
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
    String prefix = select_first([prefix, basename(vcf_file, file_ext)])

    command {
        curl ~{helper_functions_script} > clinical_helper_functions.py
        curl ~{filter_clinical_variants_snv_indel_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{prefix} ~{cpu_cores} ~{memory} \
            ~{ped_uri} ~{af_threshold} ~{ac_threshold} ~{gnomad_af_threshold} ~{genome_build} ~{pass_filter} \
            ~{include_all_maternal_carrier_variants}
    }

    output {
        File mat_carrier_tsv = if include_all_maternal_carrier_variants then prefix + '_mat_carrier_variants.tsv.gz' else empty_file
        File clinvar_tsv = prefix + '_clinvar_variants.tsv.gz'
        File clinvar_vcf = prefix + '_clinvar_variants.vcf.bgz'
        File clinvar_vcf_idx = prefix + '_clinvar_variants.vcf.bgz.tbi'
        File filtered_vcf = prefix + '_clinical.vcf.bgz'
    }
}

task runClinicalFilteringInheritance {
    input {
        File vcf_file
        File ped_uri

        String? prefix  # optional, for if vcf_file has a very long filename (e.g. NIFS)
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
    String prefix = select_first([prefix, basename(vcf_file, file_ext)])

    command {
        curl ~{helper_functions_script} > clinical_helper_functions.py
        curl ~{filter_clinical_variants_snv_indel_inheritance_script} > filter_vcf.py
        python3 filter_vcf.py --vcf_file ~{vcf_file} --prefix ~{prefix} --cores ~{cpu_cores} --mem ~{memory} --ped_uri ~{ped_uri} \
            --ac_rec_threshold ~{ac_rec_threshold} --af_rec_threshold ~{af_rec_threshold} --ac_dom_threshold ~{ac_dom_threshold} --af_dom_threshold ~{af_dom_threshold} \
            --am_rec_threshold ~{am_rec_threshold} --am_dom_threshold ~{am_dom_threshold} --mpc_rec_threshold ~{mpc_rec_threshold} --mpc_dom_threshold ~{mpc_dom_threshold} \
            --gnomad_af_rec_threshold ~{gnomad_af_rec_threshold} --gnomad_af_dom_threshold ~{gnomad_af_dom_threshold} --loeuf_v2_threshold ~{loeuf_v2_threshold} --loeuf_v4_threshold ~{loeuf_v4_threshold} \
            --build ~{genome_build} --ad_alt_threshold ~{ad_alt_threshold} --include_not_genCC_OMIM ~{include_not_genCC_OMIM} --spliceAI_threshold ~{spliceAI_threshold} --rec_gene_list_tsv ~{rec_gene_list_tsv} --dom_gene_list_tsv ~{dom_gene_list_tsv}
    }

    output {
        File recessive_vcf = prefix + '_recessive.vcf.bgz'
        File recessive_vcf_idx = prefix + '_recessive.vcf.bgz.tbi'
        File dominant_tsv = prefix + '_dominant.tsv.gz'
        File recessive_tsv = prefix + '_recessive.tsv.gz'  # NEW 1/17/2025
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