version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "mergeVCFSamples.wdl" as mergeVCFSamples
import "helpers.wdl" as helpers
import "filterClinicalVariantsSingleNIFS_v0.1.wdl" as filterClinicalVariantsSingleNIFS

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

# run on sample-set level
# expect only SNV/Indel VCF, no SV VCF for NIFS
workflow mergeClinicalVariants {
    input {
        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker

        Array[File] clinvar_tsvs
        Array[File] clinvar_vcfs
        Array[File] mat_carrier_tsvs
        Array[File] final_dominant_tsvs
        Array[File] final_recessive_tsvs
        Array[File] recessive_vcfs
        Array[File] final_comphet_xlr_hom_var_mat_carrier_tsvs
        Array[File] final_merged_clinical_excel

        # merge TSVs
        RuntimeAttr? runtime_attr_merge_clinvar
        RuntimeAttr? runtime_attr_merge_dom
        RuntimeAttr? runtime_attr_merge_rec
        RuntimeAttr? runtime_attr_merge_comphets
        RuntimeAttr? runtime_attr_merge_mat_carriers
        # merge Excels
        RuntimeAttr? runtime_attr_merge_excels
        # merge VCFs
        RuntimeAttr? runtime_attr_merge_rec_vcfs
        RuntimeAttr? runtime_attr_merge_clinvar_vcfs
    }

    call helpers.mergeResultsPython as mergeClinVar {
        input:
            tsvs=clinvar_tsvs,
            hail_docker=hail_docker,
            input_size=size(clinvar_tsvs, 'GB'),
            merged_filename=cohort_prefix+'_clinvar_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_clinvar
    }

    call helpers.mergeResultsPython as mergeMaternalCarriers {
        input:
            tsvs=mat_carrier_tsvs,
            hail_docker=hail_docker,
            input_size=size(mat_carrier_tsvs, 'GB'),
            merged_filename=cohort_prefix+'_mat_carrier_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_mat_carriers
    }

    call helpers.mergeResultsPython as mergeOMIMDominant {
        input:
            tsvs=final_dominant_tsvs,
            hail_docker=hail_docker,
            input_size=size(final_dominant_tsvs, 'GB'),
            merged_filename=cohort_prefix+'_dominant.tsv.gz',
            runtime_attr_override=runtime_attr_merge_dom
    }

    call helpers.mergeResultsPython as mergeOMIMRecessive {
        input:
            tsvs=final_recessive_tsvs,
            hail_docker=hail_docker,
            input_size=size(final_recessive_tsvs, 'GB'),
            merged_filename=cohort_prefix+'_recessive.tsv.gz',
            runtime_attr_override=runtime_attr_merge_rec
    }

    # mergeVCFSamples instead of mergeVCFs
    call mergeVCFSamples.mergeVCFs as mergeOMIMRecessiveVCFs {
        input:  
            vcf_files=recessive_vcfs,
            sv_base_mini_docker=sv_base_mini_docker,
            output_vcf_name=cohort_prefix + '_recessive.vcf.gz',
            runtime_attr_override=runtime_attr_merge_rec_vcfs
    }

    call mergeVCFSamples.mergeVCFs as mergeClinVarVCFs {
        input:  
            vcf_files=clinvar_vcfs,
            sv_base_mini_docker=sv_base_mini_docker,
            output_vcf_name=cohort_prefix + '_ClinVar_variants.vcf.gz',
            runtime_attr_override=runtime_attr_merge_clinvar_vcfs
    }

    call helpers.mergeResultsPython as mergeCompHetsXLRHomVar {
        input:
            tsvs=final_comphet_xlr_hom_var_mat_carrier_tsvs,
            hail_docker=hail_docker,
            input_size=size(final_comphet_xlr_hom_var_mat_carrier_tsvs, 'GB'),
            merged_filename="~{cohort_prefix}_SNV_Indel_comp_hets_xlr_hom_var_mat_carrier.tsv.gz",
            runtime_attr_override=runtime_attr_merge_comphets
    }

    call mergeExcelsPython as mergeFinalClinicalExcels {
        input: 
            excel_files=final_merged_clinical_excel,
            hail_docker=hail_docker,
            input_size=size(final_merged_clinical_excel, 'GB'),
            merged_filename="~{cohort_prefix}merged_clinical_variants.tsv.gz",
            runtime_attr_override=runtime_attr_merge_excels
    }

    output {
        File mat_carrier_tsv = mergeMaternalCarriers.merged_tsv
        File clinvar_tsv = mergeClinVar.merged_tsv
        File clinvar_vcf = mergeClinVarVCFs.merged_vcf_file
        File clinvar_vcf_idx = mergeClinVarVCFs.merged_vcf_idx
        File recessive_vcf = mergeOMIMRecessiveVCFs.merged_vcf_file
        File recessive_vcf_idx = mergeOMIMRecessiveVCFs.merged_vcf_idx
        File recessive_tsv = mergeOMIMRecessive.merged_tsv
        File dominant_tsv = mergeOMIMDominant.merged_tsv
        File comphet_xlr_hom_var_mat_carrier_tsv = mergeCompHetsXLRHomVar.merged_tsv
        File final_merged_clinical_excel = mergeFinalClinicalExcels.merged_tsv
    }
}

task mergeExcelsPython {
     input {
        Array[String] excel_files
        String hail_docker
        String merged_filename
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

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
        cat <<EOF > merge_excels.py
        import pandas as pd
        import numpy as np
        import sys

        excel_files = pd.read_csv(sys.argv[1], header=None)[0].tolist()
        merged_filename = sys.argv[2]

        dfs = []
        tot = len(excel_files)
        merged_df = pd.DataFrame()
        for i, uri in enumerate(excel_files):
            if (i+1)%100==0:
                print(f"Loading excel {i+1}/{tot}...")
            df = pd.read_excel(uri)
            merged_df = pd.concat([merged_df, df])
        merged_df.to_csv(merged_filename, sep='\t', index=False)
        EOF

        python3 merge_excels.py ~{write_lines(excel_files)} ~{merged_filename} > stdout
    >>>

    output {
        File merged_tsv = merged_filename 
    }   
}