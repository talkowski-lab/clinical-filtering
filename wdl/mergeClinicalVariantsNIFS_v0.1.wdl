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
        Array[File] omim_dominant_tsvs
        Array[File] omim_recessive_tsvs
        Array[File] omim_recessive_vcfs
        Array[File] comphet_xlr_hom_var_mat_carrier_tsvs

        # merge TSVs
        RuntimeAttr? runtime_attr_merge_clinvar
        RuntimeAttr? runtime_attr_merge_omim_dom
        RuntimeAttr? runtime_attr_merge_omim_rec
        RuntimeAttr? runtime_attr_merge_comphets
        RuntimeAttr? runtime_attr_merge_mat_carriers
        # merge VCFs
        RuntimeAttr? runtime_attr_merge_omim_rec_vcfs
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
            tsvs=omim_dominant_tsvs,
            hail_docker=hail_docker,
            input_size=size(omim_dominant_tsvs, 'GB'),
            merged_filename=cohort_prefix+'_OMIM_dominant.tsv.gz',
            runtime_attr_override=runtime_attr_merge_omim_dom
    }

    call helpers.mergeResultsPython as mergeOMIMRecessive {
        input:
            tsvs=omim_recessive_tsvs,
            hail_docker=hail_docker,
            input_size=size(omim_recessive_tsvs, 'GB'),
            merged_filename=cohort_prefix+'_OMIM_recessive.tsv.gz',
            runtime_attr_override=runtime_attr_merge_omim_rec
    }

    # mergeVCFSamples instead of mergeVCFs
    call mergeVCFSamples.mergeVCFs as mergeOMIMRecessiveVCFs {
        input:  
            vcf_files=omim_recessive_vcfs,
            sv_base_mini_docker=sv_base_mini_docker,
            output_vcf_name=cohort_prefix + '_OMIM_recessive.vcf.gz',
            runtime_attr_override=runtime_attr_merge_omim_rec_vcfs
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
            tsvs=comphet_xlr_hom_var_mat_carrier_tsvs,
            hail_docker=hail_docker,
            input_size=size(comphet_xlr_hom_var_mat_carrier_tsvs, 'GB'),
            merged_filename="~{cohort_prefix}_SNV_Indel_comp_hets_xlr_hom_var_mat_carrier.tsv.gz",
            runtime_attr_override=runtime_attr_merge_comphets
    }

    output {
        File mat_carrier_tsv = mergeMaternalCarriers.merged_tsv
        File clinvar_tsv = mergeClinVar.merged_tsv
        File clinvar_vcf = mergeClinVarVCFs.merged_vcf_file
        File clinvar_vcf_idx = mergeClinVarVCFs.merged_vcf_idx
        File omim_recessive_vcf = mergeOMIMRecessiveVCFs.merged_vcf_file
        File omim_recessive_vcf_idx = mergeOMIMRecessiveVCFs.merged_vcf_idx
        File omim_recessive_tsv = mergeOMIMRecessive.merged_tsv
        File omim_dominant_tsv = mergeOMIMDominant.merged_tsv
        File comphet_xlr_hom_var_mat_carrier_tsv = mergeCompHetsXLRHomVar.merged_tsv
    }
}