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
workflow filterClinicalVariants {
    input {
        Array[File] annot_vcf_files

        Array[String] predicted_sex_chrom_ploidies  # XX or XY

        String cohort_prefix
        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_clinical_variants_NIFS_v0.1.py"
        String filter_clinical_variants_snv_indel_omim_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_clinical_variants_omim_NIFS_v0.1.py"
        String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_comphets_xlr_hom_var_NIFS_v0.1.py"

        String hail_docker
        String sv_base_mini_docker

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0.8
        Float af_threshold=1  # no AF filter for NIFS
        Float gnomad_af_threshold=0.05
        Float am_rec_threshold=0.56
        Float am_dom_threshold=0.56
        Float mpc_rec_threshold=2
        Float mpc_dom_threshold=2
        Float gnomad_af_rec_threshold=0.001
        Float gnomad_af_dom_threshold=0.001
        Float loeuf_v2_threshold=0.35
        Float loeuf_v4_threshold=0.6

        String genome_build='GRCh38'

        Boolean pass_filter=false
        Boolean include_not_omim=false  # NIFS-specific

        File carrier_gene_list  # NIFS-specific
        String rec_gene_list_tsv='NA'  # for filtering by gene list(s), tab-separated "gene_list_name"\t"gene_list_uri"
        String dom_gene_list_tsv='NA'

        # merge TSVs
        RuntimeAttr? runtime_attr_merge_clinvar
        RuntimeAttr? runtime_attr_merge_omim_dom
        RuntimeAttr? runtime_attr_merge_omim_rec
        RuntimeAttr? runtime_attr_merge_comphets
        RuntimeAttr? runtime_attr_merge_mat_carriers
        # merge VCFs
        RuntimeAttr? runtime_attr_merge_omim_rec_vcfs
        RuntimeAttr? runtime_attr_merge_clinvar_vcfs
        # filtering steps
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_filter_omim
        RuntimeAttr? runtime_attr_filter_comphets
    }

    scatter (pair in zip(annot_vcf_files, predicted_sex_chrom_ploidies)) {
        File vcf_file = pair.left
        String predicted_sex_chrom_ploidy = pair.right
        call filterClinicalVariantsSingleNIFS.filterClinicalVariants as filterClinicalVariantsSingleNIFS {
            input:
            predicted_sex_chrom_ploidy=predicted_sex_chrom_ploidy,
            vcf_file=vcf_file,
            filter_clinical_variants_snv_indel_script=filter_clinical_variants_snv_indel_script,
            filter_clinical_variants_snv_indel_omim_script=filter_clinical_variants_snv_indel_omim_script,
            filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            ad_alt_threshold=ad_alt_threshold,
            spliceAI_threshold=spliceAI_threshold,
            af_threshold=af_threshold,
            gnomad_af_threshold=gnomad_af_threshold,
            am_rec_threshold=am_rec_threshold,
            am_dom_threshold=am_dom_threshold,
            mpc_rec_threshold=mpc_rec_threshold,
            mpc_dom_threshold=mpc_dom_threshold,
            gnomad_af_rec_threshold=gnomad_af_rec_threshold,
            gnomad_af_dom_threshold=gnomad_af_dom_threshold,
            loeuf_v2_threshold=loeuf_v2_threshold,
            loeuf_v4_threshold=loeuf_v4_threshold,
            genome_build=genome_build,
            pass_filter=pass_filter,
            include_not_omim=include_not_omim,
            carrier_gene_list=carrier_gene_list,
            rec_gene_list_tsv=rec_gene_list_tsv,
            dom_gene_list_tsv=dom_gene_list_tsv,
            runtime_attr_filter=runtime_attr_filter,
            runtime_attr_filter_omim=runtime_attr_filter_omim,
            runtime_attr_filter_comphets=runtime_attr_filter_comphets
        }
    }

    call helpers.mergeResultsPython as mergeClinVar {
        input:
            tsvs=filterClinicalVariantsSingleNIFS.clinvar_tsv,
            hail_docker=hail_docker,
            input_size=size(filterClinicalVariantsSingleNIFS.clinvar_tsv, 'GB'),
            merged_filename=cohort_prefix+'_clinvar_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_clinvar
    }

    call helpers.mergeResultsPython as mergeMaternalCarriers {
        input:
            tsvs=filterClinicalVariantsSingleNIFS.mat_carrier_tsv,
            hail_docker=hail_docker,
            input_size=size(filterClinicalVariantsSingleNIFS.mat_carrier_tsv, 'GB'),
            merged_filename=cohort_prefix+'_mat_carrier_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_mat_carriers
    }

    call helpers.mergeResultsPython as mergeOMIMDominant {
        input:
            tsvs=filterClinicalVariantsSingleNIFS.omim_dominant_tsv,
            hail_docker=hail_docker,
            input_size=size(filterClinicalVariantsSingleNIFS.omim_dominant_tsv, 'GB'),
            merged_filename=cohort_prefix+'_OMIM_dominant.tsv.gz',
            runtime_attr_override=runtime_attr_merge_omim_dom
    }

    call helpers.mergeResultsPython as mergeOMIMRecessive {
        input:
            tsvs=filterClinicalVariantsSingleNIFS.omim_recessive_tsv,
            hail_docker=hail_docker,
            input_size=size(filterClinicalVariantsSingleNIFS.omim_recessive_tsv, 'GB'),
            merged_filename=cohort_prefix+'_OMIM_recessive.tsv.gz',
            runtime_attr_override=runtime_attr_merge_omim_rec
    }

    # mergeVCFSamples instead of mergeVCFs
    call mergeVCFSamples.mergeVCFs as mergeOMIMRecessiveVCFs {
        input:  
            vcf_files=filterClinicalVariantsSingleNIFS.omim_recessive_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            output_vcf_name=cohort_prefix + '_OMIM_recessive.vcf.gz',
            runtime_attr_override=runtime_attr_merge_omim_rec_vcfs
    }

    call mergeVCFSamples.mergeVCFs as mergeClinVarVCFs {
        input:  
            vcf_files=filterClinicalVariantsSingleNIFS.clinvar_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            output_vcf_name=cohort_prefix + '_ClinVar_variants.vcf.gz',
            runtime_attr_override=runtime_attr_merge_clinvar_vcfs
    }

    call helpers.mergeResultsPython as mergeCompHetsXLRHomVar {
        input:
            tsvs=filterClinicalVariantsSingleNIFS.comphet_xlr_hom_var_mat_carrier_tsv,
            hail_docker=hail_docker,
            input_size=size(filterClinicalVariantsSingleNIFS.comphet_xlr_hom_var_mat_carrier_tsv, 'GB'),
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