version 1.0

import "mergeVCFs.wdl" as mergeVCFs
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

workflow filterClinicalVariants {
    input {
        Array[File] annot_vcf_files
        File ped_uri
        File empty_file  # for if include_all_maternal_carrier_variants=false

        String cohort_prefix

        String helper_functions_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_clinical_helper_functions.py"       
        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_clinical_variants_v0.1.py"
        String filter_clinical_variants_snv_indel_inheritance_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_clinical_variants_inheritance_v0.1.py"

        String hail_docker
        String sv_base_mini_docker

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0.8
        Float af_threshold=0.1
        Int ac_threshold=10
        Int ac_rec_threshold=10  # TODO
        Int ac_dom_threshold=10
        Float af_rec_threshold=0.001
        Float af_dom_threshold=0.001
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
        Int families_per_chunk=500

        Boolean include_all_maternal_carrier_variants=false
        Boolean pass_filter=false
        Boolean include_not_genCC_OMIM=true  # NIFS-specific

        String rec_gene_list_tsv='NA'  # for filtering by gene list(s), tab-separated "gene_list_name"\t"gene_list_uri"
        String dom_gene_list_tsv='NA' 

        Boolean sort_after_merge=false

        # merge TSVs
        RuntimeAttr? runtime_attr_merge_clinvar
        RuntimeAttr? runtime_attr_merge_inheritance_dom
        RuntimeAttr? runtime_attr_merge_inheritance_rec
        RuntimeAttr? runtime_attr_merge_inheritance_other
        RuntimeAttr? runtime_attr_merge_mat_carriers
        # merge VCFs
        RuntimeAttr? runtime_attr_merge_inheritance_rec_vcfs
        RuntimeAttr? runtime_attr_merge_clinvar_vcfs
        # filtering steps
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_filter_inheritance
    }

    scatter (vcf_file in annot_vcf_files) {
        call filterClinicalVariants.runClinicalFiltering as runClinicalFiltering {
            input:
            vcf_file=vcf_file,
            ped_uri=ped_uri,
            helper_functions_script=helper_functions_script,
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

        call filterClinicalVariants.runClinicalFilteringInheritance as runClinicalFilteringInheritance {
            input:
            vcf_file=runClinicalFiltering.filtered_vcf,
            ped_uri=ped_uri,
            helper_functions_script=helper_functions_script,
            filter_clinical_variants_snv_indel_inheritance_script=filter_clinical_variants_snv_indel_inheritance_script,
            hail_docker=hail_docker,
            spliceAI_threshold=spliceAI_threshold,
            ac_rec_threshold=ac_rec_threshold,
            ac_dom_threshold=ac_dom_threshold,
            af_rec_threshold=af_rec_threshold,
            af_dom_threshold=af_dom_threshold,
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
            include_not_genCC_OMIM=include_not_genCC_OMIM,
            rec_gene_list_tsv=rec_gene_list_tsv,
            dom_gene_list_tsv=dom_gene_list_tsv,
            runtime_attr_override=runtime_attr_filter_inheritance
        }
    }   

    call helpers.mergeResultsPython as mergeClinVar {
        input:
            tsvs=runClinicalFiltering.clinvar_tsv,
            hail_docker=hail_docker,
            input_size=size(runClinicalFiltering.clinvar_tsv, 'GB'),
            merged_filename=cohort_prefix+'_clinvar_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_clinvar
    }

    if (include_all_maternal_carrier_variants) {
        call helpers.mergeResultsPython as mergeMaternalCarriers {
            input:
                tsvs=runClinicalFiltering.mat_carrier_tsv,
                hail_docker=hail_docker,
                input_size=size(runClinicalFiltering.mat_carrier_tsv, 'GB'),
                merged_filename=cohort_prefix+'_mat_carrier_variants.tsv.gz',
                runtime_attr_override=runtime_attr_merge_mat_carriers
        }
    }

    # Merge TSVs
    call helpers.mergeResultsPython as mergeInheritanceDominant {
        input:
            tsvs=runClinicalFilteringInheritance.dominant_tsv,
            hail_docker=hail_docker,
            input_size=size(runClinicalFilteringInheritance.dominant_tsv, 'GB'),
            merged_filename=cohort_prefix+'_dominant.tsv.gz',
            runtime_attr_override=runtime_attr_merge_inheritance_dom
    }

    call helpers.mergeResultsPython as mergeInheritanceRecessive {
        input:
            tsvs=runClinicalFilteringInheritance.recessive_tsv,
            hail_docker=hail_docker,
            input_size=size(runClinicalFilteringInheritance.recessive_tsv, 'GB'),
            merged_filename=cohort_prefix+'_recessive.tsv.gz',
            runtime_attr_override=runtime_attr_merge_inheritance_rec
    }

    call helpers.mergeResultsPython as mergeInheritanceOther {
        input:
            tsvs=runClinicalFilteringInheritance.inheritance_other_tsv,
            hail_docker=hail_docker,
            input_size=size(runClinicalFilteringInheritance.inheritance_other_tsv, 'GB'),
            merged_filename=cohort_prefix+'_inheritance_other_variants.tsv.gz',
            runtime_attr_override=runtime_attr_merge_inheritance_other
    }

    # Merge VCFs
    call mergeVCFs.mergeVCFs as mergeInheritanceRecessiveVCFs {
        input:  
            vcf_files=runClinicalFilteringInheritance.recessive_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + '_recessive',
            sort_after_merge=sort_after_merge,
            runtime_attr_override=runtime_attr_merge_inheritance_rec_vcfs
    }

    call mergeVCFs.mergeVCFs as mergeClinVarVCFs {
        input:  
            vcf_files=runClinicalFiltering.clinvar_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix + '_ClinVar_variants',
            sort_after_merge=sort_after_merge,
            runtime_attr_override=runtime_attr_merge_clinvar_vcfs
    }

    output {
        File mat_carrier_tsv = select_first([mergeMaternalCarriers.merged_tsv, empty_file])
        File clinvar_tsv = mergeClinVar.merged_tsv
        File clinvar_vcf = mergeClinVarVCFs.merged_vcf_file
        File clinvar_vcf_idx = mergeClinVarVCFs.merged_vcf_idx
        File recessive_vcf = mergeInheritanceRecessiveVCFs.merged_vcf_file
        File recessive_vcf_idx = mergeInheritanceRecessiveVCFs.merged_vcf_idx
        File recessive_tsv = mergeInheritanceRecessive.merged_tsv
        File dominant_tsv = mergeInheritanceDominant.merged_tsv
        File inheritance_other_tsv = mergeInheritanceOther.merged_tsv
    }
}