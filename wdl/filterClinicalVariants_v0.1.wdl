version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFs.wdl" as mergeVCFs
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers
import "filterClinicalVariantsSNVIndel_v0.1.wdl" as filterClinicalVariantsSNVIndel
import "filterClinicalVariantsSV_v0.1.wdl" as filterClinicalVariantsSV
import "filterClinicalCompHets_v0.1.wdl" as filterClinicalCompHets

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? gpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterClinicalVariants {
    input {
        Array[File]? annot_vcf_files  # SNV/Indels
        File? sv_vcf  # SV, NOTE: might have to run renameVCFSamples as upstream step if SV sample IDs don't match SNV/Indels!
        File? sv_vcf_idx
        File ped_uri
        File empty_file  # to output if only SVs or SNV/Indels input
        File carrier_gene_list

        # for SVs
        File annot_beds_with_header_tsv
        File gene_list_tsv
        File inheritance_uri
        File constrained_uri
        File prec_uri
        File hi_uri
        File ts_uri

        # ONLY OPTIONAL INPUTS THAT NEED TO BE PROVIDED IN THIS WORKFLOW BECAUSE THEY'RE USED BY MULTIPLE WORKFLOWS
        Int ad_alt_threshold=3
        String genome_build='GRCh38'
        String rec_gene_list_tsv='NA'  # for filtering by gene list(s), tab-separated "gene_list_name"\t"gene_list_uri"
        String dom_gene_list_tsv='NA'

        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker
        String variant_interpretation_docker  # for SVs, svtk vcf2bed

    }
    
    # SNV/Indel
    if (defined(annot_vcf_files)) {
        call filterClinicalVariantsSNVIndel.filterClinicalVariants as filterClinicalVariantsSNVIndel {
            input:
            annot_vcf_files=select_first([annot_vcf_files]),
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            hail_docker=hail_docker,
            sv_base_mini_docker=sv_base_mini_docker,
            ad_alt_threshold=ad_alt_threshold,
            genome_build=genome_build,
            rec_gene_list_tsv=rec_gene_list_tsv,
            dom_gene_list_tsv=dom_gene_list_tsv
        }
    }

    # SV
    if (defined(sv_vcf)) {
        call filterClinicalVariantsSV.filterClinicalVariantsSV as filterClinicalVariantsSV {
            input:
            vcf_file=select_first([sv_vcf]),
            vcf_idx=select_first([sv_vcf_idx]),
            ped_uri=ped_uri,
            gene_list_tsv=gene_list_tsv,
            inheritance_uri=inheritance_uri,
            constrained_uri=constrained_uri,
            prec_uri=prec_uri,
            hi_uri=hi_uri,
            ts_uri=ts_uri,
            annot_beds_with_header_tsv=annot_beds_with_header_tsv,
            cohort_prefix=cohort_prefix,
            genome_build=genome_build,
            hail_docker=hail_docker,
            variant_interpretation_docker=variant_interpretation_docker,
        }
    }

    # CompHets
    call filterClinicalCompHets.filterClinicalCompHets as filterCompHetsXLRHomVar {
        input:
            recessive_vcf=select_first([filterClinicalVariantsSNVIndel.recessive_vcf, 'NA']),
            clinvar_vcf=select_first([filterClinicalVariantsSNVIndel.clinvar_vcf, 'NA']),
            sv_flagged_vcf=select_first([filterClinicalVariantsSV.sv_flagged_vcf, 'NA']),
            ped_uri=ped_uri,
            genome_build=genome_build,
            hail_docker=hail_docker,
            carrier_gene_list=carrier_gene_list,
            ad_alt_threshold=ad_alt_threshold
            # filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
    }

    output {
        # SNV/Indels
        File clinvar_tsv = select_first([filterClinicalVariantsSNVIndel.clinvar_tsv, empty_file])
        File clinvar_vcf = select_first([filterClinicalVariantsSNVIndel.clinvar_vcf, empty_file])
        File clinvar_vcf_idx = select_first([filterClinicalVariantsSNVIndel.clinvar_vcf_idx, empty_file])
        File mat_carrier_tsv = select_first([filterClinicalVariantsSNVIndel.mat_carrier_tsv, empty_file])
        File recessive_vcf = select_first([filterClinicalVariantsSNVIndel.recessive_vcf, empty_file])
        File recessive_vcf_idx = select_first([filterClinicalVariantsSNVIndel.recessive_vcf_idx, empty_file])
        File recessive_tsv = select_first([filterClinicalVariantsSNVIndel.recessive_tsv, empty_file])
        File dominant_tsv = select_first([filterClinicalVariantsSNVIndel.dominant_tsv, empty_file])
        File inheritance_other_tsv = select_first([filterClinicalVariantsSNVIndel.inheritance_other_tsv, empty_file])

        # SVs
        File sv_merged_clinical_tsv = select_first([filterClinicalVariantsSV.sv_merged_clinical_tsv, empty_file])
        File sv_flagged_vcf = select_first([filterClinicalVariantsSV.sv_flagged_vcf, empty_file])
        File sv_flagged_vcf_idx = select_first([filterClinicalVariantsSV.sv_flagged_vcf_idx, empty_file])

        # CompHets/XLR/HomVar
        File comphet_xlr_hom_var_mat_carrier_tsv = filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv
    }
}