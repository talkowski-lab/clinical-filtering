version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "helpers.wdl" as helpers
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
        File? sv_vcf  # SV
        File? sv_vcf_idx
        File ped_uri
        File omim_uri
        File empty_file  # to output if only SVs or SNV/Indels input

        # for SVs
        File annot_beds_with_header_tsv
        Array[String] sv_gene_fields = ["PREDICTED_BREAKEND_EXONIC","PREDICTED_COPY_GAIN","PREDICTED_DUP_PARTIAL",
                "PREDICTED_INTRAGENIC_EXON_DUP","PREDICTED_INTRONIC","PREDICTED_INV_SPAN","PREDICTED_LOF","PREDICTED_MSV_EXON_OVERLAP",
                "PREDICTED_NEAREST_TSS","PREDICTED_PARTIAL_EXON_DUP","PREDICTED_PROMOTER","PREDICTED_TSS_DUP","PREDICTED_UTR"]
        Float bed_overlap_threshold=0.5
        Float gnomad_af_threshold=0.05
        Int size_threshold=500

        String cohort_prefix
        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_v0.1.py"
        String filter_clinical_variants_snv_indel_omim_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_omim_v0.1.py"
        String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_comphets_xlr_hom_var_v0.1.py"

        String hail_docker
        String sv_base_mini_docker
        String variant_interpretation_docker  # for SVs, svtk vcf2bed

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0  # NIFS-specific
        Float af_threshold=0.1  # NIFS-specific
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
        Boolean include_not_omim=true  # NIFS-specific

        String gene_list_tsv='NA'  # for filtering by gene list(s), tab-separated "gene_list_name"\t"gene_list_uri"
    }
    
    # SNV/Indel
    if (defined(annot_vcf_files)) {
        call filterClinicalVariantsSNVIndel.filterClinicalVariants as filterClinicalVariantsSNVIndel {
            input:
            annot_vcf_files=select_first([annot_vcf_files]),
            ped_uri=ped_uri,
            cohort_prefix=cohort_prefix,
            filter_clinical_variants_snv_indel_script=filter_clinical_variants_snv_indel_script,
            filter_clinical_variants_snv_indel_omim_script=filter_clinical_variants_snv_indel_omim_script,
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
            # families_per_chunk=families_per_chunk,
            pass_filter=pass_filter,
            # include_not_omim=include_not_omim,
            gene_list_tsv=gene_list_tsv
        }
    }

    # SV
    if (defined(sv_vcf)) {
        call filterClinicalVariantsSV.filterClinicalVariantsSV as filterClinicalVariantsSV {
            input:
            vcf_file=select_first([sv_vcf]),
            vcf_idx=select_first([sv_vcf_idx]),
            ped_uri=ped_uri,
            annot_beds_with_header_tsv=annot_beds_with_header_tsv,
            cohort_prefix=cohort_prefix,
            genome_build=genome_build,
            hail_docker=hail_docker,
            variant_interpretation_docker=variant_interpretation_docker,
            sv_gene_fields=sv_gene_fields,
            bed_overlap_threshold=bed_overlap_threshold,
            gnomad_af_threshold=gnomad_af_threshold,
            size_threshold=size_threshold
        }
    }

    # CompHets
    call filterClinicalCompHets.filterClinicalCompHets as filterCompHetsXLRHomVar {
        input:
            omim_recessive_vcf=select_first([filterClinicalVariantsSNVIndel.omim_recessive_vcf, 'NA']),
            clinvar_vcf=select_first([filterClinicalVariantsSNVIndel.clinvar_vcf, 'NA']),
            sv_filtered_vcf=select_first([filterClinicalVariantsSV.sv_filtered_vcf, 'NA']),
            ped_uri=ped_uri,
            omim_uri=omim_uri,
            sv_gene_fields=sv_gene_fields,
            filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
            genome_build=genome_build,
            hail_docker=hail_docker,
            ad_alt_threshold=ad_alt_threshold,
            gene_list_tsv=gene_list_tsv
    }

    output {
        # SNV/Indels
        File clinvar_tsv = select_first([filterClinicalVariantsSNVIndel.clinvar_tsv, empty_file])
        File clinvar_vcf = select_first([filterClinicalVariantsSNVIndel.clinvar_vcf, empty_file])
        File clinvar_vcf_idx = select_first([filterClinicalVariantsSNVIndel.clinvar_vcf_idx, empty_file])
        File omim_recessive_vcf = select_first([filterClinicalVariantsSNVIndel.omim_recessive_vcf, empty_file])
        File omim_recessive_vcf_idx = select_first([filterClinicalVariantsSNVIndel.omim_recessive_vcf_idx, empty_file])
        File omim_dominant_tsv = select_first([filterClinicalVariantsSNVIndel.omim_dominant_tsv, empty_file])

        # SVs
        File sv_pathogenic_tsv = select_first([filterClinicalVariantsSV.sv_pathogenic_tsv, empty_file])
        File sv_filtered_vcf = select_first([filterClinicalVariantsSV.sv_filtered_vcf, empty_file])
        File sv_filtered_vcf_idx = select_first([filterClinicalVariantsSV.sv_filtered_vcf_idx, empty_file])

        # CompHets/XLR/HomVar
        File comphet_xlr_hom_var_tsv = filterCompHetsXLRHomVar.comphet_xlr_hom_var_tsv
    }
}