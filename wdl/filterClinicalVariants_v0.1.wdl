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

        ## SV inputs that are optional and already have defaults in the imported workflow
        # String annotate_sv_from_intersect_bed_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_annotate_sv_from_intersect_bed_v0.1.py"
        # String annotate_sv_gene_level_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_annotate_sv_gene_level_v0.1.py"
        # String filter_clinical_sv_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_clinical_sv_v0.1.py"

        # Array[String] sv_gene_fields = ["PREDICTED_BREAKEND_EXONIC","PREDICTED_COPY_GAIN","PREDICTED_DUP_PARTIAL",
        #         "PREDICTED_INTRAGENIC_EXON_DUP","PREDICTED_INTRONIC","PREDICTED_INV_SPAN","PREDICTED_LOF","PREDICTED_MSV_EXON_OVERLAP",
        #         "PREDICTED_NEAREST_TSS","PREDICTED_PARTIAL_EXON_DUP","PREDICTED_PROMOTER","PREDICTED_TSS_DUP","PREDICTED_UTR"]
        # Array[String] permissive_csq_fields = ["PREDICTED_LOF", "PREDICTED_INTRAGENIC_EXON_DUP", "PREDICTED_COPY_GAIN",
        #                  "PREDICTED_PARTIAL_EXON_DUP", "PREDICTED_DUP_PARTIAL", "PREDICTED_INTRONIC",
        #                  "PREDICTED_INV_SPAN", "PREDICTED_UTR", "PREDICTED_PROMOTER", "PREDICTED_BREAKEND_EXONIC"]
        # Array[String] restrictive_csq_fields = ["PREDICTED_LOF", "PREDICTED_INTRAGENIC_EXON_DUP", "PREDICTED_COPY_GAIN"]

        # Float bed_overlap_threshold=0.5
        # Int size_threshold=500000  # in BP
        # Float dom_af_threshold=0.001
        # Float rec_af_threshold=0.01
        # Float gnomad_af_dom_threshold=0.01
        # Float gnomad_af_rec_threshold=0.01
        # Float gnomad_popmax_af_threshold=0.05
        # String gnomad_af_field='gnomad_v4.1_sv_AF'

        ## SNV/Indel inputs that are optional and already have defaults in the imported workflow
        # Float spliceAI_threshold=0.8
        # Float af_threshold=0.1  # NIFS-specific
        # Float gnomad_af_threshold=0.05
        # Float am_rec_threshold=0.56
        # Float am_dom_threshold=0.56
        # Float mpc_rec_threshold=2
        # Float mpc_dom_threshold=2
        # Float gnomad_af_rec_threshold=0.001
        # Float gnomad_af_dom_threshold=0.001
        # Float loeuf_v2_threshold=0.35
        # Float loeuf_v4_threshold=0.6
        
        # Boolean pass_filter=false
        # Boolean include_not_inheritance=true  # NIFS-specific

        # String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_clinical_variants_v0.1.py"
        # String filter_clinical_variants_snv_indel_inheritance_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_clinical_variants_inheritance_v0.1.py"
        # String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_MD/scripts/hail_filter_comphets_xlr_hom_var_v0.1.py"
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
            # filter_clinical_variants_snv_indel_script=filter_clinical_variants_snv_indel_script,
            # filter_clinical_variants_snv_indel_inheritance_script=filter_clinical_variants_snv_indel_inheritance_script,
            # spliceAI_threshold=spliceAI_threshold,
            # af_threshold=af_threshold,
            # gnomad_af_threshold=gnomad_af_threshold,
            # am_rec_threshold=am_rec_threshold,
            # am_dom_threshold=am_dom_threshold,
            # mpc_rec_threshold=mpc_rec_threshold,
            # mpc_dom_threshold=mpc_dom_threshold,
            # gnomad_af_rec_threshold=gnomad_af_rec_threshold,
            # gnomad_af_dom_threshold=gnomad_af_dom_threshold,
            # loeuf_v2_threshold=loeuf_v2_threshold,
            # loeuf_v4_threshold=loeuf_v4_threshold,
            # families_per_chunk=families_per_chunk,
            # pass_filter=pass_filter,
            # include_not_inheritance=include_not_inheritance,
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
            # annotate_sv_from_intersect_bed_script=annotate_sv_from_intersect_bed_script,
            # annotate_sv_gene_level_script=annotate_sv_gene_level_script,
            # filter_clinical_sv_script=filter_clinical_sv_script,
            # sv_gene_fields=sv_gene_fields,
            # permissive_csq_fields=permissive_csq_fields,
            # restrictive_csq_fields=restrictive_csq_fields,
            # bed_overlap_threshold=bed_overlap_threshold,
            # size_threshold=size_threshold,
            # dom_af_threshold=dom_af_threshold,
            # rec_af_threshold=rec_af_threshold,
            # gnomad_af_dom_threshold=gnomad_af_dom_threshold,
            # gnomad_af_rec_threshold=gnomad_af_rec_threshold,
            # gnomad_popmax_af_threshold=gnomad_popmax_af_threshold,
            # gnomad_af_field=gnomad_af_field
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

        # SVs
        # File sv_pathogenic_tsv = select_first([filterClinicalVariantsSV.sv_pathogenic_tsv, empty_file])
        # File sv_genomic_disorders_tsv = select_first([filterClinicalVariantsSV.sv_genomic_disorders_tsv, empty_file])
        # File sv_large_regions_tsv = select_first([filterClinicalVariantsSV.sv_large_regions_tsv, empty_file])
        # File sv_dominant_tsv = select_first([filterClinicalVariantsSV.sv_dominant_tsv, empty_file])
        # File sv_recessive_tsv = select_first([filterClinicalVariantsSV.sv_recessive_tsv, empty_file])
        File sv_merged_clinical_tsv = select_first([filterClinicalVariantsSV.sv_merged_clinical_tsv, empty_file])
        File sv_flagged_vcf = select_first([filterClinicalVariantsSV.sv_flagged_vcf, empty_file])
        File sv_flagged_vcf_idx = select_first([filterClinicalVariantsSV.sv_flagged_vcf_idx, empty_file])

        # CompHets/XLR/HomVar
        File comphet_xlr_hom_var_mat_carrier_tsv = filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv
    }
}