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

        String helper_functions_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_clinical_helper_functions.py"       
        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_filter_clinical_variants_final_v0.1.py"
        String filter_clinical_variants_snv_indel_inheritance_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_filter_clinical_variants_inheritance_final_v0.1.py"
        String filter_final_tiers_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/tier_clinical_variants_SNV_Indel.py"
        String add_phenotypes_merge_and_prettify_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/add_phenotypes_merge_and_prettify_clinical_variants_NIFS.py"

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
        Boolean merge_first_pass_filtered_vcfs=false

        File gene_phenotype_map
        Array[String] cols_for_varkey=['locus','alleles','id','vep.transcript_consequences.SYMBOL','vep.transcript_consequences.Feature','vep.transcript_consequences.Consequence','vep.transcript_consequences.HGVSc']
        Array[String] priority_cols=['fam_id', 'id', 'sex', 'locus', 'alleles', 'Tier', 'inheritance_mode',
                        'disease_title_recessive', 'disease_title_dominant', 'CLNSIG', 'CLNREVSTAT', 
                        'SYMBOL', 'HGVSc', 'HGVSp', 'IMPACT', 'Consequence', 'EXON', 'CANONICAL', 'MANE_PLUS_CLINICAL',
                        'AD_ref,AD_alt', 'proband_entry.GT', 'father_entry.GT', 'mother_entry.GT',
                        'AlphaMissense', 'REVEL', 'MPC', 'spliceAI_score', 'INTRON', 'comphet_ID',
                        'gene_list', 'cohort_AC', 'cohort_AF', 'gnomad_popmax_af', 'maternal_carrier', 'filters']
        # Rename columns in prettify step, after removing 'vep.transcript_consequences.' and 'info.' prefixes
        Map[String, String] cols_to_rename={'proband_entry.AD': 'AD_ref,AD_alt', 'am_pathogenicity': 'AlphaMissense'}

        # merge TSVs
        RuntimeAttr? runtime_attr_merge_clinvar
        RuntimeAttr? runtime_attr_merge_inheritance_dom
        RuntimeAttr? runtime_attr_merge_inheritance_rec
        RuntimeAttr? runtime_attr_merge_inheritance_other
        RuntimeAttr? runtime_attr_merge_mat_carriers
        # merge VCFs
        RuntimeAttr? runtime_attr_merge_filtered_vcfs
        RuntimeAttr? runtime_attr_merge_inheritance_rec_vcfs
        RuntimeAttr? runtime_attr_merge_clinvar_vcfs
        # filtering steps
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_filter_inheritance
        RuntimeAttr? runtime_attr_filter_tiers
        # merge and prettify TSVs
        RuntimeAttr? runtime_attr_merge_prettify
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
    if (merge_first_pass_filtered_vcfs) {
        call mergeVCFs.mergeVCFs as mergeFirstPassFilteredVCFs {
            input:
                vcf_files=runClinicalFiltering.filtered_vcf,
                sv_base_mini_docker=sv_base_mini_docker,
                cohort_prefix=cohort_prefix + '_first_pass_filtered',
                sort_after_merge=sort_after_merge,
                runtime_attr_override=runtime_attr_merge_filtered_vcfs
        }     
    }

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

    call filterClinicalVariants.splitByInheritance as splitClinVarByInheritance {
        input:
            input_tsv=mergeClinVar.merged_tsv,
            hail_docker=hail_docker,
            inheritance_code_col='vep.transcript_consequences.inheritance_code',
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersInheritanceOther {
        input:
            input_tsv=mergeInheritanceOther.merged_tsv,
            inheritance_type='other',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersClinVarDominant {
        input:
            input_tsv=splitClinVarByInheritance.dominant_tsv,
            inheritance_type='dominant',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersClinVarRecessive {
        input:
            input_tsv=splitClinVarByInheritance.recessive_tsv,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersClinVarOther {
        input:
            input_tsv=splitClinVarByInheritance.other_tsv,
            inheritance_type='other',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersDominant {
        input:
            input_tsv=mergeInheritanceDominant.merged_tsv,
            inheritance_type='dominant',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersRecessive {
        input:
            input_tsv=mergeInheritanceRecessive.merged_tsv,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call filterClinicalVariants.addPhenotypesMergeAndPrettifyOutputs as addPhenotypesMergeAndPrettifyOutputs {
        input:
            input_uris=[
                # finalFilteringTiersCompHet.filtered_tsv,  # ORDER MATTERS (CompHet output first)
                finalFilteringTiersRecessive.filtered_tsv,
                finalFilteringTiersDominant.filtered_tsv,
                finalFilteringTiersClinVarRecessive.filtered_tsv,
                finalFilteringTiersClinVarDominant.filtered_tsv,
                finalFilteringTiersClinVarOther.filtered_tsv,
                finalFilteringTiersInheritanceOther.filtered_tsv],
            gene_phenotype_map=gene_phenotype_map,
            cols_for_varkey=cols_for_varkey,
            priority_cols=priority_cols,
            cols_to_rename=cols_to_rename,
            add_phenotypes_merge_and_prettify_script=add_phenotypes_merge_and_prettify_script,
            prefix=cohort_prefix,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_merge_prettify
    }

    output {
        File mat_carrier_tsv = select_first([mergeMaternalCarriers.merged_tsv, empty_file])
        File clinvar_tsv = mergeClinVar.merged_tsv
        File clinvar_vcf = mergeClinVarVCFs.merged_vcf_file
        File clinvar_vcf_idx = mergeClinVarVCFs.merged_vcf_idx
        File first_pass_filtered_vcf = select_first([mergeFirstPassFilteredVCFs.merged_vcf_file, empty_file])
        File first_pass_filtered_vcf_idx = select_first([mergeFirstPassFilteredVCFs.merged_vcf_idx, empty_file])
        File recessive_vcf = mergeInheritanceRecessiveVCFs.merged_vcf_file
        File recessive_vcf_idx = mergeInheritanceRecessiveVCFs.merged_vcf_idx
        File recessive_tsv = mergeInheritanceRecessive.merged_tsv
        File dominant_tsv = mergeInheritanceDominant.merged_tsv
        File inheritance_other_tsv = mergeInheritanceOther.merged_tsv
    
        # After tiering (does not include ClinVar separate output for now)
        File final_recessive_tsv = finalFilteringTiersRecessive.filtered_tsv
        File final_dominant_tsv = finalFilteringTiersDominant.filtered_tsv
        File final_inheritance_other_tsv = finalFilteringTiersInheritanceOther.filtered_tsv

        # Merged and prettified
        File final_merged_clinical_tsv = addPhenotypesMergeAndPrettifyOutputs.merged_output
    }
}

task finalFilteringTiers {
    input {
        File input_tsv

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
    String prefix = basename(input_tsv, file_ext)

    command <<<
    set -eou pipefail
    curl ~{filter_final_tiers_script} > tier.py
    python3 tier.py -i ~{input_tsv} -p ~{prefix} -t ~{inheritance_type}
    >>>

    output {
        File filtered_tsv = prefix + '_tiers.tsv'
    }
}