version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFs.wdl" as mergeVCFs
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers
import "filterClinicalVariantsTasks_v0.1.wdl" as filterClinicalVariants
import "filterClinicalCompHets_v0.1.wdl" as filterClinicalCompHets

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
        String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_filter_comphets_xlr_hom_var_v0.1.py"
        String filter_final_tiers_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/tier_clinical_variants_SNV_Indel.py"
        String add_phenotypes_merge_and_prettify_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/add_phenotypes_merge_and_prettify_clinical_variants.py"

        String hail_docker
        String sv_base_mini_docker

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0.8
        Float af_threshold=0.1
        Int ac_threshold=10
        Int ac_rec_threshold=10  # TODO
        Int ac_dom_threshold=3
        Float af_rec_threshold=0.05
        Float af_dom_threshold=0.01
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

        File omim_uri  # All OMIM genes
        File gene_phenotype_map
        File carrier_gene_list  
        Array[String] cols_for_varkey=['locus','alleles','id','vep.transcript_consequences.SYMBOL','vep.transcript_consequences.Feature','vep.transcript_consequences.Consequence','vep.transcript_consequences.HGVSc']
        Array[String] priority_cols=['fam_id', 'id', 'sex', 'trio_status', 'ID', 'Tier', 'inheritance_mode',
                        'disease_title_dominant', 'disease_title_recessive', 'CLNSIG', 'CLNSIGCONF', 'CLNREVSTAT', 'CLNGENE',
                        'OMIM_Gene', 'SYMBOL', 'HGVSc', 'HGVSp', 'IMPACT', 'Consequence', 'EXON', 'CANONICAL_OR_MANE_PLUS_CLINICAL',
                        'AD_ref,AD_alt', 'transmission', 'mendel_code', 'proband_entry.GT', 'father_entry.GT', 'mother_entry.GT',
                        'comphet_ID', 'AlphaMissense', 'REVEL', 'MPC', 'spliceAI_score', 'INTRON', 
                        'gene_list', 'cohort_AC', 'cohort_AF', 'gnomad_popmax_af', 'maternal_carrier', 'filters']
        # Rename columns in prettify step, after removing 'vep.transcript_consequences.' and 'info.' prefixes
        Map[String, String] cols_to_rename={'proband_entry.AD': 'AD_ref,AD_alt', 'am_pathogenicity': 'AlphaMissense', 'GENEINFO': 'CLNGENE'}

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
        # split by sample
        RuntimeAttr? runtime_attr_split_output_by_sample
    }

    scatter (vcf_file in annot_vcf_files) {
        call filterClinicalVariants.runClinicalFiltering as runClinicalFiltering {
            input:
            vcf_file=vcf_file,
            ped_uri=ped_uri,
            prefix=cohort_prefix,
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
            prefix=cohort_prefix,
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

    # CompHets
    call filterClinicalCompHets.filterClinicalCompHets as filterCompHetsXLRHomVar {
        input:
            recessive_vcf=select_first([mergeInheritanceRecessiveVCFs.merged_vcf_file, 'NA']),
            clinvar_vcf=select_first([mergeClinVarVCFs.merged_vcf_file, 'NA']),
            cohort_prefix=cohort_prefix,
            sv_base_mini_docker=sv_base_mini_docker,
            ped_uri=ped_uri,
            genome_build=genome_build,
            hail_docker=hail_docker,
            carrier_gene_list=carrier_gene_list,
            ad_alt_threshold=ad_alt_threshold,
            filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script
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

    call finalFilteringTiers as finalFilteringTiersCompHet {
        input:
            input_tsv=filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call filterClinicalVariants.addPhenotypesMergeAndPrettifyOutputs as addPhenotypesMergeAndPrettifyOutputs {
        input:
            input_uris=[
                finalFilteringTiersCompHet.filtered_tsv,  # ORDER MATTERS (CompHet output first)
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
            omim_uri=omim_uri,
            add_phenotypes_merge_and_prettify_script=add_phenotypes_merge_and_prettify_script,
            prefix=cohort_prefix,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_merge_prettify
    }

    call splitMergedOutputBySample {
        input:
            input_tsv=addPhenotypesMergeAndPrettifyOutputs.merged_output,
            hail_docker=hail_docker,
            helper_functions_script=helper_functions_script,
            runtime_attr_override=runtime_attr_split_output_by_sample
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
        File comphet_xlr_hom_var_mat_carrier_tsv = filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv

        # After tiering (does not include ClinVar separate output for now)
        File final_recessive_tsv = finalFilteringTiersRecessive.filtered_tsv
        File final_dominant_tsv = finalFilteringTiersDominant.filtered_tsv
        File final_inheritance_other_tsv = finalFilteringTiersInheritanceOther.filtered_tsv
        File final_comphet_xlr_hom_var_mat_carrier_tsv = finalFilteringTiersCompHet.filtered_tsv

        # Merged and prettified
        File final_merged_clinical_tsv = splitMergedOutputBySample.final_merged_clinical_tsv
        File final_merged_clinical_excel = splitMergedOutputBySample.final_merged_clinical_excel

        # Separate excel for each sample 
        Array[File] final_clinical_sample_excels = splitMergedOutputBySample.sample_excels
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
        File filtered_tsv = prefix + '_tiers.tsv.gz'
    }
}

task splitMergedOutputBySample {
    input {
        File input_tsv
        String hail_docker
        String helper_functions_script

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
    curl ~{helper_functions_script} > clinical_helper_functions.py
    cat <<EOF > split_by_sample.py
    import pandas as pd
    import numpy as np
    import os
    import argparse
    from clinical_helper_functions import sort_final_merged_output_by_tiers

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='input_tsv', help='Input TSV to annotate with phenotypes')

    args = parser.parse_args()
    input_uri = args.input_tsv

    merged_df = pd.concat(pd.read_csv(input_uri, sep='\t', chunksize=100_000))
    # Sort by tier (in helper functions script)
    merged_df = sort_final_merged_output_by_tiers(merged_df)
    # Replace SPACE_{i} columns with empty column names (added in addPhenotypesMergeAndPrettifyOutputs task)
    space_cols = merged_df.columns[merged_df.columns.str.contains('SPACE_')].tolist()
    merged_df = merged_df.rename({col: '' for col in space_cols}, axis=1)

    # Save merged after sort_final_merged_output_by_tiers
    merged_df.to_csv(os.path.basename(input_uri).split('.tsv')[0] + '.final.tsv.gz', sep='\t', index=False)
    # Also export full Excel
    try:
        merged_df.to_excel(os.path.basename(input_uri).split('.tsv')[0] + '.final.xlsx', index=False)
    except Exception as e:
        print(str(e))
        # Empty Excel output if too large
        pd.DataFrame().to_excel(os.path.basename(input_uri).split('.tsv')[0] + '.final.xlsx', index=False)

    # Save each sample to separate Excel output
    all_samples = merged_df.id.unique()
    tot_n_samples = all_samples.size

    for i, sample_id in enumerate(all_samples):
        output_filename = f"{sample_id}.final.merged.clinical.variants.xlsx"
        print(f"Exporting Excel for sample {i+1}/{tot_n_samples}...")
        merged_df[merged_df.id==sample_id].to_excel(output_filename, index=False)
    EOF
    python3 split_by_sample.py -i ~{input_tsv}
    >>>

    output {
        File final_merged_clinical_tsv = basename(input_tsv, file_ext) + '.final.tsv.gz'
        File final_merged_clinical_excel = basename(input_tsv, file_ext) + '.final.xlsx'
        Array[File] sample_excels = glob('*final.merged.clinical.variants.xlsx')
    }
}