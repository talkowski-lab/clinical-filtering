version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "mergeVCFSamples.wdl" as mergeVCFSamples
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

# run on sample level
# expect only SNV/Indel VCF, no SV VCF for NIFS
workflow filterClinicalVariants {
    input {
        File vcf_file
        File empty_file  # for if include_all_maternal_carrier_variants=false

        String? confirmation_vcf  # NIFS-specific
        String? maternal_vcf  # NIFS-specific

        String sample_id  # NIFS-specific
        String? confirmation_sample_id  # NIFS-specific
        String? maternal_sample_id  # NIFS-specific

        String predicted_sex_chrom_ploidy  # XX or XY, NIFS-specific
        Float xgenotyping_nomat_fetal_fraction_estimate  # NIFS-specific

        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_filter_clinical_variants_final_v0.1.py"
        String filter_clinical_variants_snv_indel_inheritance_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_filter_clinical_variants_inheritance_final_v0.1.py"
        String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_filter_comphets_xlr_hom_var_NIFS_v0.1.py"
        String filter_final_tiers_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/tier_clinical_variants_NIFS.py"
        String add_phenotypes_merge_and_prettify_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/add_phenotypes_merge_and_prettify_clinical_variants_NIFS.py"
        String flag_from_confirmation_maternal_vcf_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/flag_clinical_variants_from_confirmation_maternal_NIFS.py"
        String helper_functions_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/ECS_small_variants_test_CLNSIGCONF/scripts/hail_clinical_helper_functions.py"

        String hail_docker
        String sv_base_mini_docker

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0.8
        Float af_threshold=0.1
        Int ac_threshold=10  
        Int ac_rec_threshold=10  
        Int ac_dom_threshold=3  
        Float af_rec_threshold=0.05  
        Float af_dom_threshold=0.05  

        Float gnomad_af_threshold=0.05
        Float am_rec_threshold=0.56
        Float am_dom_threshold=0.56
        Float mpc_rec_threshold=2
        Float mpc_dom_threshold=2
        Float gnomad_af_rec_threshold=0.001
        Float gnomad_af_dom_threshold=0.001
        Float loeuf_v2_threshold=0.35
        Float loeuf_v4_threshold=0.6

        Int ECNT_threshold=6
        Float ncount_over_proband_DP_threshold=0.05
        Int GQ_threshold=20

        String genome_build='GRCh38'

        Boolean pass_filter=false
        Boolean include_not_genCC_OMIM=false  # NIFS-specific
        Boolean include_all_maternal_carrier_variants=false

        File gene_phenotype_map
        File carrier_gene_list  # NIFS-specific, TODO: not actually NIFS-specific anymore?
        File sample_hpo_uri  # NIFS-specific
        File gene_hpo_uri  # NIFS-specific
        File hpo_id_to_name_uri  # NIFS-specific
        String hpo_id_col = 'Anomalies with HPO codes (Screening)'
        String phenotype_col = 'Anomalies on PG03 at Eligibility Screening'
        String rec_gene_list_tsv='NA'  # for filtering by gene list(s), tab-separated "gene_list_name"\t"gene_list_uri"
        String dom_gene_list_tsv='NA'

        # ALL NIFS-specific, for addPhenotypesMergeAndPrettifyOutputs task
        Array[String] cols_for_varkey=['locus','alleles','id','vep.transcript_consequences.SYMBOL','vep.transcript_consequences.Feature','vep.transcript_consequences.Consequence','vep.transcript_consequences.HGVSc']
        Array[String] priority_cols=['fam_id', 'sex', 'Fetal_Fraction', 'Case_Pheno',
                        'locus', 'alleles', 'Tier', 'inheritance_mode', 'HGVSc_symbol',
                        'Pheno_Overlapping_HPO_IDs', 'disease_title_recessive', 'disease_title_dominant',
                        'CLNSIG', 'CLNREVSTAT', 'SYMBOL', 'HGVSc', 'HGVSp', 'IMPACT', 'Consequence', 'EXON',  
                        'CANONICAL', 'MANE_PLUS_CLINICAL', 'AD_ref,AD_alt', 'proband_entry.GT', 'mother_entry.GT', 
                        'AlphaMissense', 'REVEL', 'MPC', 'spliceAI_score', 'INTRON', 'comphet_ID',
                        'gene_list', 'cohort_AC', 'cohort_AF', 'cohort_AN', 'gnomad_popmax_af', 'GAF',
                        'maternal_carrier', 'filters']
        # Rename columns in prettify step, after removing 'vep.transcript_consequences.' and 'info.' prefixes
        Map[String, String] cols_to_rename={'proband_entry.AD': 'AD_ref,AD_alt', 'am_pathogenicity': 'AlphaMissense'}
        Array[String] static_cols = ['fam_id','id','Fetal_Fraction','sex','Case_Pheno']
        Array[String] static_cols_to_combine = ['fam_id', 'sex', 'Fetal_Fraction']  # will be '/'-separated in output

        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_filter_inheritance
        RuntimeAttr? runtime_attr_filter_comphets
        RuntimeAttr? runtime_attr_filter_tiers
        RuntimeAttr? runtime_attr_merge_prettify
        RuntimeAttr? runtime_attr_flag_conf_mat_gt
    }

    call makeDummyPed {
        input:
        predicted_sex_chrom_ploidy=predicted_sex_chrom_ploidy,
        vcf_file=vcf_file,
        hail_docker=hail_docker,
        genome_build=genome_build
    }

    call filterClinicalVariants.runClinicalFiltering as runClinicalFiltering {
        input:
        vcf_file=vcf_file,
        prefix=sample_id,
        ped_uri=makeDummyPed.ped_uri,
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
        prefix=sample_id,
        ped_uri=makeDummyPed.ped_uri,
        helper_functions_script=helper_functions_script,
        filter_clinical_variants_snv_indel_inheritance_script=filter_clinical_variants_snv_indel_inheritance_script,
        hail_docker=hail_docker,
        spliceAI_threshold=spliceAI_threshold,
        ac_rec_threshold=ac_rec_threshold,
        af_rec_threshold=af_rec_threshold,
        ac_dom_threshold=ac_dom_threshold,
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

    call filterClinicalVariants.filterCompHetsXLRHomVar as filterCompHetsXLRHomVar {
        input:
            snv_indel_vcf=runClinicalFilteringInheritance.recessive_vcf,
            clinvar_vcf=runClinicalFiltering.clinvar_vcf,
            sv_vcf='NA',
            prefix=sample_id,
            ped_uri=makeDummyPed.ped_uri,
            helper_functions_script=helper_functions_script,
            filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
            genome_build=genome_build,
            hail_docker=hail_docker,
            ad_alt_threshold=ad_alt_threshold,
            carrier_gene_list=carrier_gene_list,
            runtime_attr_override=runtime_attr_filter_comphets
    }

    call filterClinicalVariants.splitByInheritance as splitClinVarByInheritance {
        input:
            input_tsv=runClinicalFiltering.clinvar_tsv,
            hail_docker=hail_docker,
            inheritance_code_col='vep.transcript_consequences.inheritance_code',
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersInheritanceOther {
        input:
            input_tsv=runClinicalFilteringInheritance.inheritance_other_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='other',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersClinVarDominant {
        input:
            input_tsv=splitClinVarByInheritance.dominant_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='dominant',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersClinVarRecessive {
        input:
            input_tsv=splitClinVarByInheritance.recessive_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersClinVarOther {
        input:
            input_tsv=splitClinVarByInheritance.other_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='other',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersDominant {
        input:
            input_tsv=runClinicalFilteringInheritance.dominant_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='dominant',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersRecessive {
        input:
            input_tsv=runClinicalFilteringInheritance.recessive_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call finalFilteringTiers as finalFilteringTiersCompHet {
        input:
            input_tsv=filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv,
            ECNT_threshold=ECNT_threshold,
            ncount_over_proband_DP_threshold=ncount_over_proband_DP_threshold,
            GQ_threshold=GQ_threshold,
            inheritance_type='recessive',
            hail_docker=hail_docker,
            filter_final_tiers_script=filter_final_tiers_script,
            runtime_attr_override=runtime_attr_filter_tiers
    }

    call filterClinicalVariants.addPhenotypesMergeAndPrettifyOutputs as addPhenotypesMergeAndPrettifyOutputs {
        input:
            input_uris=[finalFilteringTiersCompHet.filtered_tsv,  # ORDER MATTERS (CompHet output first)
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
            prefix=sample_id,
            sample_id=sample_id,
            xgenotyping_nomat_fetal_fraction_estimate=xgenotyping_nomat_fetal_fraction_estimate,
            sample_hpo_uri=sample_hpo_uri,
            gene_hpo_uri=gene_hpo_uri,
            hpo_id_to_name_uri=hpo_id_to_name_uri,
            hpo_id_col=hpo_id_col,
            phenotype_col=phenotype_col,
            hail_docker=hail_docker,
            runtime_attr_override=runtime_attr_merge_prettify
    }

    call flagFromConfirmationMaternalVCF {
        input:
            input_tsv=addPhenotypesMergeAndPrettifyOutputs.merged_output,
            confirmation_vcf=select_first([confirmation_vcf, 'NA']),
            maternal_vcf=select_first([maternal_vcf, 'NA']),
            confirmation_sample_id=select_first([confirmation_sample_id, 'NA']),
            maternal_sample_id=select_first([maternal_sample_id, 'NA']),
            hail_docker=hail_docker,
            genome_build=genome_build,
            helper_functions_script=helper_functions_script,
            flag_from_confirmation_maternal_vcf_script=flag_from_confirmation_maternal_vcf_script,
            static_cols=static_cols,
            static_cols_to_combine=static_cols_to_combine,
            runtime_attr_override=runtime_attr_flag_conf_mat_gt
    }

    output {
        File clinvar_tsv = runClinicalFiltering.clinvar_tsv
        File clinvar_vcf = runClinicalFiltering.clinvar_vcf
        File clinvar_vcf_idx = runClinicalFiltering.clinvar_vcf_idx
        File recessive_vcf = runClinicalFilteringInheritance.recessive_vcf
        File recessive_vcf_idx = runClinicalFilteringInheritance.recessive_vcf_idx
        File recessive_tsv = runClinicalFilteringInheritance.recessive_tsv  # NEW 1/17/2025
        File dominant_tsv = runClinicalFilteringInheritance.dominant_tsv
        File comphet_xlr_hom_var_mat_carrier_tsv = filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv

        # After tiering
        File final_recessive_tsv = finalFilteringTiersRecessive.filtered_tsv  # NEW 2/25/2025
        File final_dominant_tsv = finalFilteringTiersDominant.filtered_tsv  # NEW 2/10/2025
        File final_comphet_xlr_hom_var_mat_carrier_tsv = finalFilteringTiersCompHet.filtered_tsv  # NEW 2/10/2025
        File final_inheritance_other_tsv = finalFilteringTiersInheritanceOther.filtered_tsv  # NEW 3/28/2025

        File final_merged_clinical_excel = flagFromConfirmationMaternalVCF.flagged_excel
    }
}

task makeDummyPed {
    input {
        String predicted_sex_chrom_ploidy
        File vcf_file
        String hail_docker
        String genome_build

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
    String out_ped = basename(vcf_file, file_ext) + '.ped'

    command <<<
    set -eou pipefail
    cat <<EOF > make_ped.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os

    vcf_file = sys.argv[1]
    out_ped = sys.argv[2]
    genome_build = sys.argv[3]
    cores = sys.argv[4]
    mem = int(np.floor(float(sys.argv[5])))
    predicted_sex_chrom_ploidy = sys.argv[6]

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )

    mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'], 
        reference_genome=genome_build, array_elements_required=False, call_fields=[])
    samples = mt.s.collect()
    probands = [s for s in samples if '_fetal' in s]
    mothers = [s for s in samples if '_maternal' in s]
    fam_ids = [s.split('_fetal')[0] if s in probands else s.split('_maternal')[0] for s in samples]

    # check for one proband and one mother
    if len(probands) != 1:
        raise Exception(f"You have {len(probands)} proband samples!")
    if len(mothers) != 1:
        raise Exception(f"You have {len(mothers)} maternal samples!")
    proband = probands[0]
    mother = mothers[0]

    ped = pd.DataFrame({
        'family_id': fam_ids,
        'sample_id': samples,
        'paternal_id': [0 for _ in range(len(samples))],
        'maternal_id': [0 for _ in range(len(samples))],
        'sex': [0 for _ in range(len(samples))],
        'phenotype': [0 for _ in range(len(samples))],
    })
    ped.index = ped.sample_id

    # set mother's sex to 2
    ped.loc[mother, 'sex'] = 2

    # use predicted_sex_chrom_ploidy for proband sex
    if predicted_sex_chrom_ploidy == 'XX':
        ped.loc[proband, 'sex'] = 2
    elif predicted_sex_chrom_ploidy == 'XY':
        ped.loc[proband, 'sex'] = 1
    # NEW 1/17/2025: set maternal_id to mother for proband (treat as duo instead of two singletons)
    ped.loc[proband, 'maternal_id'] = mother

    ped.to_csv(out_ped, sep='\t', index=False)
    EOF

    python3 make_ped.py ~{vcf_file} ~{out_ped} ~{genome_build} ~{cpu_cores} ~{memory} ~{predicted_sex_chrom_ploidy} 
    >>>

    output {
        File ped_uri = out_ped
    }
}

task finalFilteringTiers {
    input {
        File input_tsv
        Int ECNT_threshold
        Float ncount_over_proband_DP_threshold
        Int GQ_threshold

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
    python3 tier.py -i ~{input_tsv} -p ~{prefix} \
        --ECNT-threshold ~{ECNT_threshold} --NCount-over-proband-DP-threshold ~{ncount_over_proband_DP_threshold} \
        --GQ-threshold ~{GQ_threshold} -t ~{inheritance_type}
    >>>

    output {
        File filtered_tsv = prefix + '_tiers.tsv.gz'
    }
}

task flagFromConfirmationMaternalVCF {
    input {
        File input_tsv
        String confirmation_vcf
        String maternal_vcf

        String confirmation_sample_id
        String maternal_sample_id
        String hail_docker
        String helper_functions_script
        String flag_from_confirmation_maternal_vcf_script
        String genome_build

        Array[String] static_cols
        Array[String] static_cols_to_combine

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
    curl ~{helper_functions_script} > clinical_helper_functions.py
    curl ~{flag_from_confirmation_maternal_vcf_script} > add_GT_flags.py
    
    python3 add_GT_flags.py -i ~{input_tsv} -c ~{confirmation_vcf} -m ~{maternal_vcf} -p ~{prefix} \
        --build ~{genome_build} --conf-id ~{confirmation_sample_id} --mat-id ~{maternal_sample_id} \
        --static-cols ~{static_cols} --static-cols-to-combine ~{static_cols_to_combine}
    >>>

    output {
        File flagged_excel = prefix + '.conf.mat.flag.xlsx'
    }
}
