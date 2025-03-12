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

        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_NIFS_v0.1.py"
        String filter_clinical_variants_snv_indel_omim_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_omim_NIFS_v0.1.py"
        String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_comphets_xlr_hom_var_NIFS_v0.1.py"
        String filter_final_tiers_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/tier_clinical_variants_NIFS.py"
        String add_phenotypes_merge_and_prettify_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/add_phenotypes_merge_and_prettify_clinical_variants_NIFS.py"
        String flag_from_confirmation_maternal_vcf_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/flag_clinical_variants_from_confirmation_maternal_NIFS.py"
        String helper_functions_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_clinical_helper_functions.py"

        String hail_docker
        String sv_base_mini_docker

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0.8
        Float af_threshold=0.1
        Int ac_threshold=10  
        Int ac_rec_threshold=10  
        Int ac_dom_threshold=3  
        Float af_rec_threshold=0.1  
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
        Boolean include_not_omim=false  # NIFS-specific
        Boolean include_all_maternal_carrier_variants=false

        File gene_phenotype_map  # NIFS-specific for now (2/27/2025)
        File carrier_gene_list  # NIFS-specific, TODO: not actually NIFS-specific anymore?
        String rec_gene_list_tsv='NA'  # for filtering by gene list(s), tab-separated "gene_list_name"\t"gene_list_uri"
        String dom_gene_list_tsv='NA'

        # ALL NIFS-specific, for addPhenotypesMergeAndPrettifyOutputs task
        Array[String] dup_exclude_cols=['info.CSQ','Tier','variant_source']
        Array[String] cols_for_varkey=['locus','alleles','id','vep.transcript_consequences.SYMBOL','vep.transcript_consequences.Feature','vep.transcript_consequences.Consequence','vep.transcript_consequences.HGVSc']
        Array[String] float_cols=['vep.transcript_consequences.cDNA_position', 'vep.transcript_consequences.CDS_position', 'vep.transcript_consequences.Protein_position']
        Array[String] priority_cols=['id', 'is_female', 'fam_id',
                        'Tier', 'inheritance_mode', 'CLNSIG', 'CLNREVSTAT','locus', 'alleles',  # disease_title_recessive, disease_title_dominant inserted here
                        'HGVSc_symbol', 'HGVSc', 'HGVSp', 'Consequence', 'filters', 
                        'CANONICAL', 'MANE_PLUS_CLINICAL', 'gene_list', 'maternal_carrier',
                        'proband_entry.GT', 'proband_entry.AD', 'mother_entry.GT', 'mother_entry.AD', 
                        'am_pathogenicity', 'spliceAI_score', 'gnomad_popmax_af', 'cohort_AC', 'cohort_AF', 'comphet_ID']
        
        RuntimeAttr? runtime_attr_filter
        RuntimeAttr? runtime_attr_filter_omim
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

    call filterClinicalVariants.runClinicalFilteringOMIM as runClinicalFilteringOMIM {
        input:
        vcf_file=runClinicalFiltering.filtered_vcf,
        ped_uri=makeDummyPed.ped_uri,
        helper_functions_script=helper_functions_script,
        filter_clinical_variants_snv_indel_omim_script=filter_clinical_variants_snv_indel_omim_script,
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
        include_not_omim=include_not_omim,
        rec_gene_list_tsv=rec_gene_list_tsv,
        dom_gene_list_tsv=dom_gene_list_tsv,
        runtime_attr_override=runtime_attr_filter_omim
    }

    call filterClinicalVariants.filterCompHetsXLRHomVar as filterCompHetsXLRHomVar {
        input:
            snv_indel_vcf=runClinicalFilteringOMIM.recessive_vcf,
            clinvar_vcf=runClinicalFiltering.clinvar_vcf,
            sv_vcf='NA',
            ped_uri=makeDummyPed.ped_uri,
            helper_functions_script=helper_functions_script,
            filter_comphets_xlr_hom_var_script=filter_comphets_xlr_hom_var_script,
            genome_build=genome_build,
            hail_docker=hail_docker,
            ad_alt_threshold=ad_alt_threshold,
            carrier_gene_list=carrier_gene_list,
            runtime_attr_override=runtime_attr_filter_comphets
    }

    call splitByInheritance as splitClinVarByInheritance {
        input:
            input_tsv=runClinicalFiltering.clinvar,
            hail_docker=hail_docker,
            inheritance_code_col='vep.transcript_consequences.OMIM_inheritance_code',
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

    call finalFilteringTiers as finalFilteringTiersDominant {
        input:
            input_tsv=runClinicalFilteringOMIM.dominant,
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
            input_tsv=runClinicalFilteringOMIM.recessive_tsv,
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

    call addPhenotypesMergeAndPrettifyOutputs {
        input:
            input_uris=[finalFilteringTiersCompHet.filtered_tsv,  # ORDER MATTERS (CompHet output first)
                finalFilteringTiersRecessive.filtered_tsv,
                finalFilteringTiersDominant.filtered_tsv,
                finalFilteringTiersClinVarRecessive.filtered_tsv,
                finalFilteringTiersClinVarDominant.filtered_tsv],
            gene_phenotype_map=gene_phenotype_map,
            dup_exclude_cols=dup_exclude_cols,
            cols_for_varkey=cols_for_varkey,
            float_cols=float_cols,
            priority_cols=priority_cols,
            add_phenotypes_merge_and_prettify_script=add_phenotypes_merge_and_prettify_script,
            output_filename=sample_id + '.merged.clinical.variants.tsv',
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
            flag_from_confirmation_maternal_vcf_script=flag_from_confirmation_maternal_vcf_script,
            runtime_attr_override=runtime_attr_flag_conf_mat_gt
    }

    output {
        File clinvar_tsv = runClinicalFiltering.clinvar
        File clinvar_vcf = runClinicalFiltering.clinvar_vcf
        File clinvar_vcf_idx = runClinicalFiltering.clinvar_vcf_idx
        File recessive_vcf = runClinicalFilteringOMIM.recessive_vcf
        File recessive_vcf_idx = runClinicalFilteringOMIM.recessive_vcf_idx
        File recessive_tsv = runClinicalFilteringOMIM.recessive_tsv  # NEW 1/17/2025
        File dominant_tsv = runClinicalFilteringOMIM.dominant
        File comphet_xlr_hom_var_mat_carrier_tsv = filterCompHetsXLRHomVar.comphet_xlr_hom_var_mat_carrier_tsv

        File final_recessive_tsv = finalFilteringTiersRecessive.filtered_tsv  # NEW 2/25/2025
        File final_dominant_tsv = finalFilteringTiersDominant.filtered_tsv  # NEW 2/10/2025
        File final_comphet_xlr_hom_var_mat_carrier_tsv = finalFilteringTiersCompHet.filtered_tsv  # NEW 2/10/2025

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

task splitByInheritance {
    input {
        File input_tsv

        String inheritance_code_col
        String hail_docker
        
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
    cat <<EOF > split_by_inheritance.py
    import datetime
    import pandas as pd
    import numpy as np
    import sys
    import ast
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Parse arguments')
    parser.add_argument('-i', dest='input_uri', help='Input TSV')
    parser.add_argument('-c', dest='inheritance_code_col', help='Column containing the (numeric) inheritance code (e.g. vep.transcript_consequences.OMIM_inheritance_code)')
    parser.add_argument('--file-ext', dest='file_ext', help='File extension (.tsv or .tsv.gz)')

    args = parser.parse_args()
    input_uri = args.input_uri
    inheritance_code_col = args.inheritance_code_col
    file_ext = args.file_ext

    df = pd.read_csv(input_uri, sep='\t')

    rec_df = df[(df[inheritance_code_col].astype(str).str.contains('2')) | 
                (df[inheritance_code_col].astype(str).str.contains('4'))]
    dom_df = df[(df[inheritance_code_col].astype(str).str.contains('1')) | 
                (df[inheritance_code_col].astype(str).str.contains('3'))]

    rec_df.loc[:, 'variant_category'] = rec_df['variant_category'] + '_recessive'
    dom_df.loc[:, 'variant_category'] = dom_df['variant_category'] + '_dominant'

    rec_df.to_csv(os.path.basename(input_uri).split(file_ext)[0] + '.recessive.tsv', sep='\t', index=False)
    dom_df.to_csv(os.path.basename(input_uri).split(file_ext)[0] + '.dominant.tsv', sep='\t', index=False)
    EOF
    
    python3 split_by_inheritance.py -i ~{input_tsv} -c ~{inheritance_code_col} --file-ext ~{file_ext}
    >>>

    output {
        File recessive_tsv = basename(input_tsv, file_ext) + '.recessive.tsv'
        File dominant_tsv = basename(input_tsv, file_ext) + '.dominant.tsv'
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
    String output_filename = basename(input_tsv, file_ext) + '.tiers.tsv'

    command <<<
    set -eou pipefail
    curl ~{filter_final_tiers_script} > tier.py
    python3 tier.py -i ~{input_tsv} -o ~{output_filename} \
        --ECNT-threshold ~{ECNT_threshold} --NCount-over-proband-DP-threshold ~{ncount_over_proband_DP_threshold} \
        --GQ-threshold ~{GQ_threshold} -t ~{inheritance_type}
    >>>

    output {
        File filtered_tsv = output_filename
    }
}

task addPhenotypesMergeAndPrettifyOutputs {
    input {
        Array[File] input_uris
        File gene_phenotype_map  # Expects TSV with gene_symbol, disease_title_recessive, disease_title_dominant columns

        Array[String] dup_exclude_cols  # Columns to exclude when calculating duplicate rows to drop
        Array[String] cols_for_varkey  # Columns to use to create unique string for each row
        Array[String] float_cols  # Columns to convert from float to int to str for uniform formatting across inputs
        Array[String] priority_cols  # Columns to prioritize/put at front of output
        
        String add_phenotypes_merge_and_prettify_script
        String output_filename
        String hail_docker

        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(input_uris, 'GB')
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
    set -eou pipefail
    curl ~{add_phenotypes_merge_and_prettify_script} > add_phenotypes_merge_and_prettify.py

    python3 add_phenotypes_merge_and_prettify.py -i ~{sep="," input_uris} -o ~{output_filename} -p ~{gene_phenotype_map} \
        --exclude-cols "~{sep=',' dup_exclude_cols}" --cols-for-varkey "~{sep=',' cols_for_varkey}" \
        --float-cols "~{sep=',' float_cols}" --priority-cols "~{sep=',' priority_cols}"
    >>>

    output {
        File merged_output = output_filename
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
        String flag_from_confirmation_maternal_vcf_script
        String genome_build

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

    String output_filename = basename(input_tsv, '.tsv') + '.conf.mat.flag.xlsx'  # Excel file
    command <<<
    set -eou pipefail
    curl ~{flag_from_confirmation_maternal_vcf_script} > add_GT_flags.py
    
    python3 add_GT_flags.py -i ~{input_tsv} -c ~{confirmation_vcf} -m ~{maternal_vcf} -o ~{output_filename} \
        --build ~{genome_build} --conf-id ~{confirmation_sample_id} --mat-id ~{maternal_sample_id}
    >>>

    output {
        File flagged_excel = output_filename
    }

}