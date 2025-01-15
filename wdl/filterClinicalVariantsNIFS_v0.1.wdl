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
        # File ped_uri  # make a dummy ped with all singletons for NIFS

        Array[String] predicted_sex_chrom_ploidies  # XX or XY

        String cohort_prefix
        String filter_clinical_variants_snv_indel_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_NIFS_v0.1.py"
        String filter_clinical_variants_snv_indel_omim_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_variants_omim_NIFS_v0.1.py"
        String filter_comphets_xlr_hom_var_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_comphets_xlr_hom_var_NIFS_v0.1.py"

        String hail_docker
        String sv_base_mini_docker

        Int ad_alt_threshold=3
        Float spliceAI_threshold=0.95
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

    # mergeVCFSamples instead of mergeVCFs
    call mergeVCFSamples.mergeVCFs as mergeOMIMRecessive {
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
        File omim_recessive_vcf = mergeOMIMRecessive.merged_vcf_file
        File omim_recessive_vcf_idx = mergeOMIMRecessive.merged_vcf_idx
        File omim_dominant_tsv = mergeOMIMDominant.merged_tsv
        File comphet_xlr_hom_var_mat_carrier_tsv = mergeCompHetsXLRHomVar.merged_tsv
    }
}

task runClinicalFiltering {
    input {
        File vcf_file
        File ped_uri

        String filter_clinical_variants_snv_indel_script
        String hail_docker
        String genome_build

        Float af_threshold
        Float gnomad_af_threshold
        Boolean pass_filter

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
    String prefix = basename(vcf_file, file_ext) + '_filtered'

    command {
        curl ~{filter_clinical_variants_snv_indel_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{prefix} ~{cpu_cores} ~{memory} \
            ~{ped_uri} ~{af_threshold} ~{gnomad_af_threshold} ~{genome_build} ~{pass_filter}
    }

    output {
        File clinvar = prefix + '_clinvar_variants.tsv.gz'
        File clinvar_vcf = prefix + '_clinvar_variants.vcf.bgz'
        File filtered_vcf = prefix + '_clinical.vcf.bgz'
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

    ped = pd.DataFrame({
        'family_id': fam_ids,
        'sample_id': samples,
        'paternal_id': [0 for _ in range(len(samples))],
        'maternal_id': [0 for _ in range(len(samples))],
        'sex': [0 for _ in range(len(samples))],
        'phenotype': [0 for _ in range(len(samples))],
    })
    ped.index = ped.sample_id

    # set mothers' sex to 2
    ped.loc[mothers, 'sex'] = 2

    # use predicted_sex_chrom_ploidy for proband sex
    for proband in probands:
        if predicted_sex_chrom_ploidy == 'XX':
            ped.loc[proband, 'sex'] = 2
        elif predicted_sex_chrom_ploidy == 'XY':
            ped.loc[proband, 'sex'] = 1

    ped.to_csv(out_ped, sep='\t', index=False)
    EOF

    python3 make_ped.py ~{vcf_file} ~{out_ped} ~{genome_build} ~{cpu_cores} ~{memory} ~{predicted_sex_chrom_ploidy} 
    >>>

    output {
        File ped_uri = out_ped
    }
}

task runClinicalFilteringOMIM {
    input {
        File vcf_file
        File ped_uri

        String filter_clinical_variants_snv_indel_omim_script
        String hail_docker
        String genome_build
        
        Int ad_alt_threshold
        Float spliceAI_threshold
        Float am_rec_threshold
        Float am_dom_threshold
        Float mpc_rec_threshold
        Float mpc_dom_threshold
        Float gnomad_af_rec_threshold
        Float gnomad_af_dom_threshold
        Float loeuf_v2_threshold
        Float loeuf_v4_threshold

        Boolean include_not_omim
        String rec_gene_list_tsv
        String dom_gene_list_tsv
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
    String prefix = basename(vcf_file, file_ext) + '_filtered'

    command {
        curl ~{filter_clinical_variants_snv_indel_omim_script} > filter_vcf.py
        python3 filter_vcf.py ~{vcf_file} ~{prefix} ~{cpu_cores} ~{memory} ~{ped_uri} \
            ~{am_rec_threshold} ~{am_dom_threshold} ~{mpc_rec_threshold} ~{mpc_dom_threshold} \
            ~{gnomad_af_rec_threshold} ~{gnomad_af_dom_threshold} ~{loeuf_v2_threshold} ~{loeuf_v4_threshold} \
            ~{genome_build} ~{ad_alt_threshold} ~{include_not_omim} ~{spliceAI_threshold} ~{rec_gene_list_tsv} ~{dom_gene_list_tsv}
    }

    output {
        File omim_recessive_vcf = prefix + '_OMIM_recessive.vcf.bgz'
        File omim_dominant = prefix + '_OMIM_dominant.tsv.gz'
    }
}

task filterCompHetsXLRHomVar {
    input {
        String snv_indel_vcf
        String clinvar_vcf
        String sv_vcf
        File ped_uri
        File carrier_gene_list
        String omim_uri
        String rec_gene_list_tsv

        Int ad_alt_threshold

        Array[String] sv_gene_fields
        String genome_build

        String filter_comphets_xlr_hom_var_script
        String hail_docker
        
        RuntimeAttr? runtime_attr_override
    }
    String variant_types_ = if (snv_indel_vcf!='NA') then 'SV_SNV_Indel' else 'SV'
    String variant_types = if (sv_vcf!='NA') then variant_types_ else 'SNV_Indel'
    Map[String, Array[String]] vcf_files = {'SV_SNV_Indel': [snv_indel_vcf, clinvar_vcf, sv_vcf], 'SV': [sv_vcf], 'SNV_Indel': [snv_indel_vcf, clinvar_vcf]}

    Float input_size = size(vcf_files[variant_types], 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0

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

    String vcf_file = if (variant_types=='SV') then sv_vcf else snv_indel_vcf
    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String prefix = basename(vcf_file, file_ext) + '_filtered'

    command {
        curl ~{filter_comphets_xlr_hom_var_script} > filter_vcf.py
        python3 filter_vcf.py ~{snv_indel_vcf} ~{clinvar_vcf} ~{sv_vcf} ~{ped_uri} ~{prefix} ~{omim_uri} \
            ~{sep=',' sv_gene_fields} ~{genome_build} ~{cpu_cores} ~{memory} ~{ad_alt_threshold} ~{rec_gene_list_tsv} ~{carrier_gene_list}
    }

    output {
        File comphet_xlr_hom_var_mat_carrier_tsv = glob('*_comp_hets_xlr_hom_var_mat_carrier.tsv.gz')[0]
    }
}