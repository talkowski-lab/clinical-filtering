version 1.0

import "mergeVCFs.wdl" as mergeVCFs
import "helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow filterClinicalVariantsSV {
    input {
        File vcf_file  # NOTE: might have to run renameVCFSamples as upstream step if SV sample IDs don't match SNV/Indels!
        File vcf_idx
        File ped_uri

        File gene_list_tsv
        File omim_uri
        File constrained_uri
        File prec_uri
        File hi_uri
        File ts_uri

        File annot_beds_with_header_tsv  # no header, 3 columns: name, uri, match_svtype (True or False)

        String cohort_prefix
        String genome_build='GRCh38'
        String hail_docker
        String variant_interpretation_docker

        String annotate_sv_from_intersect_bed_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_annotate_sv_from_intersect_bed_v0.1.py"
        String annotate_sv_gene_level_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_annotate_sv_gene_level_v0.1.py"
        String filter_clinical_sv_script = "https://raw.githubusercontent.com/talkowski-lab/clinical-filtering/refs/heads/main/scripts/hail_filter_clinical_sv_v0.1.py"

        Array[String] permissive_csq_fields = ["PREDICTED_LOF", "PREDICTED_INTRAGENIC_EXON_DUP", "PREDICTED_COPY_GAIN",
                         "PREDICTED_PARTIAL_EXON_DUP", "PREDICTED_DUP_PARTIAL", "PREDICTED_INTRONIC",
                         "PREDICTED_INV_SPAN", "PREDICTED_UTR", "PREDICTED_PROMOTER", "PREDICTED_BREAKEND_EXONIC"]
        Array[String] restrictive_csq_fields = ["PREDICTED_LOF", "PREDICTED_INTRAGENIC_EXON_DUP", "PREDICTED_COPY_GAIN"]

        Float bed_overlap_threshold=0.5
        Int size_threshold=500000  # in BP
        Float dom_af_threshold=0.001
        Float rec_af_threshold=0.01
        Float gnomad_af_dom_threshold=0.01
        Float gnomad_af_rec_threshold=0.01
        Float gnomad_popmax_af_threshold=0.05
        String gnomad_af_field='gnomad_v4.1_sv_AF'

        RuntimeAttr? runtime_attr_bcftools
        RuntimeAttr? runtime_attr_annotate
        RuntimeAttr? runtime_attr_filter_vcf
    }

    call vcfToBed {
        input:
        vcf_file=vcf_file,
        cohort_prefix=cohort_prefix,
        variant_interpretation_docker=variant_interpretation_docker
    }

    call removeGenotypes {
        input:
        vcf_file=vcf_file,
        genome_build=genome_build,
        variant_interpretation_docker=variant_interpretation_docker,
        runtime_attr_override=runtime_attr_bcftools
    }

    scatter (arr in read_tsv(annot_beds_with_header_tsv)) {
        String annot_name = arr[0]
        File ref_bed_with_header = arr[1]
        String match_svtype = arr[2]
    
        call intersectBed {
            input:
            bed_file=vcfToBed.bed_output,
            ref_bed_with_header=ref_bed_with_header,
            match_svtype=match_svtype,
            cohort_prefix=cohort_prefix,
            bed_overlap_threshold=bed_overlap_threshold,
            variant_interpretation_docker=variant_interpretation_docker
        }

        call annotateVCFWithBed {
            input:
            vcf_file=removeGenotypes.no_gt_vcf_file,
            intersect_bed=intersectBed.intersect_bed,
            ref_bed_with_header=ref_bed_with_header,
            genome_build=genome_build,
            hail_docker=hail_docker,
            annot_name=annot_name,
            annotate_sv_from_intersect_bed_script=annotate_sv_from_intersect_bed_script,
            runtime_attr_override=runtime_attr_annotate
        }    
    }

    call combineBedAnnotations {
        input:
        preannotated_vcf=vcf_file,
        preannotated_vcf_idx=vcf_idx,
        annotated_vcfs=annotateVCFWithBed.annotated_vcf,
        annotated_vcfs_idx=annotateVCFWithBed.annotated_vcf_idx,
        genome_build=genome_build,
        variant_interpretation_docker=variant_interpretation_docker,
        runtime_attr_override=runtime_attr_bcftools
    }

    call annotateGeneLevelVCF {
        input:
        vcf_file=combineBedAnnotations.combined_vcf,
        gene_list_tsv=gene_list_tsv,
        omim_uri=omim_uri,
        constrained_uri=constrained_uri,
        prec_uri=prec_uri,
        hi_uri=hi_uri,
        ts_uri=ts_uri,
        permissive_csq_fields=permissive_csq_fields,
        restrictive_csq_fields=restrictive_csq_fields,
        size_threshold=size_threshold,
        dom_af_threshold=dom_af_threshold,
        rec_af_threshold=rec_af_threshold,
        gnomad_af_dom_threshold=gnomad_af_dom_threshold,
        gnomad_af_rec_threshold=gnomad_af_rec_threshold,
        gnomad_popmax_af_threshold=gnomad_popmax_af_threshold,
        gnomad_af_field=gnomad_af_field,
        genome_build=genome_build,
        hail_docker=hail_docker,
        annotate_sv_gene_level_script=annotate_sv_gene_level_script,
        runtime_attr_override=runtime_attr_annotate
    }

    call filterVCF {
        input:
        vcf_file=annotateGeneLevelVCF.annotated_vcf,
        ped_uri=ped_uri,
        genome_build=genome_build,
        hail_docker=hail_docker,
        filter_clinical_sv_script=filter_clinical_sv_script,
        runtime_attr_override=runtime_attr_filter_vcf
    }

    output {
        File sv_pathogenic_tsv = filterVCF.sv_pathogenic_tsv
        File sv_genomic_disorders_tsv = filterVCF.sv_genomic_disorders_tsv
        File sv_large_regions_tsv = filterVCF.sv_large_regions_tsv
        File sv_dominant_tsv = filterVCF.sv_dominant_tsv
        File sv_recessive_tsv = filterVCF.sv_recessive_tsv
        File sv_merged_clinical_tsv = filterVCF.sv_merged_clinical_tsv
        File sv_flagged_vcf = annotateGeneLevelVCF.annotated_vcf
        File sv_flagged_vcf_idx = annotateGeneLevelVCF.annotated_vcf_idx
    }
}

task vcfToBed {
    input{
        File vcf_file
        String cohort_prefix
        String variant_interpretation_docker
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
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -eou pipefail
        svtk vcf2bed --no-samples ~{vcf_file} ~{cohort_prefix}.bed.gz
    >>>

    output {
        File bed_output = "~{cohort_prefix}.bed.gz"
    }
}

task intersectBed {
    input {
        File bed_file
        File ref_bed_with_header
        String match_svtype
        String cohort_prefix
        String variant_interpretation_docker

        Float bed_overlap_threshold
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([bed_file, ref_bed_with_header], 'GB')
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
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String ref_bed_with_header_str = basename(ref_bed_with_header, '.bed')

    command <<<
        set -eou pipefail
        tail -n +2 ~{ref_bed_with_header} > ref.bed 
        if [ "~{match_svtype}" = "True" ]; then
            # Assume 5th column of bed_file is SVTYPE and last column of ref.bed is SVTYPE
            bedtools intersect -wao -f ~{bed_overlap_threshold} -r -a ~{bed_file} -b ref.bed | awk '{if ($5 == $(NF-1)) print}' | bgzip > ~{cohort_prefix}_~{ref_bed_with_header_str}.bed.gz
        else 
            bedtools intersect -wao -f ~{bed_overlap_threshold} -r -a ~{bed_file} -b ref.bed | bgzip > ~{cohort_prefix}_~{ref_bed_with_header_str}.bed.gz
        fi
    >>>

    output {
        File intersect_bed = "~{cohort_prefix}_~{ref_bed_with_header_str}.bed.gz"
    }
}

task annotateVCFWithBed {
    input {
        File vcf_file
        File intersect_bed
        File ref_bed_with_header
        String annot_name
        String genome_build
        String hail_docker
        String annotate_sv_from_intersect_bed_script
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_file, intersect_bed], 'GB')
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
    curl ~{annotate_sv_from_intersect_bed_script} > annotate_vcf.py
    python3 annotate_vcf.py ~{vcf_file} ~{intersect_bed} ~{ref_bed_with_header} ~{genome_build} \
    ~{annot_name} ~{cpu_cores} ~{memory}
    >>>

    output {
        File annotated_vcf = basename(intersect_bed, '.bed.gz') + '.vcf.bgz'
        File annotated_vcf_idx = basename(intersect_bed, '.bed.gz') + '.vcf.bgz.tbi'
    }
}

task combineBedAnnotations {
    input {
        File preannotated_vcf
        File preannotated_vcf_idx
        Array[File] annotated_vcfs
        Array[File] annotated_vcfs_idx
        String genome_build
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(preannotated_vcf, 'GB') + size(annotated_vcfs, 'GB')
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
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String file_ext = if sub(basename(preannotated_vcf), '.vcf.gz', '')!=basename(preannotated_vcf) then '.vcf.gz' else '.vcf.bgz'
    String merged_no_gt_vcf = basename(preannotated_vcf, file_ext) + '.no.GTs.combined.annotations.vcf.gz'
    
    command <<<
    set -eou pipefail
    VCFS="~{write_lines(annotated_vcfs)}"
    cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
    # merge annotations
    bcftools merge --no-version -Oz -o ~{merged_no_gt_vcf} --file-list vcfs_sorted.list
    tabix ~{merged_no_gt_vcf}
    # add genotypes back
    bcftools merge --no-version -Oz -o ~{basename(preannotated_vcf, file_ext) + '.combined.annotations.vcf.gz'} \
        ~{merged_no_gt_vcf} ~{preannotated_vcf}
    tabix ~{basename(preannotated_vcf, file_ext) + '.combined.annotations.vcf.gz'}
    >>>

    output {
        File combined_vcf = basename(preannotated_vcf, file_ext) + '.combined.annotations.vcf.gz'
        File combined_vcf_idx = basename(preannotated_vcf, file_ext) + '.combined.annotations.vcf.gz.tbi'
    }    
}

task annotateGeneLevelVCF {
    input {
        File vcf_file
        File gene_list_tsv
        File omim_uri
        File constrained_uri
        File prec_uri
        File hi_uri
        File ts_uri

        Array[String] permissive_csq_fields
        Array[String] restrictive_csq_fields
        
        Int size_threshold
        Float dom_af_threshold
        Float rec_af_threshold
        Float gnomad_af_dom_threshold
        Float gnomad_af_rec_threshold
        Float gnomad_popmax_af_threshold

        String gnomad_af_field
        String genome_build
        String hail_docker
        String annotate_sv_gene_level_script
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([vcf_file], 'GB')
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
    String output_filename = basename(vcf_file, file_ext) + '.annot.genes.vcf.bgz'

    command <<<
    set -eou pipefail
    curl ~{annotate_sv_gene_level_script} > annotate_vcf.py
    python3 annotate_vcf.py -i ~{vcf_file} -o ~{output_filename} -l ~{gene_list_tsv} -s ~{size_threshold} \
        --omim ~{omim_uri} --cores ~{cpu_cores} --mem ~{memory} --build ~{genome_build} \
        --permissive-csq-fields ~{sep=',' permissive_csq_fields} --restrictive-csq-fields ~{sep=',' restrictive_csq_fields} \
        --constrained-uri ~{constrained_uri} --prec-uri ~{prec_uri} --hi-uri ~{hi_uri} --ts-uri ~{ts_uri} \
        --dom-af ~{dom_af_threshold} --rec-af ~{rec_af_threshold} \
        --gnomad-dom-af ~{gnomad_af_dom_threshold} --gnomad-rec-af ~{gnomad_af_rec_threshold} \
        --gnomad-af-field ~{gnomad_af_field} --gnomad-popmax-af ~{gnomad_popmax_af_threshold}
    >>>

    output {
        File annotated_vcf = output_filename
        File annotated_vcf_idx = output_filename + '.tbi'
    }
}

task removeGenotypes { 
        input {
        File vcf_file
        String genome_build
        String variant_interpretation_docker
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
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String no_gt_vcf = basename(vcf_file, file_ext) + '.no.GTs.vcf.gz'

    command {
        set -eou pipefail
        bcftools view -G -Oz -o ~{no_gt_vcf} --no-version ~{vcf_file}
        tabix ~{no_gt_vcf}
    }

    output {
        File no_gt_vcf_file = no_gt_vcf
        File no_gt_vcf_idx = no_gt_vcf + '.tbi'
    }
}

task renameVCFSamples {
    input {
        File vcf_file
        File sample_map_tsv
        String variant_interpretation_docker
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
        docker: variant_interpretation_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '_renamed_samples.vcf.gz'
    command {
        bcftools reheader -s ~{sample_map_tsv} -o ~{output_filename} ~{vcf_file}
    }

    output {
        File output_vcf = output_filename
    }
}

task filterVCF {
    input {
        File vcf_file
        File ped_uri
        String genome_build
        String hail_docker
        String filter_clinical_sv_script
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

    command <<<
    set -eou pipefail
    curl ~{filter_clinical_sv_script} > filter_vcf.py
    python3 filter_vcf.py -i ~{vcf_file} --ped ~{ped_uri} --cores ~{cpu_cores} --mem ~{memory} --build ~{genome_build}
    >>>

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    output {
        File sv_pathogenic_tsv = basename(vcf_file, file_ext) + '_path_variants.tsv.gz'
        File sv_genomic_disorders_tsv = basename(vcf_file, file_ext) + '_GD_variants.tsv.gz'
        File sv_large_regions_tsv = basename(vcf_file, file_ext) + '_large_regions_variants.tsv.gz'
        File sv_dominant_tsv = basename(vcf_file, file_ext) + '_dominant_variants.tsv.gz'
        File sv_recessive_tsv = basename(vcf_file, file_ext) + '_recessive_variants.tsv.gz'
        File sv_merged_clinical_tsv = basename(vcf_file, file_ext) + '_merged_variants.tsv.gz'
    }
}