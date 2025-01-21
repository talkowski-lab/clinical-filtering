# clinical-filtering
Repository for clinical filtering for SNV/Indels and/or SVs.

Each ```v0.x``` is associated with a published Dockstore workflow. 

**```filterClinicalVariantsSNVIndel_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/clinical-filtering/filter-clinical-variants-SNV-Indel-v01
- SNV/Indel clinical filtering.
  
**```filterClinicalVariantsSV_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/clinical-filtering/filter-clinical-variants-SV-v01
- SV clinical filtering.

**```filterClinicalCompHets_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/clinical-filtering/filter-clinical-comphets-xlr-hom-var-v01
- Outputs comphets, XLR, and homvar in proband from SNV/Indel and/or SV input.

**```filterClinicalVariants_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/clinical-filtering/filter-clinical-variants-v01
- Combines ```filterClinicalVariantsSNVIndel_v0.1.wdl```,```filterClinicalVariantsSV_v0.1.wdl```, and ```filterClinicalCompHets_v0.1.wdl``` into one workflow.

**```filterClinicalVariantsNIFS_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/clinical-filtering/filter-clinical-variants-NIFS-v01
- SNV/Indel only. Run at *sample-set* level.
- Combines ```filterClinicalVariants_v0.1.wdl``` and ```filterClinicalCompHets_v0.1.wdl``` into one workflow, with some differences in filters.

**```filterClinicalVariantsSingleNIFS_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/clinical-filtering/filter-clinical-variants-single-NIFS-v01
- SNV/Indel only. Run at *sample* level.

**```mergeClinicalVariantsNIFS_v0.1.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/clinical-filtering/merge-clinical-variants-NIFS-v01
- Merge runs of ```filterClinicalVariantsSingleNIFS_v0.1.wdl``` (so that when new samples need to be run, not all samples have to be rerun).

**TODO**: Generally want to add scripts to the Docker images for each task/workflow instead of having them as WDL inputs, but not a priority.
