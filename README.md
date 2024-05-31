# Robust Radiotherapy planning

This project involves optimizing for static radiotherapy plans (dose coloring) in the presemce of biological uncertainty.

he code in this repository requires the download of the data files from https://zenodo.org/records/11372979 to be placed in the subfolder
./src/Julia/Data 
where the files should be placed in subdirectories as follows

Liver subdirectory:
liverEx_2.mat

Patient1 subdirectory:

Patient1_Visit1_16beams_refpointpercent50_value26_notincludingdeadvoxels_20230905.mat
Patient1_Visit1_Merged_ptv_brain_withExt_withTelemetry_16beams.mat 

Patient4 subdirectory:

Patient4_Visit1_16beams_refpointpercent50_value26_notincludingdeadvoxels_20230719.mat
Patient4_Visit1_MergedPTCT_structs_20210428_withExt_0102_withDosimetry_16bm.mat
Patient4_Visit1_MergedPTCT_structs_20210428_withExt_0102_withDosimetry_16bm_withDoses_20231115.mat