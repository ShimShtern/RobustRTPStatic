
global planC
%fileN = 'C:\Users\Shim\Dropbox (MIT)\technion\research\Robust Radiotherapy\data\CERR\MergedPTCT_structs_20210428.mat';
FileName='Patient1_Visit1_Merged_ptv_brain';
FileDir='C:\Users\Shim\Dropbox (MIT)\technion\research\Robust Radio Therapy (All)\Robust Radiotherapy\data\ACRIN-FMISO-Brain\Segmented data\Patient 1\';
fileFullName = [FileDir,FileName,'.mat'];
planC = loadPlanC(fileFullName);
planC = updatePlanFields(planC);
%planC = quality_assure_planC(movCerrFile,planC);
indexS = planC{end};

scanInfoindex=indexS.scan;
sTypes=convertCharsToStrings({planC{scanInfoindex}(:).scanType});
scanNum=find(startsWith(sTypes,"CT"));
scan3M = getScanArray(scanNum,planC);
strMask3M = getPatientOutline(scan3M, 1:size(scan3M,3),0);

structureName = 'Exterior';
% true - to indicate it is a "uniform CT scan" ?
planC = maskToCERRStructure(strMask3M, true, scanNum, structureName, planC);
planC = updatePlanFields(planC);
fileN = 'C:\Users\Shim\Dropbox (MIT)\technion\research\Robust Radio Therapy (All)\Robust Radiotherapy\data\ACRIN-FMISO-Brain\Segmented data\Patient 1\Patient1_Visit1_Merged_ptv_brain.mat';
FileNameNew=[FileDir,FileName,'_withExt.mat'];
save_planC(planC)