
global planC
fileN = 'C:\Users\Shim\Dropbox (MIT)\technion\research\Robust Radiotherapy\data\CERR\MergedPTCT_structs_20210428.mat';
planC = loadPlanC(fileN);
planC = updatePlanFields(planC);
%planC = quality_assure_planC(movCerrFile,planC);
indexS = planC{end};

scanNum = 2;
scan3M = getScanArray(scanNum,planC);
strMask3M = getPatientOutline(scan3M, 1:size(scan3M,3),0);

structureName = "Exterior";
% true - to indicate it is a "uniform CT scan" ?
planC = maskToCERRStructure(strMask3M, true, scanNum, structureName, planC);
