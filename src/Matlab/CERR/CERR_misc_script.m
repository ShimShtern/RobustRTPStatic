%% The following code can be used to convert DICOM / RTOG files in a directory (including any files in subdirectories) to CERR format, and place the output in a specified folder
%% Define data source and destination paths 
    
% Structures are stored in the planC cell array and can be accessed as follows:
% Image Registration via command
% line:Wrapper around Plastimatch
% [basePlanC, movPlanC] =
% register_scans(basePlanC, movPlanC,
% baseScanNum, movScanNum, algorithm)

%planC = openPlanCFromFile('CTAC_3_07042021_new_cerr.mat');
global planC
planC = openPlanCFromFile('FusedCTFMISO.mat');

planC = openPlanCFromFile('');

%planC = loadPlanC('FusedCTFMISO.mat');
planC = updatePlanFields(planC);
%planC = openPlanCFromFile('fusedImagesWithStruct07042021.mat');

indexS = planC{end};
structureNum = 1;

%    before using the following code. For further information on planC, see structure of planC
%Influence matrix (The A matrix in eq (1). The size of this matrix is the number of voxels by the number of pencil beams.) For a given IMNumber, potentially defined as follows,
IMNumber = size(planC{indexS.IM},2); %to use the last IM structure created and for structureNumb which corresponds to the desired structure number in planC, the influence matrix can be found in this way:
%dosim = planC{indexS.IM}(IMNumber);
IMNumber
planC{indexS.IM}(IMNumber).IMDosimetry
influenceM = getSingleGlobalInfluenceM(planC{indexS.IM}(IMNumber).IMDosimetry, structureNum);
%influenceM
suvType = 'BW';    
planC = calc_suv(1,planC,suvType); %PET is scan 1 in this file

%%% display a particular slice in PET
slice = 20
figure;
imagesc(planC{indexS.scan}(1).scanArray(:,:,slice))
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% get XYZ coordinates of voxel centers in PET and CT
[PT.x.vals, PT.y.vals, PT.z.vals] = getScanXYZVals(planC{indexS.scan}(1));
[CT.x.vals, CT.y.vals, CT.z.vals] = getScanXYZVals(planC{indexS.scan}(2));



