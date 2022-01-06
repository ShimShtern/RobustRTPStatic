function [CTOmf,Mask_PTV_not_dead,Max_Delta] = ct_suv_to_omf(CTindex,PTindex,strructNumBrain,structNumPTV)
global planC
%if ~exist('fname','var')
%    fname = 'C:\Users\Shim\Dropbox (MIT)\technion\research\Robust Radiotherapy\data\ACRIN-FMISO-Brain\Segmented data\Patient 4\MergedPTCT_structs_20210428_withExt_0102_withDosimetry.mat'; 
%end
if ~exist('CTindex','var')
    CTindex = 1;
end
if ~exist('PTindex','var')
    PTindex = 2;
end
if ~exist('strructNumBrain','var')
    strructNumBrain = 5;
end
if ~exist('structNumPTV','var')
    structNumPTV = 4;
end

%planC = loadPlanC(fname);
planC = updatePlanFields(planC);
indexS = planC{end};
suvType = 'BW';
planC = calc_suv(PTindex,planC,suvType); 
[CT, PT]= convert_suv_to_CT(planC,PTindex,CTindex);
Mask = getUniformStr(strructNumBrain,planC)>0; %mask of brain;
Mask_PTV = getUniformStr(structNumPTV,planC)>0; %mask of ptv;
[CTOmf, CT_Mask_deadcells] = suv_to_omf(CT.suv,Mask,Mask_PTV);
Mask_PTV_not_dead=logical(Mask_PTV.*(1-CT_Mask_deadcells));
Max_Delta=calc_max_omf_diff(Mask_PTV_not_dead,CTOmf);

%Mask_PTV_PT =getUniformStr(1,planC)>0; %mask of ptv on PET;
%new_suv_PT=planC{indexS.scan}(PTindex).scanArray.*Mask_PTV_PT;
%figure; imagesc(new_suv_PT(:,:,21))
end