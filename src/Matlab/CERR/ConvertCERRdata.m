global planC

%fname="C:\Users\Shim\Dropbox (MIT)\technion\research\Robust Radiotherapy\data\ACRIN-FMISO-Brain\Segmented data\Patient4\MergedPTCT_structs_20210428_withExt_0102_withDosimetry2.mat";
%dir='C:\Users\Shim\Dropbox (MIT)\technion\research\Robust Radio Therapy (All)\Robust Radiotherapy\data\ACRIN-FMISO-Brain\Segmented data\Patient1\';
%fname=[dir,'Patient1_Visit1_Merged_ptv_brain_withExt_withTelemetry_16beams.mat'];
StructNumPTV=4;
StructNumBrain=3;
StructNumExterior=5;

fname = "C:\Users\Shimrit\Dropbox (MIT)\technion\research\Robust Radio Therapy (All)\Robust Radiotherapy\data\ACRIN-FMISO-Brain\Segmented data\Patient4\Patient4_Visit1_MergedPTCT_structs_20210428_withExt_0102_withDosimetry_16bm.mat";
StructNumPTV=4;
StructNumBrain=5;
StructNumChiasm=6;
StructNumExterior=7;


planC = loadPlanC(fname);
planC = updatePlanFields(planC);
indexS = planC{end};
CTindex=1;

CTscan=1;
PETScan=2;
[omf_Vec,Mask_PTV_not_dead,Max_Delta] = ct_suv_to_omf(CTscan,PETScan,StructNumBrain,StructNumPTV);
omf_Vec=omf_Vec(Mask_PTV_not_dead);
[Dij,V]=CreateIM(StructNumPTV,[StructNumBrain,StructNumChiasm,StructNumExterior]);% for Patient 4
V{4}=setdiff(V{4},V{2}); %Exterior-Brain
V{4}=setdiff(V{4},V{1}); %Exterior-PTV   
V{2}=setdiff(V{2},V{1}); %Brain-PTV   
%[Dij,V]=CreateIM(StructNumPTV, [StructNumBrain,StructNumExterior]); %for Patient 1
Dij=getGlobalInfluenceM(planC{indexS.IM}(CTscan).IMDosimetry,[StructNumPTV,StructNumBrain,StructNumChiasm,StructNumExterior]);


%%
V{1}=find(Mask_PTV_not_dead)'; %update to just include living voxels
T_voxels = V{1}';%V{1} always PTV
T_voxel_num=numel(T_voxels); 
sz_voxels=size(planC{indexS.scan}(CTindex).scanArray);%size of voxel in xyz
[x,y,z]= ind2sub(sz_voxels,T_voxels);
threeD_T_voxels=[x,y,z];

ntriplets = T_voxel_num*27 ; %upper bound on number of neigbors
I = zeros (ntriplets, 1) ;
J = zeros (ntriplets, 1) ;
X = ones (ntriplets, 1) ;
ntriplets=0;
for i=1:T_voxel_num
    i
    for j=i+1:T_voxel_num
        if abs(x(i)-x(j))<=1 && abs(y(i)-y(j))<=1 && abs(z(i)-z(j))<=1
            ntriplets=ntriplets+1;
            I(ntriplets)=i;
            J(ntriplets)=j;
            ntriplets=ntriplets+1;
            J(ntriplets)=i;
            I(ntriplets)=j;
        end
    end
end
I((ntriplets+1):end)=[];
J((ntriplets+1):end)=[];
X((ntriplets+1):end)=[];
neighbors_Mat=sparse(I,J,X,T_voxel_num,T_voxel_num);

%save([dir,'Patient1_Visit1_16beams_refpoint0.1.mat'],'Dij','V', 'neighbors_Mat', 'omf_Vec','-v7.3')
save('Patient4_Visit1_16beams_refpoint5_notincludingdeadvoxels_20220810.mat','Dij','V', 'neighbors_Mat', 'omf_Vec','-v7.3')
