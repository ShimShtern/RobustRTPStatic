function [D,V]=CreateIM(Target_S_num, Organ_S_nums)
%%
% Target_S_num -structure num of target
% strOrgan_S_nums -array with structure nums of organs 
% CT_IMNumber - number of SC scan
%%
global planC
indexS = planC{length(planC)};

V = cell(numel(Organ_S_nums)+1,1);
% create target IM
IMNumber=1; %assuming there is only one ifluence matrix created
influenceM = getSingleGlobalInfluenceM(planC{indexS.IM}(IMNumber).IMDosimetry, Target_S_num);
V{1}=1:size(influenceM,1);
V{1}=V{1}(any(influenceM'));
D = sparse(influenceM);
%create organ IM
for i=1:numel(Organ_S_nums)
    influenceM = getSingleGlobalInfluenceM(planC{indexS.IM}(IMNumber).IMDosimetry, Organ_S_nums(i));
    V{i+1}=(1:size(influenceM,1));%+V{i}(end);
    V{i+1}=V{i+1}(any(influenceM'));
    D = max(D,influenceM);
end
D = sparse(D);

