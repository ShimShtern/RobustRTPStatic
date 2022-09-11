 function [omf, Mask_PTV_deadcells] = suv_to_omf(suvM, Mask,Mask_PTV)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
EPS = 1e-3;

K = 3; %constant in formula of Titz and Jeraj (2008), and choice of value that seems to have been used by
% Saka, Rardin, Langer (2014)
m = 3;

% normalize SUV - convert to "unitless" T/B 
Mask=(Mask>0);
Mask_PTV=(Mask_PTV>0);
MaskNoPTV = ((Mask).*(1-(Mask_PTV)))>0; %nullify PTV voxels
size(suvM)
size(Mask)
suvMv = suvM(Mask);
suvMnoPTV = suvM(MaskNoPTV);
refPt = prctile(suvMnoPTV,5); % most oxygenated voxel %was 5% for Patient4
suvMv = suvMv/refPt;
min(suvMv,[],'all')
max(suvMv,[],'all')
%specific for FMISO
A = 10.9;
B = 10.7; %A-min(suvM,[],'all')-EPS;
C = 2.5;
% conversion P02=(A-Uptake)*C/(Uptake-A+B), where A=10.9 and B=10.7, Lindblom, Dasu,
% Uhurdin et al. (2017)
pO2=zeros(size(Mask));
pO2(Mask) = (A-suvMv)*C./(suvMv-A+B);
%pO2(suvM==inf)=0;

min(pO2(Mask),[],'all')
max(pO2,[],'all')

oer = (m*pO2+K)./(pO2+K);
maxOer = max(oer(Mask),[],'all');
omf = zeros(size(oer));
omf(Mask) = oer(Mask)/maxOer;
omf(MaskNoPTV)=0;
Mask_PTV_deadcells=(suvM<=refPt);
sum(Mask_PTV_deadcells(:).*Mask_PTV(:))
omf(logical(Mask_PTV_deadcells.*Mask_PTV))=1;

end

