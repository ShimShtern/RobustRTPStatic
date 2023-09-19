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
%reference value from "Defining normoxia, physoxia and hypoxia in
%tumours—implications for treatment response" MCKEOWN 2014
refPt = prctile(suvMnoPTV,50);%.1); % most oxygenated voxel %was 5% for Patient4
ref_pO2=26; 
%refPt = mean(suvMnoPTV)
% conversion P02=(A-Uptake)*C/(Uptake-A+B)
%specific for FMISO Dasu (2012) From Lindblom et al. (2018)
A = 10.9;
B = 10.7; 
C = 2.5;
%A=1.5; %[Chakoyan et al 2017]
%B=0.7;
%C=9.9;

%ref_pO2=60; % Uhurdin et al. (2017)
suv_for_ref_pO2=A-(ref_pO2*B)/(ref_pO2+C);
%suv_for_ref_pO2=1; %changed 21/06/2023
suvMv = suvMv/refPt*suv_for_ref_pO2;
min(suvMv,[],'all')
max(suvMv,[],'all')

pO2=zeros(size(Mask));
pO2(Mask) = (A-suvMv)*C./(suvMv-A+B);
%pO2(suvM==inf)=0;

min(pO2(Mask),[],'all')
max(pO2(Mask),[],'all')
min(pO2(Mask_PTV),[],'all')
max(pO2(Mask_PTV),[],'all')

oer = (m*pO2+K)./(pO2+K);
maxOer = max(oer(MaskNoPTV),[],'all');
omf = zeros(size(oer));
omf(Mask) = oer(Mask)/maxOer;
omf(MaskNoPTV)=0;
deadref=prctile(suvMnoPTV,5)
Mask_PTV_deadcells=(suvM<=deadref);
sum(Mask_PTV_deadcells(:).*Mask_PTV(:))
omf(logical(Mask_PTV_deadcells.*Mask_PTV))=1;

end

