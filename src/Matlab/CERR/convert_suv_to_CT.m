function [CT, PT]= convert_suv_to_CT(planC,PTindex,CTindex)
%%calc SUV converted to CT vovels.
%input:
%PlanC - global structure
%PTindex - scan index of PET
%CTindex - scan index of CT
%output
%CT - CT structure including interpulated suv values
%PT - PET structure

indexS = planC{end}; %structure containing fields of PlanC
%getting conversion between pixels and slices to xyz values
[PT.x.vals, PT.y.vals, PT.z.vals] = getScanXYZVals(planC{indexS.scan}(PTindex));
[CT.x.vals, CT.y.vals, CT.z.vals] = getScanXYZVals(planC{indexS.scan}(CTindex));
%was used for normalization, was removed
%max_suv=max(max(max(planC{indexS.scan}(1).scanArray))); %this is supposed to be suv, but not [0,1]
%min_suv=min(min(min(planC{indexS.scan}(1).scanArray)));

%assumes CT overlaps at most 2 PT voxels in each dimension
for dim=['x','y','z']
    PT.(dim).n=length(PT.(dim).vals);%compute number of voxelds along axis
    CT.(dim).n=length(CT.(dim).vals);
    PT.(dim).diff=abs(PT.(dim).vals(2)-PT.(dim).vals(1))/2; %compute width of voxel along axis
    CT.(dim).diff=abs(CT.(dim).vals(2)-CT.(dim).vals(1))/2;
    CT.(dim).Ratio=sparse(zeros(CT.(dim).n,PT.(dim).n)); %gives the ratio of each PET pixel in each CT pixel
    for i=1:CT.(dim).n
        %find two voxels that the CT voxels intersects in this dimension
        upper_voxel=find((CT.(dim).vals(i)+CT.(dim).diff<=PT.(dim).vals+PT.(dim).diff).*(CT.(dim).vals(i)+CT.(dim).diff>=PT.(dim).vals-PT.(dim).diff));
        lower_voxel=find((CT.(dim).vals(i)-CT.(dim).diff<=PT.(dim).vals+PT.(dim).diff).*(CT.(dim).vals(i)-CT.(dim).diff>=PT.(dim).vals-PT.(dim).diff));
        %compute the ratio of each of the PET voxels in the CT voxels
        if upper_voxel==1
            CT.(dim).Ratio(i,1)=1;
        elseif lower_voxel==PT.(dim).n
            CT.(dim).Ratio(i,PT.(dim).n)=1;
        elseif lower_voxel==upper_voxel
            CT.(dim).Ratio(i,lower_voxel)=1;
        elseif ~isempty(upper_voxel)
            CT.(dim).Ratio(i,lower_voxel)=((PT.(dim).vals(lower_voxel)+PT.(dim).diff)-(CT.(dim).vals(i)-CT.(dim).diff))/(2*CT.(dim).diff);
            CT.(dim).Ratio(i,upper_voxel)=((CT.(dim).vals(i)+CT.(dim).diff)-(PT.(dim).vals(upper_voxel)-PT.(dim).diff))/(2*CT.(dim).diff);
        end
        %saves the mapping from the PET to CT voxels
        CT.(dim).Mapping{i}=find(CT.(dim).Ratio(i,:));
    end
end
%create the CT suv array
%initialize
CT.suv=zeros(size(planC{indexS.scan}(CTindex).scanArray));
for CTslice=1:size(planC{indexS.scan}(CTindex).scanArray,3) %go over CT slices
    PTslices=CT.z.Mapping{CTslice};
    %for each CT slice compute the suvs of both PET slices and iterpulate
    %them
    if length(PTslices)>1
        suv1=planC{indexS.scan}(PTindex).scanArray(:,:,PTslices(1));
        suv2=planC{indexS.scan}(PTindex).scanArray(:,:,PTslices(2));
        suv = suv1* full(CT.z.Ratio(CTslice,PTslices(1)))+suv2* full(CT.z.Ratio(CTslice,PTslices(2)));
    elseif length(PTslices)==1    
        suv=planC{indexS.scan}(PTindex).scanArray(:,:,PTslices(1));
    else
        CT.suv(:,:,CTslice)=Inf;
        continue;
    end
    %for each pixel in the slice compute the interpulation based on the x
    %and y PET ratios
    for x_ind=1:CT.x.n
        for y_ind=1:CT.y.n
            x_nonzero=CT.x.Mapping{x_ind};
            y_nonzero=CT.y.Mapping{y_ind};
            if isempty(x_nonzero) || isempty(y_nonzero)
                CT.suv(x_ind,y_ind,CTslice)=inf;
            else    
                CT.suv(x_ind,y_ind,CTslice)=(full(CT.x.Ratio(x_ind,x_nonzero))*suv(x_nonzero,y_nonzero))*full(CT.y.Ratio(y_ind,y_nonzero)');
            end
        end   
    end
end
end