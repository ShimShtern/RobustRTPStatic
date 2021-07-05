%max distance from Neighbors
function Max_Delta=calc_max_omf_diff(Mask_PTV_not_dead,CTOmf)

PTV_inds=find(Mask_PTV_not_dead);
Max_Delta=zeros(size(CTOmf));
delta_space=-1:1:1;
[X,Y,Z] = meshgrid(delta_space,delta_space,delta_space);
for i=PTV_inds'
    [x,y,z]=ind2sub(size(CTOmf),i);    
    xs=x+reshape(X,length(X(:)),1);
    ys=y+reshape(Y,length(Y(:)),1);
    zs=z+reshape(Z,length(Z(:)),1);
    neigbor_inds=sub2ind(size(CTOmf),xs,ys,zs);
    in_ptv=find(Mask_PTV_not_dead(neigbor_inds));
    if ~isempty(in_ptv)
        inds=neigbor_inds(in_ptv);
        max_delta_omf=max(abs(CTOmf(i)-CTOmf(inds)));
        Max_Delta(inds)=max_delta_omf;
    end
end
