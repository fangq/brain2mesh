function [node2,elem2,face2] = brain2mesh(seg,res,wh,options,maxvol,ratio)
% Surface-based brain meshing using volumetric data
%
% [node,elem,face] = brain2mesh(seg,res)
%  or
% [node,elem,face] = brain2mesh(seg,res,options)
% [node,elem,face] = brain2mesh(seg,res,options,maxvol)
% [node,elem,face] = brain2mesh(seg,res,options,maxvol,ratio)
%
% input:
%      seg: a 4-D containing the volumetric segmented information
%      res: default is 1, indicates the length of a voxel in mm.
%             e.g. if MRI has a 0.5x0.5x0.5 mm^3 resolution, use 0.5%
%      options: (optional) radius of the Delaunay sphere (element size)
%             Default is [1.7 1.7 2.0 2.5 3] for WM, GM, CSF, Bone and Scalp, respectively.
%             If values are too large, mesh generation may fail.
%      maxvol: (optional) default is 40, indicates the volumetric maximum size of elements
%      ratio: (optional) radius-edge ratio, default is 1.414. Lower values increase 
%             quality better, but results in denser meshes
%
% output:
%      node: output, node coordinates of the tetrahedral mesh
%      elem: output, element list of the tetrahedral mesh
%      face: output, mesh surface element list of the tetrahedral mesh 
%             the last column denotes the boundary ID


if nargin >= 2
    options = []; maxvol = []; ratio = 1.414; wh = 0;
elseif nargin >=3
    maxvol = []; ratio = 1.414; options = [];
elseif nargin >= 4
    ratio = 1.414; maxvol = [];
elseif narging >= 5
    ratio = 1.414;
else
    error('At least two inputs are required');
end

if isempty(res)
    res = 1;
end

dim_seg = size(seg,4);

default_opt = [1.7 1.7 2 2.5 3];
THRESH = 0.5;
MAX_NODE = 40000/res;
for i = 1:5
    opt(i).maxnode = MAX_NODE;
end
if isempty(options)
    for k = 1:5
        opt(k).radbound = default_opt(k)/res;
    end
elseif length(options) == 4
    for k = 1:3
        opt(k).radbound = options(k)/res;
    end
    opt(4).radbound = 2.5;
    opt(5).radbound = options(4)/res;
else
    for k = 1:length(options)
        opt(k).radbound = options(k)/res;
    end
end

dim = size(seg);
seg(:,:,:,1) = imfill(seg(:,:,:,1),'holes');
p_wm = seg(:,:,:,1);
p_pial = p_wm+seg(:,:,:,2);
p_pial = max(p_pial,max_filter(p_wm,3,1));
p_pial = imfill(p_pial,'holes');
p_csf = p_pial+seg(:,:,:,3);
p_csf(p_csf>1) = 1;
p_csf = max(p_csf,max_filter(p_pial,3,1));

if dim_seg > 4
    p_bone = p_csf + seg(:,:,:,4);
    p_bone(p_bone>1) = 1;
    p_bone = max(p_bone,max_filter(p_csf,3,1));
end
if dim_seg == 4
    p_skin = p_csf + seg(:,:,:,4);
    p_skin(p_skin>1) = 1;
    p_skin = max(p_skin,max_filter(p_csf,3,1));
elseif dim_seg >3
    p_skin = p_bone + seg(:,:,:,5);
    p_skin(p_skin>1) = 1;
    p_skin = max(p_skin,max_filter(p_bone,3,1));
end

[wm_n,wm_f] = v2s(p_wm,THRESH,opt(1),'cgalsurf');
[pial_n,pial_f] = v2s(p_pial,THRESH,opt(2),'cgalsurf');
[csf_n,csf_f] = v2s(p_csf,THRESH,opt(3),'cgalsurf');
if dim_seg > 4
    [bone_n,bone_f] = v2s(p_bone,THRESH,opt(4),'cgalsurf');
end
if dim_seg > 3
    [skin_n,skin_f] = v2s(p_skin,THRESH,opt(5),'cgalsurf');
end

if dim_seg > 4
    ISO2MESH_TETGENOPT = '-A';
    [bone_n,el_bone] = surf2mesh(bone_n,bone_f,[],[],1.0,30,[],[]);
    [bone_f] = volface(el_bone(:,1:4));
    bone_f=removedupelem(bone_f);
    [bone_n,bone_f]=removeisolatednode(bone_n,bone_f);
end

for w = 1:2
    if (w==2) && (exist('label_elem'))
        continue;
    end
    if (w==2) && (~exist('label_elem'))
        [bone_n,bone_f] = surfboolean(bone_n(:,1:3),bone_f(:,1:3),'decouple',skin_n(:,1:3),skin_f(:,1:3));
        [csf_n,csf_f] = surfboolean(csf_n(:,1:3),csf_f(:,1:3),'decouple',bone_n(:,1:3),bone_f(:,1:3));
        [pial_n,pial_f] = surfboolean(pial_n(:,1:3),pial_f(:,1:3),'decouple',csf_n(:,1:3),csf_f(:,1:3));
        [wm_n,wm_f] = surfboolean(wm_n(:,1:3),wm_f(:,1:3),'decouple',pial_n(:,1:3),pial_f(:,1:3));
    end
    [final_n,final_f] = surfboolean(wm_n(:,1:3),wm_f(:,1:3),'resolve',pial_n,pial_f);
    [final_n,final_f] = surfboolean(final_n,final_f,'resolve',csf_n,csf_f);
    if dim_seg > 4
        [final_n,final_f] = surfboolean(final_n,final_f,'resolve',bone_n,bone_f);
    end
    if dim_seg > 3
        [final_n,final_f] = surfboolean(final_n,final_f,'resolve',skin_n,skin_f);
    end
    
    if (wh == 0)
        dim2 = min(csf_n);
        ISO2MESH_TETGENOPT = '-A';
        [no,f]=meshabox([-1 -1 dim2(3)+4.1],[dim(1)+1 dim(2)+1 dim(3)+1],500);
        [no,f]=removeisolatednode(no,f);
        [final_n,final_f] = surfboolean(no,f(:,[1 3 2]),'first',final_n,final_f);
    end
    
    ISO2MESH_TETGENOPT = '-A ';
    [final_n,final_e,final_f] = surf2mesh(final_n,final_f,[],[],1.0,30,[],[]);
    
    if (wh == 0)
        [max_node, M] = max(final_n);
        k = find(final_e(:,1:4)==M(3),1);
        final_e = final_e(final_e(:,5)~=final_e(rem(k,length(final_e(:,1))),5),:);
        [final_n,final_e]=removeisolatednode(final_n,final_e);
    end
    
    label = unique(final_e(:,5));
    for i = 1:length(label)
        label_elem(i,1) = find(final_e(:,5) == label(i),1);
        label_centroid(i,:) = meshcentroid(final_n,final_e(label_elem(i),1:4));
    end
    
    if dim_seg > 3
        ISO2MESH_TETGENOPT = '-A';
        [no_skin,el_skin] = surf2mesh(skin_n,skin_f,[],[],1.0,30,[],[]);
        for i = 1:length(unique(el_skin(:,5)))
            vol_skin(i) = sum(elemvolume(no_skin,el_skin(el_skin(:,5)==i,1:4)));
        end
        [Y,I] = max(vol_skin);
        if (length(unique(el_skin(:,5)))>1)
            no_air = no_skin; el_air = el_skin(el_skin(:,5)~=I,:);
            [no_air,el_air]=removeisolatednode(no_air,el_air);
        end
        el_skin = el_skin(el_skin(:,5)==I,:);
        [no_skin,el_skin]=removeisolatednode(no_skin,el_skin);
        
        [f_skin] = volface(el_skin(:,1:4));
        f_skin=removedupelem(f_skin);
    end
    
    if (~exist('label_elem'))&& (w==1)
        fprintf('Initial meshing procedure failed. The option parameter might need to be adjusted. \n')
        fprintf('Activating alternative meshing pathway... \n')
        pause(3)
        continue;
    end
    
    for i = 1:length(label_elem(:,1))
        if intriangulation(wm_n,wm_f(:,1:3),label_centroid(i,:))%inside_test(no_wm,el_wm(:,1:4),label_centroid(i,:))
            label_label(i,1) = 1;
        elseif intriangulation(pial_n,pial_f(:,1:3),label_centroid(i,:))%inside_test(no_pial,el_pial(:,1:4),label_centroid(i,:))
            label_label(i,1) = 2;
        elseif intriangulation(csf_n,csf_f(:,1:3),label_centroid(i,:))
            label_label(i,1) = 3;
        elseif (exist('bone_n'))
            if intriangulation(bone_n,bone_f(:,1:3),label_centroid(i,:))
                label_label(i,1) = 4;
            elseif intriangulation(no_skin,f_skin(:,1:3),label_centroid(i,:))
                label_label(i,1) = 5;
            else
                label_label(i,1) = 6;
            end
        elseif (exist('no_skin'))
            if intriangulation(no_skin,f_skin(:,1:3),label_centroid(i,:))
                label_label(i,1) = 5;
            else
                label_label(i,1) = 6;
            end
        else
            label_label(i,1) = 6;
        end
    end
    
    for i = 1:length(final_e(:,1))
        final_e(i,5) = label_label(final_e(i,5));
    end
    
    new_label = unique(final_e(:,5));
    face = [];
    for i = 1:length(new_label)
        face = [face; volface(final_e(final_e(:,5)==new_label(i),1:4))];
    end
    face = sort(face,2);
    face = unique(face,'rows');
    [node,face] = removeisolatednode(final_n,face);
    
    ISO2MESH_TETGENOPT = sprintf('-A -Rmpq%fa%i',ratio,maxvol);
    node(:,4) = 5*ones(length(node(:,1)),1)/res;
    [node2,elem2,face2] = surf2mesh(node,face,[],[],1.0,30,[],[]);
    
    label2 = unique(elem2(:,5));
    for i = 1:length(label2)
        label_elem2(i,1) = find(elem2(:,5) == label(i),1);
        label_centroid2(i,:) = meshcentroid(node2,elem2(label_elem2(i),1:4));
    end
    
    for i = 1:length(label_elem2(:,1))
        if intriangulation(wm_n,wm_f(:,1:3),label_centroid2(i,:))
            label_label2(i,1) = 1;
        elseif intriangulation(pial_n,pial_f(:,1:3),label_centroid2(i,:))
            label_label2(i,1) = 2;
        elseif intriangulation(csf_n,csf_f(:,1:3),label_centroid2(i,:))
            label_label2(i,1) = 3;
        elseif (exist('bone_n'))
            if intriangulation(bone_n,bone_f(:,1:3),label_centroid2(i,:))
                label_label2(i,1) = 4;
            elseif intriangulation(no_skin,f_skin(:,1:3),label_centroid2(i,:))
                label_label2(i,1) = 5;
            else
                label_label2(i,1) = 6;
            end
        elseif (exist('no_skin'))
            if intriangulation(no_skin,f_skin(:,1:3),label_centroid2(i,:))
                label_label2(i,1) = 5;
            else
                label_label2(i,1) = 6;
            end
        else
            label_label2(i,1) = 6;
        end
    end
    
    for i = 1:length(elem2(:,1))
        elem2(i,5) = label_label2(elem2(i,5));
    end
    
end
