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
%      seg: a structure with fields (wm,gm,csf,skull,scalp)
%           or a 4D array for which the first tissue volume is the most
%           outter one. Currently assumed scenarios are as follow: 
%           size(seg,4) == 5 assumes 1-Scalp, 2-Skull, 3-CSF, 4. GM, 5. WM
%           size(seg,4) == 4 assumes 1-Scalp, 2-CSF, 3-GM., 4-WM
%           size(seg,4) == 3 assumes 1-CSF, 2-GM, 3-WM
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
%      elem: output, element list of the tetrahedral mesh / the last column denotes the boundary ID
%      face: output, mesh surface element list of the tetrahedral mesh 
%      
% Tissue ID for the outputs are as follow:
% 0-Air/background, 1-Scalp, 2-Skull, 3-CSF, 4-GM, 5-WM
% 
% Methodology behind this function is presented in:
% Anh Phong Tran and Qianqian Fang, "Fast and high-quality tetrahedral mesh generation \
% from neuroanatomical scans,". In: arXiv pre-print (August 2017). arXiv:1708.08954v1 [physics.med-ph]
%


%% Handling the inputs
if nargin >= 2
    options = []; maxvol = []; ratio = 1.414; wh = 0;
elseif nargin >=3
    maxvol = []; ratio = 1.414; options = [];
elseif nargin >= 4
    ratio = 1.414; maxvol = [];
elseif nargin >= 5
    ratio = 1.414;
else
    error('At least two inputs are required');
end

if isempty(res)
    res = 1;
end

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

if isstruct(seg)
    seg2 = seg;
elseif size(seg,4) == 5
    seg2.scalp = seg(:,:,:,1);
    seg2.skull = seg(:,:,:,2);
    seg2.csf = seg(:,:,:,3);
    seg2.gm = seg(:,:,:,4);
    seg2.wm = seg(:,:,:,5); 
elseif size(seg,4) == 4
    seg2.scalp = seg(:,:,:,1);
    seg2.csf = seg(:,:,:,2);
    seg2.gm = seg(:,:,:,3);
    seg2.wm = seg(:,:,:,4);
elseif size(seg,4) == 3
    seg2.csf = seg(:,:,:,1);
    seg2.gm = seg(:,:,:,2);
    seg2.wm = seg(:,:,:,3);
else
    fprintf('This seg input is currently not supported \n')
end

%% Pre-processing steps to create separations between the tissues in the
%% volume space
dim = size(seg2.wm);
seg2.wm = imfill(seg2.wm,'holes');
p_wm = seg2.wm;
p_pial = p_wm+seg2.gm;
p_pial = max(p_pial,max_filter(p_wm,3,1));
p_pial = imfill(p_pial,'holes');
p_csf = p_pial+seg2.csf;
p_csf(p_csf>1) = 1;
p_csf = max(p_csf,max_filter(p_pial,3,1));
if isfield(seg2,'skull') && isfield(seg2,'scalp')
    p_bone = p_csf + seg2.skull;
    p_bone(p_bone>1) = 1;
    p_bone = max(p_bone,max_filter(p_csf,3,1));
    p_skin = p_bone + seg2.scalp;
    p_skin(p_skin>1) = 1;
    p_skin = max(p_skin,max_filter(p_bone,3,1));
    expandedGM = p_pial - seg2.wm - seg2.gm;
	expandedCSF = p_csf - seg2.wm - seg2.gm - seg2.csf - expandedGM;
	expandedSkull = p_bone - seg2.wm - seg2.gm - seg2.csf - seg2.skull - expandedCSF - expandedGM;
	expandedGM = (expandedGM+max_filter(expandedGM,3,1))/2;
	expandedCSF = (expandedCSF+max_filter(expandedCSF,3,1))/2;
	expandedSkull = (expandedSkull+max_filter(expandedSkull,3,1))/2;
elseif isfield(seg2,'scalp') && ~isfield(seg2,'skull')
    p_skin = p_csf + seg2.scalp;
    p_skin(p_skin>1) = 1;
    p_skin = max(p_skin,max_filter(p_csf,3,1));
elseif isfield(seg2,'skull') && ~isfield(seg2,'scalp')
    p_bone = p_csf + seg2.skull;
    p_bone(p_bone>1) = 1;
    p_bone = max(p_bone,max_filter(p_csf,3,1));
end

%% Grayscale/Binary extractions of the surface meshes for the different
%% tissues
[wm_n,wm_f] = v2s(p_wm,THRESH,opt(1),'cgalsurf');
[pial_n,pial_f] = v2s(p_pial,THRESH,opt(2),'cgalsurf');
[csf_n,csf_f] = v2s(p_csf,THRESH,opt(3),'cgalsurf');
if isfield(seg2,'skull')
    [bone_n,bone_f] = v2s(p_bone,THRESH,opt(4),'cgalsurf');
    ISO2MESH_TETGENOPT = '-A';
    [bone_n,el_bone] = surf2mesh(bone_n,bone_f,[],[],1.0,30,[],[],0,'tetgen1.5');
    [bone_f] = volface(el_bone(:,1:4));
    bone_f=removedupelem(bone_f);
    [bone_n,bone_f]=removeisolatednode(bone_n,bone_f);
end
if isfield(seg2,'scalp')
    [skin_n,skin_f] = v2s(p_skin,THRESH,opt(5),'cgalsurf');
end


%% Main loop for the meshing pipeline to combine the individual surface
%% meshes or each of the tissues and to generate the detailed 3D tetrahedral
%% mesh of the brain/head
for w = 1:2
    %% If the first pass fails, a second pass is called using the decoupled function
    %% to eliminate intersections between surface meshes
    if (w==2) && (exist('label_elem','var'))
        continue;
    end
    if (w==2) && (~exist('label_elem','var'))
        [bone_n,bone_f] = surfboolean(bone_n(:,1:3),bone_f(:,1:3),'decouple',skin_n(:,1:3),skin_f(:,1:3));
        [csf_n,csf_f] = surfboolean(csf_n(:,1:3),csf_f(:,1:3),'decouple',bone_n(:,1:3),bone_f(:,1:3));
        [pial_n,pial_f] = surfboolean(pial_n(:,1:3),pial_f(:,1:3),'decouple',csf_n(:,1:3),csf_f(:,1:3));
        [wm_n,wm_f] = surfboolean(wm_n(:,1:3),wm_f(:,1:3),'decouple',pial_n(:,1:3),pial_f(:,1:3));
    end
    [final_n,final_f] = surfboolean(wm_n(:,1:3),wm_f(:,1:3),'resolve',pial_n,pial_f);
    [final_n,final_f] = surfboolean(final_n,final_f,'resolve',csf_n,csf_f);
    if isfield(seg2,'skull')
        [final_n,final_f] = surfboolean(final_n,final_f,'resolve',bone_n,bone_f);
    end
    if isfield(seg2,'scalp')
        [final_n,final_f] = surfboolean(final_n,final_f,'resolve',skin_n,skin_f);
    end
    
    %% If the whole head option is deactivated, the cut is made at the base of the brain using a box cutting
    if (wh == 0)
        dim2 = min(csf_n);
        ISO2MESH_TETGENOPT = '-A';
        [no,f]=meshabox([-1 -1 dim2(3)+4.1],[dim(1)+1 dim(2)+1 dim(3)+1],500);
        [no,f]=removeisolatednode(no,f);
        [final_n,final_f] = surfboolean(no,f(:,[1 3 2]),'first',final_n,final_f);
    end
    
    %% Generates a coarse tetrahedral mesh of the combined tissues
    ISO2MESH_TETGENOPT = '-A ';
    [final_n,final_e,final_f] = surf2mesh(final_n,final_f,[],[],1.0,30,[],[],0,'tetgen1.5');
    
    %% Removes the elements that are part of the box, but not the brain/head
    if (wh == 0)
        [max_node, M] = max(final_n);
        k = find(final_e(:,1:4)==M(3),1);
        final_e = final_e(final_e(:,5)~=final_e(rem(k,length(final_e(:,1))),5),:);
        [final_n,final_e]=removeisolatednode(final_n,final_e);
    end
    
    %% Here the labels created through the coarse mesh generated through Tetgen are saved
    %% with the centroid of one of the elements for intriangulation seg(:,:,:,1)testing later
    label = unique(final_e(:,5));
    for i = 1:length(label)
        label_elem(i,1) = find(final_e(:,5) == label(i),1);
        label_centroid(i,:) = meshcentroid(final_n,final_e(label_elem(i),1:4));
    end
    
    %% This step separates the scalp from the air
    if isfield(seg2,'scalp')
        ISO2MESH_TETGENOPT = '-A';
        [no_skin,el_skin] = surf2mesh(skin_n,skin_f,[],[],1.0,30,[],[],0,'tetgen1.5');
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
    
    %% When the label_elem does not exist, it often indicates a failure at the generation of a coarse
    %% tetrahedral mesh. The alternative meshing pathway using decoupling is then called to make a
    %% second attempt at creating the combined tetrahedral mesh.
    if (~exist('label_elem','var'))&& (w==1)
        fprintf('Initial meshing procedure failed. The option parameter might need to be adjusted. \n')
        fprintf('Activating alternative meshing pathway... \n')
        pause(3)
        continue;
    end
    
    %% The labels are given to each of the tissues
    %% WM(1) - GM(2) - CSF(3) - Bone(4) - Scalp(5) - Air(6)
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
    
    %% This step consolidates adjacent labels of the same tissue
    new_label = unique(final_e(:,5));
    face = [];
    for i = 1:length(new_label)
        face = [face; volface(final_e(final_e(:,5)==new_label(i),1:4))];
    end
    face = sort(face,2);
    face = unique(face,'rows');
    [node,face] = removeisolatednode(final_n,face);
    
    %% The final mesh is generated here with the desired properties
    ISO2MESH_TETGENOPT = sprintf('-A -Rmpq%fa%i',ratio,maxvol);
    node(:,4) = 5*ones(length(node(:,1)),1)/res;
    [node2,elem2,face2] = surf2mesh(node,face,[],[],1.0,30,[],[],0,'tetgen1.5');
    
    label2 = unique(elem2(:,5));
    for i = 1:length(label2)
        label_elem2(i,1) = find(elem2(:,5) == label(i),1);
        label_centroid2(i,:) = meshcentroid(node2,elem2(label_elem2(i),1:4));
    end
    
    %% The labeling process is repeated for the final mesh
    %% WM(1) - GM(2) - CSF(3) - Bone(4) - Scalp(5) - Air(6)
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
    
    if isfield(seg2,'skull') && isfield(seg2,'scalp')
        centroid = ceil(meshcentroid(node2(:,1:3),elem2(:,1:4)));
        for i = 1:length(elem2(:,1))
            if (expandedGM(centroid(i,1),centroid(i,2),centroid(i,3))>0.5) && (elem2(i,5) == 2)
                elem2(i,5) = 3;
            elseif (expandedCSF(centroid(i,1),centroid(i,2),centroid(i,3))>0.5) && (elem2(i,5) == 3)
                elem2(i,5) = 4;
            elseif (expandedSkull(centroid(i,1),centroid(i,2),centroid(i,3))>0.5) && (elem2(i,5) == 4)
                elem2(i,5) = 5;
            end
        end
    end
	
    elem2(:,5) = 6 - elem2(:,5);
end
end