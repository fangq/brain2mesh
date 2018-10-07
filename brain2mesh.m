function [brain_n,brain_el,brain_f] = brain2mesh(seg,cfg)
%
% Brain2mesh: a one-liner for human brain 3D mesh generation
% 
% == Format == 
% [node,elem,face] = brain2mesh(seg,cfg); 
%     or 
% [node,elem,face] = brain2mesh(seg)
% 
% == Input ==
%      seg: pre-segmented brain volume (supporting both probalistic tissue 
%           segmentation and labeled volume). Two formats are accepted
%           1. a structure with subfields (wm,gm,csf,skull,scalp)
%              e.g.: seg.wm, seg.gm, seg.csf represents the white-matter, 
%              gray-matter and csf segmentaions, respectively,
%                  or 
%           2. a 4D array for with tissued sorted in outer-to-inner order
%              the 4th dimension of the array can 3-6, with the following assumptions
%              size(seg,4) == 6 assumes 1-Scalp, 2-Skull, 3-CSF, 4-GM, 5-WM, 6-air pockets
%              size(seg,4) == 5 assumes 1-Scalp, 2-Skull, 3-CSF, 4-GM, 5-WM
%              size(seg,4) == 4 assumes 1-Scalp, 2-CSF, 3-GM, 4-WM
%              size(seg,4) == 3 assumes 1-CSF, 2-GM, 3-WM
%
%      cfg: a struct defining the options for the resulting tetrahedral mesh
%          default values are applied if field is not defined
%          cfg.samplingwm /.samplinggm/.samplingcsf/.samplingskull/.samplingscalp:
%             Radius of the Delaunay sphere used in the sampling the surfaces.
%             Default values are 1.7, 1.7, 2, 2.5 and 3, respectively (reference values for 1x1x1mm^3)
%             Scale proportionally for denser volumes. Lower values correspond to denser, higher
%             fidelity surface extraction, but also results in denser meshes.
%          cfg.maxnode: [120000] - when  the value cfg.sampling__ creates surfaces that are too 
%             dense. This limits the maximum of number of nodes extracted for a given surface.
%          cfg.maxvol: [30] indicates the volumetric maximum size of elements
%             Lowering this value helps with obtaining a denser tetrahedral
%             mesh. For dense meshes, values close to 3-5 are recommended.
%          cfg.sizefield: [3] - mesh sizing field used at the mesh generation step.
%             Higher values can help lowering the final mesh density.         
%          cfg.ratio: [1.414] radius-edge ratio. Lower values increase 
%             the quality of tetrahedral elements, but results in denser meshes
%          cfg.relabeling: [0] or 1 - This step removes most of the assumptions 
%             created by the layered meshing workflow. Currently only works if all five tissue types are present.
%             When deactivated, a 1 voxel length gap is assumed between each of the tissue layers.
%          cfg.skullair: 0 or [1]. Within the skull layer, self-contained entities can be created. 
%             By default, these entities are merged to the skull because they can be ambiguously be
%             vessels or air. When the option 1 is chosen, these entities
%             are labeled as air/background instead.
%          cfg.wh: [0] or 1. by default, the mesh is truncated a few pixel in the 
%             +z direction below the lowest pixel containing CSF to focus on the brain areas. 
%             A value of 1 gives a complete head mesh. 
% 
% == Outputs ==
%      node: node coordinates of the tetrahedral mesh
%      elem: element list of the tetrahedral mesh / the last column denotes the boundary ID
%      face: mesh surface element list of the tetrahedral mesh 
%      
% Tissue ID for the outputs are as follow:
% 0-Air/background, 1-Scalp, 2-Skull, 3-CSF, 4-GM, 5-WM, 6-air pockets
% 
% == Methodology ==
% The underlying methodology behind this function is described in:
% Anh Phong Tran and Qianqian Fang, "Fast and high-quality tetrahedral mesh generation \
% from neuroanatomical scans,". In: arXiv pre-print (August 2017). arXiv:1708.08954v1 [physics.med-ph]
%

%% Handling the inputs
if nargin == 1
elseif nargin ~= 2
    error('Number of input parameters should be either 1 or 2') 
end

if (~exist('cfg'))
    radbd=struct('wm',1.7,'gm',1.7,'csf',2.0,'skull',2.5,'scalp',3.0);
    cfg = struct('ratio',1.414,'wh',0,'maxvol',30,'sizefield',3,'relabeling',0,'radbound',radbd,'skullair',1,'maxnode',120000);
end

if (~isfield(cfg,'ratio'))
    cfg.ratio = 1.414;
end
if (~isfield(cfg,'wh'))
    cfg.wh = 0;
end
if (~isfield(cfg,'maxvol'))
    cfg.maxvol = 30;
end
if (~isfield(cfg,'sizefield'))
    cfg.sizefield = 3;
end
if (~isfield(cfg,'relabeling'))
    cfg.relabeling = 0;
end
if (~isfield(cfg,'samplingwm'))
    cfg.samplingwm = 1.7;
end
if (~isfield(cfg,'samplinggm'))
    cfg.samplinggm = 1.7;
end
if (~isfield(cfg,'samplingcsf'))
    cfg.samplingcsf = 2.0;
end
if (~isfield(cfg,'samplingskull'))
    cfg.samplingskull = 2.5;
end
if (~isfield(cfg,'samplingscalp'))
    cfg.samplingscalp = 3.0;
end
if (~isfield(cfg,'skullair'))
    cfg.skullair = 1;
end
if (~isfield(cfg,'maxnode'))
    cfg.maxnode = 120000;
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

THRESH = 0.5;
for i = 1:5
    opt(i).maxnode = cfg.maxnode; 
end
opt(1).radbound = cfg.samplingwm; opt(2).radbound = cfg.samplinggm; 
opt(3).radbound = cfg.samplingcsf; opt(4).radbound = cfg.samplingskull; 
opt(5).radbound = cfg.samplingscalp; 

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
expandedGM = p_pial - seg2.wm - seg2.gm;
expandedCSF = p_csf - seg2.wm - seg2.gm - seg2.csf - expandedGM;
expandedGM = max_filter(expandedGM,3,1);
expandedCSF = max_filter(expandedCSF,3,1);
if isfield(seg2,'skull') && isfield(seg2,'scalp')
    p_bone = p_csf + seg2.skull;
    p_bone(p_bone>1) = 1;
    p_bone = max(p_bone,max_filter(p_csf,3,1));
    p_skin = p_bone + seg2.scalp;
    p_skin(p_skin>1) = 1;
    p_skin = max(p_skin,max_filter(p_bone,3,1));
	expandedSkull = p_bone - seg2.wm - seg2.gm - seg2.csf - seg2.skull - expandedCSF - expandedGM;
    expandedSkull = max_filter(expandedSkull,3,1);
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
    [bone_node,el_bone] = surf2mesh(bone_n,bone_f,[],[],1.0,30,[],[],0,'tetgen1.5');
    for i = 1:length(unique(el_bone(:,5)))
        vol_bone(i) = sum(elemvolume(bone_node,el_bone(el_bone(:,5)==i,1:4)));
    end
    [~,I] = max(vol_bone);
    if (length(unique(el_bone(:,5)))>1)
        no_air2 = bone_node; el_air2 = el_bone(el_bone(:,5)~=I,:);
        [no_air2,el_air2]=removeisolatednode(no_air2,el_air2);
        f_air2 = volface(el_air2(:,1:4));
    end
    bone_n2 = bone_node;
    [bone_f2] = volface(el_bone(:,1:4));
    bone_f2 = removedupelem(bone_f2);
    [bone_n2,bone_f2]=removeisolatednode(bone_n2,bone_f2);
    if cfg.skullair == 0
        bone_n = bone_n2; bone_f = bone_f2;
    end
end
if isfield(seg2,'scalp')
    [skin_n,skin_f] = v2s(p_skin,THRESH,opt(5),'cgalsurf');
end

%% Main loop for the meshing pipeline to combine the individual surface
%% meshes or each of the tissues and to generate the detailed 3D tetrahedral
%% mesh of the brain/head
for loop = 1:2
    %% If the first pass fails, a second pass is called using the decoupled function
    %% to eliminate intersections between surface meshes
    if (loop==2) && (exist('label_elem'))
        continue;
    end
    if (loop==2) && (~exist('label_elem'))
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
    if (cfg.wh == 0)
        dim2 = min(csf_n);
        ISO2MESH_TETGENOPT = '-A';
        [no,f]=meshabox([-1 -1 dim2(3)+4.1],[dim(1)+1 dim(2)+1 dim(3)+1],500);
        [no,f]=removeisolatednode(no,f);
        [final_n,final_f] = surfboolean(no,f(:,[1 3 2]),'first',final_n,final_f);
    end
    
    %% Generates a coarse tetrahedral mesh of the combined tissues
    ISO2MESH_TETGENOPT = '-A ';
    [final_n,final_e,~] = surf2mesh(final_n,final_f,[],[],1.0,30,[],[],0,'tetgen1.5');

    %% Removes the elements that are part of the box, but not the brain/head
    if (cfg.wh == 0)
        [~, M] = max(final_n);
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

    if isfield(seg2,'scalp')
        ISO2MESH_TETGENOPT = '-A';
        [no_skin,el_skin] = surf2mesh(skin_n,skin_f,[],[],1.0,30,[],[],0,'tetgen1.5');
        for i = 1:length(unique(el_skin(:,5)))
            vol_skin(i) = sum(elemvolume(no_skin,el_skin(el_skin(:,5)==i,1:4)));
        end
        [~,I] = max(vol_skin);
        if (length(unique(el_skin(:,5)))>1)
            no_air = no_skin; el_air = el_skin(el_skin(:,5)~=I,:);
            [no_air,el_air]=removeisolatednode(no_air,el_air);
            f_air = volface(el_air(:,1:4));
            f_air = removedupelem(f_air);
        end
        el_skin = el_skin(el_skin(:,5)==I,:);
        [no_skin,el_skin]=removeisolatednode(no_skin,el_skin);

        [f_skin] = volface(el_skin(:,1:4));
        f_skin=removedupelem(f_skin);
    end
    
    %% When the label_elem does not exist, it often indicates a failure at the generation of a coarse
    %% tetrahedral mesh. The alternative meshing pathway using decoupling is then called to make a
    %% second attempt at creating the combined tetrahedral mesh.
    if (~exist('label_elem'))&& (loop==1)
        fprintf('Initial meshing procedure failed. The option parameter might need to be adjusted. \n')
        fprintf('Activating alternative meshing pathway... \n')
        pause(2)
        continue;
    end
    
    %% The labels are given to each of the tissues
    %% WM(1) - GM(2) - CSF(3) - Bone(4) - Scalp(5) - Air(6)
    for i = 1:length(label_elem(:,1))
         if (exist('bone_n') && exist('no_air2'))
            if intriangulation(no_air2,f_air2(:,1:3),label_centroid(i,:))
                label_label(i,1) = 6;
                continue;
            end
        end
        if (exist('no_skin') && exist('no_air'))
            if intriangulation(no_air,f_air(:,1:3),label_centroid(i,:))
                label_label(i,1) = 6;
                continue;
            end
        end
        if intriangulation(wm_n,wm_f(:,1:3),label_centroid(i,:))
            label_label(i,1) = 1;
        elseif intriangulation(pial_n,pial_f(:,1:3),label_centroid(i,:))
            label_label(i,1) = 2;  
        elseif intriangulation(csf_n,csf_f(:,1:3),label_centroid(i,:))
            label_label(i,1) = 3; 
        elseif (exist('bone_n2'))
            if intriangulation(bone_n2,bone_f2(:,1:3),label_centroid(i,:))
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
    ISO2MESH_TETGENOPT = sprintf('-A -Rmpq%fa%i',cfg.ratio,cfg.maxvol);
    node(:,4) = cfg.sizefield*ones(length(node(:,1)),1);
    [brain_n,brain_el,brain_f] = surf2mesh(node,face,[],[],1.0,10,[],[],0,'tetgen1.5');

    label2 = unique(brain_el(:,5)); 
    for i = 1:length(label2)
        label_brain_el(i,1) = find(brain_el(:,5) == label2(i),1);
        label_centroid2(i,:) = meshcentroid(brain_n,brain_el(label_brain_el(i),1:4));
    end
    
    %% The labeling process is repeated for the final mesh
    %% WM(1) - GM(2) - CSF(3) - Bone(4) - Scalp(5) - Air(6)
    for i = 1:length(label_brain_el(:,1))
         if (exist('bone_n') && exist('no_air2'))
            if intriangulation(no_air2,f_air2(:,1:3),label_centroid2(i,:))
                label_label2(i,1) = 6;
                continue;
            end
        end
        if (exist('no_skin') && exist('no_air'))
            if intriangulation(no_air,f_air(:,1:3),label_centroid2(i,:))
                label_label2(i,1) = 6;
                continue;
            end
        end
        if intriangulation(wm_n,wm_f(:,1:3),label_centroid2(i,:))
            label_label2(i,1) = 1;
        elseif intriangulation(pial_n,pial_f(:,1:3),label_centroid2(i,:))
            label_label2(i,1) = 2;  
        elseif intriangulation(csf_n,csf_f(:,1:3),label_centroid2(i,:))
            label_label2(i,1) = 3; 
        elseif (exist('bone_n2'))
            if intriangulation(bone_n2,bone_f2(:,1:3),label_centroid2(i,:))
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

    for i = 1:length(brain_el(:,1))
        brain_el(i,5) = label_label2(brain_el(i,5));
    end
end
%% Relabeling step to remove layered assumptions
if cfg.relabeling == 1 && (isfield(seg2,'skull') && isfield(seg2,'scalp')) 
    centroid = meshcentroid(brain_n(:,1:3),brain_el(:,1:4));centroid = ceil(centroid);
    tag = zeros(length(brain_el(:,1)),1);
    facenb = faceneighbors(brain_el(:,1:4));
    for i = 1:length(brain_el(:,1))
        if (expandedGM(centroid(i,1),centroid(i,2),centroid(i,3))>0.5) && (brain_el(i,5) == 2)
            if seg2.scalp(centroid(i,1),centroid(i,2),centroid(i,3)) > 0.5
                brain_el(i,5) = 5;
            elseif seg2.skull(centroid(i,1),centroid(i,2),centroid(i,3)) > 0.5
                brain_el(i,5) = 4;
            else
                brain_el(i,5) = 3;
            end
            tag(i) = 1;
            for j = 1:4
                if facenb(i,j) > 0
                    tag(facenb(i,j),1) = 1; 
                end
            end
        elseif (expandedCSF(centroid(i,1),centroid(i,2),centroid(i,3))>0.5) && (brain_el(i,5) == 3)
            if seg2.scalp(centroid(i,1),centroid(i,2),centroid(i,3)) > 0.5
                brain_el(i,5) = 5;
            else
                brain_el(i,5) = 4;
            end
            brain_el(i,5) = 4;
            tag(i) = 1;
            for j = 1:4
                if facenb(i,j) > 0
                tag(facenb(i,j),1) = 1; 
                end
            end
        elseif (expandedSkull(centroid(i,1),centroid(i,2),centroid(i,3))>0.5) && (brain_el(i,5) == 4)
            brain_el(i,5) = 5;
            tag(i) = 1;
            for j = 1:4
                if facenb(i,j) > 0
                tag(facenb(i,j),1) = 1; 
                end
            end
        end
    end

    labels = zeros(length(brain_el(:,1)),4);
    labels2 = zeros(length(brain_el(:,1)),6);
    for i = 1:length(brain_el(:,1))
        for j = 1:4
            if facenb(i,j) > 0
                labels(i,j) = brain_el(facenb(i,j),5);
                labels2(i,labels(i,j)) = labels2(i,labels(i,j)) + 1; 
            else
                labels(i,j) = 0;
            end
        end
    end

    [labels(:,5),labels(:,6)] = max(labels2,[],2);
    for i = 1:length(brain_el(:,1))
        if tag(i) == 1
            if (labels(i,5) > 2) && (brain_el(i,5) ~= labels(i,6))
                brain_el(i,5) = labels(i,6);
            end
        end
    end
end
brain_el(:,5) = 6 - brain_el(:,5);
end