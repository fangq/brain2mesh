function [brain_n,brain_el,brain_f] = brain2mesh(seg,varargin)
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
%          cfg.radbound.{wm,gm,csf,skull,scalp}
%             Radius of the Delaunay sphere used in the sampling the surfaces.
%             Default values are 1.7, 1.7, 2, 2.5 and 3, respectively (reference values for 1x1x1mm^3)
%             Scale proportionally for denser volumes. Lower values correspond to denser, higher
%             fidelity surface extraction, but also results in denser meshes.
%          cfg.maxnode: [100000] - when  the value cfg.sampling__ creates surfaces that are too 
%             dense. This limits the maximum of number of nodes extracted for a given surface.
%          cfg.maxvol: [100] indicates the volumetric maximum size of elements
%             Lowering this value helps with obtaining a denser tetrahedral
%             mesh. For dense meshes, values close to 3-5 are recommended.
%          cfg.smooth: [0] - number of iterations to smooth each tissue surface
%          cfg.ratio: [1.414] radius-edge ratio. Lower values increase 
%             the quality of tetrahedral elements, but results in denser meshes
%          cfg.dorelabel: [0] or 1 - This step removes most of the assumptions 
%             created by the layered meshing workflow. Currently only works if all five tissue types are present.
%             When deactivated, a 1 voxel length gap is assumed between each of the tissue layers.
%          cfg.doairseg: 0 or [1]. Within the skull layer, additional segmentations can be found. 
%             By default, these regions are merged to the skull because they can be ambiguously be
%             vessels or air. When the option 1 is chosen, these regions are labeled as air instead.
%          cfg.dotruncate: 0 or [1]. by default, the mesh is truncated a few pixel in the 
%             +z direction below the lowest pixel containing CSF to focus on the brain areas. 
%             A value of 0 gives a complete head mesh. 
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
if nargin == 0
    help brain2mesh
    return;
end

if(~exist('v2m','file'))
    error('Missing dependency. You must download and addpath to brain2mesh.m, URL: https://github.com/fangq/brain2mesh')
end
if(~exist('imfill','file'))
    error('Missing dependency. You must install MATLAB image processing toolbox')
end
if(~exist('intriangulation','file'))
    error('Missing dependency. You must download and addpath to intriangulation.m, URL: https://www.mathworks.com/matlabcentral/fileexchange/43381-intriangulation-vertices-faces-testp-heavytest')
end

density=struct('wm',2,'gm',2,'csf',5,'skull',4,'scalp',8);
adaptiveness=struct('wm',1,'gm',1,'csf',1,'skull',1,'scalp',1);
cfg=varargin2struct(varargin{:});

radbound=jsonopt('radbound',density,cfg);
distbound=jsonopt('distbound',adaptiveness,cfg);
qratio=jsonopt('ratio',1.414,cfg);
%sizefield=jsonopt('sizefield',100,cfg);
maxvol=jsonopt('maxvol',100.414,cfg);
maxnode=jsonopt('maxnode',100000,cfg);
dotruncate=jsonopt('dotruncate',1,cfg);
dorelabel=jsonopt('dorelabel',0,cfg);
doairseg=jsonopt('doairseg',1,cfg);
threshold=jsonopt('threshold',0.5,cfg);
smooth=jsonopt('smooth',10,cfg);

segname=fieldnames(density);

if isstruct(seg)
    seg2 = seg;
elseif(ndim(seg)==4)
    for i=1:size(seg,4)
        seg2.(segname{i})=seg(:,:,:,i);
    end
else
    fprintf('This seg input is currently not supported \n')
end

opt=struct;

for i = 1:size(fieldnames(seg2))
    opt(i).maxnode = maxnode; 
    if(isfield(radbound,segname{i}))
        opt(i).radbound= radbound.(segname{i});
    end
    %opt(i).distbound= distbound.(segname{i});
end

cube3=true(3,3,3);

%% Pre-processing steps to create separations between the tissues in the
% volume space

dim = size(seg2.wm);
seg2.wm = imfill(seg2.wm,'holes');
p_wm = seg2.wm;
p_pial = p_wm+seg2.gm;
p_pial = max(p_pial,imdilate(p_wm,cube3));
p_pial = imfill(p_pial,'holes');
p_csf = p_pial+seg2.csf;
p_csf(p_csf>1) = 1;
p_csf = max(p_csf,imdilate(p_pial,cube3));

expandedGM = p_pial - seg2.wm - seg2.gm;
expandedCSF = p_csf - seg2.wm - seg2.gm - seg2.csf - expandedGM;
expandedGM = imdilate(expandedGM,cube3);
expandedCSF = imdilate(expandedCSF,cube3);

if isfield(seg2,'skull') && isfield(seg2,'scalp')
    p_bone = p_csf + seg2.skull;
    p_bone(p_bone>1) = 1;
    p_bone = max(p_bone,imdilate(p_csf,cube3));
    p_skin = p_bone + seg2.scalp;
    p_skin(p_skin>1) = 1;
    p_skin = max(p_skin,imdilate(p_bone,cube3));
	expandedSkull = p_bone - seg2.wm - seg2.gm - seg2.csf - seg2.skull - expandedCSF - expandedGM;
    expandedSkull = imdilate(expandedSkull,cube3);
elseif isfield(seg2,'scalp') && ~isfield(seg2,'skull')
    p_skin = p_csf + seg2.scalp;
    p_skin(p_skin>1) = 1;
    p_skin = max(p_skin,imdilate(p_csf,cube3));
elseif isfield(seg2,'skull') && ~isfield(seg2,'scalp')
    p_bone = p_csf + seg2.skull;
    p_bone(p_bone>1) = 1;
    p_bone = max(p_bone,imdilate(p_csf,cube3));
end

%% Grayscale/Binary extractions of the surface meshes for the different
% tissues

[wm_n,wm_f] = v2s(p_wm,threshold,opt(1),'cgalsurf');
[pial_n,pial_f] = v2s(p_pial,threshold,opt(2),'cgalsurf');
[csf_n,csf_f] = v2s(p_csf,threshold,opt(3),'cgalsurf');

% [wm_n,wm_f]=meshcheckrepair(wm_n,wm_f(:,1:3),'isolated');
% [pial_n,pial_f]=meshcheckrepair(pial_n,pial_f(:,1:3),'isolated');
% [csf_n,csf_f]=meshcheckrepair(csf_n,csf_f(:,1:3),'isolated');
% 
% wm_n=sms(wm_n,wm_f,smooth,0.5,'lowpass');
% pial_n=sms(pial_n,pial_f,smooth,0.5,'lowpass');
% csf_n=sms(csf_n,csf_f,smooth,0.5,'lowpass');

if isfield(seg2,'skull')
    optskull=struct('radbound',radbound.skull,'maxnode',maxnode);
    [bone_n,bone_f] = v2s(p_bone,threshold,optskull,'cgalsurf');
    %[bone_n,bone_f]=meshcheckrepair(bone_n,bone_f(:,1:3),'isolated');
    %bone_n=sms(bone_n,bone_f,smooth,0.5,'lowpass');

    [bone_node,el_bone] = s2m(bone_n,bone_f,1.0,maxvol,'tetgen1.5',[],[],'-A');
    for i = 1:length(unique(el_bone(:,5)))
        vol_bone(i) = sum(elemvolume(bone_node,el_bone(el_bone(:,5)==i,1:4)));
    end
    [maxval,I] = max(vol_bone);
    if (length(unique(el_bone(:,5)))>1)
        no_air2 = bone_node; el_air2 = el_bone(el_bone(:,5)~=I,:);
        [no_air2,el_air2]=removeisolatednode(no_air2,el_air2);
        f_air2 = volface(el_air2(:,1:4));
    end
    bone_n2 = bone_node;
    [bone_f2] = volface(el_bone(:,1:4));
    bone_f2 = removedupelem(bone_f2);
    [bone_n2,bone_f2]=removeisolatednode(bone_n2,bone_f2);
    if doairseg == 0
        bone_n = bone_n2; bone_f = bone_f2;
    end
end
if isfield(seg2,'scalp')
    optscalp=struct('radbound',radbound.scalp,'maxnode',maxnode);
    [skin_n,skin_f] = v2s(p_skin,threshold,optscalp,'cgalsurf');
    %[skin_n,skin_f]=meshcheckrepair(skin_n,skin_f(:,1:3),'isolated');
    %skin_n=sms(skin_n,skin_f,smooth,0.5,'lowpass');
end

%% Main loop for the meshing pipeline to combine the individual surface
% meshes or each of the tissues and to generate the detailed 3D tetrahedral
% mesh of the brain/head

for loop = 1:2
    %% If the first pass fails, a second pass is called using the decoupled function
    % to eliminate intersections between surface meshes
    if (loop==2) && (exist('label_elem'))
        continue;
    end
    if (loop==2) && (~exist('label_elem'))
        [bone_n,bone_f] = surfboolean(bone_n(:,1:3),bone_f(:,1:3),'decouple',skin_n(:,1:3),skin_f(:,1:3));
        [csf_n,csf_f] = surfboolean(csf_n(:,1:3),csf_f(:,1:3),'decouple',bone_n(:,1:3),bone_f(:,1:3));
        [pial_n,pial_f] = surfboolean(pial_n(:,1:3),pial_f(:,1:3),'decouple',csf_n(:,1:3),csf_f(:,1:3));
        [wm_n,wm_f] = surfboolean(wm_n(:,1:3),wm_f(:,1:3),'decouple',pial_n(:,1:3),pial_f(:,1:3));
    end
    [surf_n,surf_f] = surfboolean(wm_n(:,1:3),wm_f(:,1:3),'resolve',pial_n,pial_f);
    [surf_n,surf_f] = surfboolean(surf_n,surf_f,'resolve',csf_n,csf_f);
    if isfield(seg2,'skull')
        [surf_n,surf_f] = surfboolean(surf_n,surf_f,'resolve',bone_n,bone_f);
    end
    if isfield(seg2,'scalp')
        [surf_n,surf_f] = surfboolean(surf_n,surf_f,'resolve',skin_n,skin_f);
    end
    final_surf_n=surf_n;
    final_surf_f=surf_f;
    %% If the whole head option is deactivated, the cut is made at the base of the brain using a box cutting
    if (dotruncate==1)
        dim=max(surf_n);
        dim2 = min(csf_n);
        [nbox,fbox,ebox]=meshabox([-1 -1 dim2(3)+4.1],[dim(1)+1 dim(2)+1 dim(3)+1],500);
        fbox=volface(ebox);
        [nbox,fbox]=removeisolatednode(nbox,fbox);
        [final_surf_n,final_surf_f] = surfboolean(nbox,fbox(:,[1 3 2]),'first',surf_n,surf_f);
    end
    
    %% Generates a coarse tetrahedral mesh of the combined tissues
    try
        [final_n,final_e] = s2m(final_surf_n,final_surf_f,1.0,maxvol,'tetgen1.5',[],[],'-A');
    catch
        fprintf('volumetric mesh generation failed, returning the intermediate surface model only');
        brain_n=final_surf_n;
        brain_f=final_surf_f;
        brain_el=[];
        return;
    end

    %% Removes the elements that are part of the box, but not the brain/head
    if (dotruncate == 1)
        [maxval, M] = max(final_n);
        k = find(final_e(:,1:4)==M(3),1);
        final_e = final_e(final_e(:,5)~=final_e(rem(k,length(final_e(:,1))),5),:);
        [final_n,final_e]=removeisolatednode(final_n,final_e);
    end

    %% Here the labels created through the coarse mesh generated through Tetgen are saved
    % with the centroid of one of the elements for intriangulation seg(:,:,:,1)testing later
    [label, label_elem] = unique(final_e(:,5)); 
    label_centroid=meshcentroid(final_n,final_e(label_elem,1:4));

    if isfield(seg2,'scalp')
        [no_skin,el_skin] = s2m(skin_n,skin_f,1.0,maxvol,'tetgen1.5',[],[],'-A');
        for i = 1:length(unique(el_skin(:,5)))
            vol_skin(i) = sum(elemvolume(no_skin,el_skin(el_skin(:,5)==i,1:4)));
        end
        [maxval,I] = max(vol_skin);
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
    % tetrahedral mesh. The alternative meshing pathway using decoupling is then called to make a
    % second attempt at creating the combined tetrahedral mesh.
    if (~exist('label_elem'))&& (loop==1)
        fprintf('Initial meshing procedure failed. The option parameter might need to be adjusted. \n')
        fprintf('Activating alternative meshing pathway... \n')
        pause(2)
        continue;
    end
    
    %% The labels are given to each of the tissues
    % WM(1) - GM(2) - CSF(3) - Bone(4) - Scalp(5) - Air(6)
    newlabel = zeros(length(label_elem),1);
    if (exist('bone_n') && exist('no_air2'))
        newlabel= intriangulation(no_air2,f_air2(:,1:3),label_centroid);
    end
    if (exist('no_skin') && exist('no_air'))
        newlabel= newlabel | intriangulation(no_air,f_air(:,1:3),label_centroid);
    end
    newlabel=double(newlabel);
    idx=find(newlabel==0);
    
    newtag=zeros(length(idx),1);
    newtag=intriangulation(wm_n,wm_f(:,1:3),label_centroid(idx,:))*6;
    newtag=max(newtag,intriangulation(pial_n,pial_f(:,1:3),label_centroid(idx,:))*5);
    newtag=max(newtag,intriangulation(csf_n,csf_f(:,1:3),label_centroid(idx,:))*4);
    if(exist('bone_n2'))
        newtag=max(newtag,intriangulation(bone_n2,bone_f2(:,1:3),label_centroid(idx,:))*3);
    end
    if(exist('no_skin'))
        newtag=max(newtag,intriangulation(no_skin,f_skin(:,1:3),label_centroid(idx,:))*2);
    end
    newlabel(idx)=newtag;
    newlabel=7-newlabel;
    final_e(:,5)=newlabel(final_e(:,5));

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
    cmdopt = sprintf('-A -pq%fa%f',qratio,maxvol);
    %node(:,4) = sizefield(:).*ones(length(node(:,1)),1);
    [brain_n,brain_el,brain_f] = s2m(node,face,1.0,maxvol,'tetgen1.5',[],[],cmdopt);

    [label2, label_brain_el] = unique(brain_el(:,5)); 
    label_centroid2=meshcentroid(brain_n,brain_el(label_brain_el,1:4));
    
    %% The labeling process is repeated for the final mesh
    % WM(1) - GM(2) - CSF(3) - Bone(4) - Scalp(5) - Air(6)
    newlabel = zeros(length(label_brain_el),1);
    if (exist('bone_n') && exist('no_air2'))
        newlabel= intriangulation(no_air2,f_air2(:,1:3),label_centroid2);
    end
    if (exist('no_skin') && exist('no_air'))
        newlabel= newlabel | intriangulation(no_air,f_air(:,1:3),label_centroid2);
    end
    newlabel=double(newlabel);
    idx=find(newlabel==0);
    
    newtag=zeros(length(idx),1);
    newtag=intriangulation(wm_n,wm_f(:,1:3),label_centroid2(idx,:))*6;
    newtag=max(newtag,intriangulation(pial_n,pial_f(:,1:3),label_centroid2(idx,:))*5);
    newtag=max(newtag,intriangulation(csf_n,csf_f(:,1:3),label_centroid2(idx,:))*4);
    if(exist('bone_n2'))
        newtag=max(newtag,intriangulation(bone_n2,bone_f2(:,1:3),label_centroid2(idx,:))*3);
    end
    if(exist('no_skin'))
        newtag=max(newtag,intriangulation(no_skin,f_skin(:,1:3),label_centroid2(idx,:))*2);
    end
    newlabel(idx)=newtag;
    newlabel=7-newlabel;
    brain_el(:,5)=newlabel(brain_el(:,5));
end

%% Relabeling step to remove layered assumptions
if dorelabel == 1 && (isfield(seg2,'skull') && isfield(seg2,'scalp')) 
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
