function [landmarks, curves, initpoints]=brain1020(node, face, initpoints, perc1, perc2, varargin)
%
% landmarks=brain1020(node, face)
%   or
% landmarks=brain1020(node, face, [], perc1, perc2)
% landmarks=brain1020(node, face, initpoints)
% landmarks=brain1020(node, face, initpoints, perc1, perc2)
% [landmarks, curves, initpoints]=brain1020(node, face, initpoints, perc1, perc2, options)
%
% compute 10-20-like scalp landmarks with user-specified density on a head mesh
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% == Input ==
%    node: full head mesh node list
%    face: full head mesh element list- a 3-column array defines face list
%          for the exterior (scalp) surface; a 4-column array defines the
%          tetrahedral mesh of the full head.
%    initpoints:(optional) one can provide the 3-D coordinates of the below
%          5 landmarks: nz, iz, lpa, rpa, cz
%          (Cz). initpoints can be a struct with the above landmark names
%          as subfield, or a 5x3 array definining these points in the above
%          mentioned order (one can use the output landmarks as initpoints)
%    perc1:(optional) the percentage of geodesic distance towards the rim of 
%          the landmarks; this is the first number of the 10-20 or 10-10 or
%          10-5 systems, in this case, it is 10 (for 10%). default is 10.
%    perc2:(optional) the percentage of geodesic distance twoards the center 
%          of the landmarks; this is the 2nd number of the 10-20 or 10-10 or
%          10-5 systems, which are 20, 10, 5, respectively, default is 20
%    options: one can add additional 'name',value pairs to the function
%          call to provide additional control. Supported optional names
%          include
%           'display' : [1] or 0, if set to 1, plot landmarks and curves
%           'cztol' : [1e-6], the tolerance for searching cz that bisects
%                   saggital and coronal reference curves
%           'maxcziter' : [10] the maximum number of iterations to update
%                   cz to bisect both cm and sm curves
%           'baseplane' : [1] or 0, if set to 1, create the reference
%                   curves along the primary control points (nz,iz,lpa,rpa)
%           'minangle' : [0] if set to a positive number, this specifies
%                   the minimum angle (radian) between adjacent segments in
%                   the reference curves to avoid sharp turns (such as the
%                   dipping near ear cavities), this parameter will be
%                   passed to polylinesimplify to simplify the curve first.
%                   Please be note that the landmarks generated with
%                   simplified curves may not land exactly on the surface.
%
% == Output ==
%    landmarks: a structure storing all computed landmarks. The subfields
%          include two sections: 
%          1) 'nz','iz','lpa','rpa','cz': individual 3D positions defining
%             the 5 principle reference points: nasion (nz), inion (in),
%             left-pre-auricular-point (lpa), right-pre-auricular-point 
%             (rpa) and vertex (cz)
%          2) landmarks along specific cross-sections, each cross section
%             may contain more than 1 position. The cross-sections are
%             named in the below format:
%             2.1: a fieldname starting from 'c', 's' and 'a' indicates
%                  the cut is along coronal, saggital and axial directions,
%                  respectively;
%             2.2: a name starting from 'pa' indicates the cut is along the
%                  axial plane acrossing principle reference points
%             2.3: the following letter 'm', 'a','p' suggests the 'medial',
%                  'anterior' and 'posterior', respectively
%             2.4: the last letter 'l' or 'r' suggests the 'left' and
%                  'right' side, respectively
%             2.5: non-medial coronal cuts are divided into two groups, the
%                  anterior group (ca{l,r}) and the posterior group
%                  (cp{lr}), with a number indicates the node spacing 
%                  stepping away from the medial plane.
%
%             for example, landmarks.cm refers to the landmarks along the
%                  medial-coronal plane, anterior-to-posteior order
%             similarly, landmarks.cpl_3 refers to the landmarks along the
%                  coronal (c) cut plane located in the posterior-left side
%                  of the head, with 3 saggital landmark spacing from the
%                  medial-coronal reference curve.
%    curves: a structure storing all computed cross-section curves. The
%             subfields are named similarly to landmarks, except that
%             landmarks stores the 10-? points, and curves stores the
%             detailed cross-sectional curves
%    initpoints: a 5x3 array storing the principle reference points in the
%             orders of 'nz','iz','lpa','rpa','cz'
%
% 
% == Example ==
% See brain2mesh/examples/SPM_example_brain.m for an example
% https://github.com/fangq/brain2mesh/blob/master/examples/SPM_example_brain.m
%
% == Dependency ==
% This function requires a pre-installed Iso2Mesh Toolbox
%  Download URL: http://github.com/fangq/iso2mesh
%  Website: http://iso2mesh.sf.net
%
% == Reference ==
% If you use this function in your publication, the authors of this toolbox
% apprecitate if you can cite the below paper
%
%  Anh Phong Tran, Shijie Yan and Qianqian Fang, "Improving model-based
%  fNIRS analysis using mesh-based anatomical and light-transport models,"
%  Neurophotonics, 7(1), 015008, URL: http://dx.doi.org/10.1117/1.NPh.7.1.015008
%
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

if(nargin<2)
    error('one must provide a head-mesh to call this function');
end

if(isempty(node) || isempty(face) || size(face,2)<=2 || size(node,2)<3)
    error('input node must have 3 columns, face must have at least 3 columns');
end

if(nargin<5)
    perc2=20;
    if(nargin<4)
        perc1=10;
        if(nargin<3)
            initpoints=[];
        end
    end
end

% parse user options
opt=varargin2struct(varargin{:});

showplot=jsonopt('display',1,opt);
baseplane=jsonopt('baseplane',1,opt);
tol=jsonopt('cztol',1e-6,opt);
dosimplify=jsonopt('minangle',0,opt);
maxcziter=jsonopt('maxcziter',10,opt);

% convert initpoints input to a 5x3 array
if(isstruct(initpoints))
    initpoints=struct('nz', initpoints.nz(:).','iz', initpoints.iz(:).',...
                     'lpa',initpoints.lpa(:).','rpa',initpoints.rpa(:).',...
		             'cz', initpoints.cz(:).');
    landmarks=initpoints;
    if(exist('struct2array','file'))
        initpoints=struct2array(initpoints);
        initpoints=reshape(initpoints(:),3,length(initpoints(:))/3)';
    else
        initpoints=[initpoints.nz(:).'; initpoints.iz(:).';initpoints.lpa(:).';initpoints.rpa(:).';initpoints.cz(:).'];
    end
end

% convert tetrahedral mesh into a surface mesh
if(size(face,2)>=4)
    face=volface(face(:,1:4));
end

% remove nodes not located in the surface
[node,face]=removeisolatednode(node,face);

% if initpoints is not sufficient, ask user to interactively select nz, iz, lpa, rpa and cz first
if(isempty(initpoints) || size(initpoints,1)<5)
    hf=figure;
    plotmesh(node,face);
    set(hf,'userdata',initpoints);
    if(~isempty(initpoints))
        hold on;
        plotmesh(initpoints,'gs', 'LineWidth',4);
    end
    idx=size(initpoints,1)+1;
    landmarkname={'Nasion','Inion','Left-pre-auricular-point','Right-pre-auricular-point','Vertex/Cz','Done'};
    title(sprintf('Rotate the mesh, select data cursor, click on P%d: %s',idx, landmarkname{idx}));
    rotate3d('on');
    set(datacursormode(hf),'UpdateFcn',@myupdatefcn);
end

% wait until all 5 points are defined
if(exist('hf','var'))
    try
        while(size(get(hf,'userdata'),1)<5)
            pause(0.1);
        end
    catch
        error('user aborted');
    end
    datacursormode(hf,'off');
    initpoints=get(hf,'userdata');
    close(hf);
end

if(showplot)
    disp(initpoints);
end

% save input initpoints to landmarks output, cz is not finalized
if(size(initpoints,1)>=5)
    landmarks=struct('nz', initpoints(1,:),'iz', initpoints(2,:),...
                     'lpa',initpoints(3,:),'rpa',initpoints(4,:),...
		     'cz', initpoints(5,:));
end

% at this point, initpoints contains {nz, iz, lpa, rpa, cz0}
% plot the head mesh
if(showplot)
    figure;
    hp=plotmesh(node,face,'facealpha',0.6,'facecolor',[1 0.8 0.7]);
    if(~isoctavemesh)
        set(hp,'linestyle','none');
    end
    camlight;
    lighting gouraud
    hold on;
end

lastcz=[1 1 1]*inf;
cziter=0;

%% Find cz that bisects cm and sm curves within a tolerance, using UI 10-10 approach
while(norm(initpoints(5,:)-lastcz)>tol && cziter<maxcziter)
    %% Step 1: nz, iz and cz0 to determine saggital reference curve
    
    nsagg=slicesurf(node, face, initpoints([1,2,5],:));
    %% Step 1.1: get cz1 as the mid-point between iz and nz
    [slen, nsagg]=polylinelen(nsagg, initpoints(1,:), initpoints(2,:), initpoints(5,:));
    if(dosimplify)
        [nsagg, slen]=polylinesimplify(nsagg,dosimplify);
    end
    [idx, weight, cz]=polylineinterp(slen, sum(slen)*0.5, nsagg);
    initpoints(5,:)=cz(1,:);

    %% Step 1.2: lpa, rpa and cz1 to determine coronal reference curve, update cz1
    curves.cm=slicesurf(node, face, initpoints([3,4,5],:));
    [len, curves.cm]=polylinelen(curves.cm, initpoints(3,:), initpoints(4,:), initpoints(5,:));
    if(dosimplify)
        [curves.cm, len]=polylinesimplify(curves.cm,dosimplify);
    end
    [idx, weight, coro]=polylineinterp(len, sum(len)*0.5, curves.cm);
    lastcz=initpoints(5,:);
    initpoints(5,:)=coro(1,:);
    cziter=cziter+1;
    if(showplot)
        fprintf('cz iteration %d error %e\n',cziter, norm(initpoints(5,:)-lastcz));
    end
end

% set the finalized cz to output
landmarks.cz=initpoints(5,:);
[idx, weight, coro]=polylineinterp(len, sum(len)*(perc1:perc2:(100-perc1))*0.01, curves.cm);
landmarks.cm=coro;                 % t7, c3, cz, c4, t8

if(showplot)
    disp(initpoints);
end

%% Step 2: subdivide saggital (sm) and coronal (cm) ref curves

curves.sm=slicesurf(node, face, initpoints([1,2,5],:));
[slen, curves.sm]=polylinelen(curves.sm, initpoints(1,:), initpoints(2,:), initpoints(5,:));
if(dosimplify)
    [curves.sm, slen]=polylinesimplify(curves.sm,dosimplify);
end
[idx, weight, sagg]=polylineinterp(slen, sum(slen)*(perc1:perc2:(100-perc1))*0.01, curves.sm);
landmarks.sm=sagg;           % fpz, fz, cz, pz, oz

%% Step 3: fpz, t7 and oz to determine left 10% axial reference curve

[landmarks.aal, curves.aal, landmarks.apl, curves.apl]=slicesurf3(node,face,landmarks.sm(1,:), landmarks.cm(1,:), landmarks.sm(end,:),perc2*2);

%% Step 4: fpz, t8 and oz to determine right 10% axial reference curve

[landmarks.aar, curves.aar, landmarks.apr, curves.apr]=slicesurf3(node,face,landmarks.sm(1,:), landmarks.cm(end,:),landmarks.sm(end,:), perc2*2);


%% show plots of the landmarks
if(showplot)
    plotmesh(curves.sm,'r-','LineWidth',1);
    plotmesh(curves.cm,'g-','LineWidth',1);
    plotmesh(curves.aal,'k-','LineWidth',1);
    plotmesh(curves.aar,'k-','LineWidth',1);
    plotmesh(curves.apl,'b-','LineWidth',1);
    plotmesh(curves.apr,'b-','LineWidth',1);

    plotmesh(landmarks.sm,'ro','LineWidth',2);
    plotmesh(landmarks.cm,'go','LineWidth',2);
    plotmesh(landmarks.aal,'ko','LineWidth',2);
    plotmesh(landmarks.aar,'mo','LineWidth',2);
    plotmesh(landmarks.apl,'ko','LineWidth',2);
    plotmesh(landmarks.apr,'mo','LineWidth',2);
end

%% Step 5: computing all anterior coronal cuts, moving away from the medial cut (cm) toward frontal

idxcz=closestnode(landmarks.sm,landmarks.cz);

skipcount=max(floor(10/perc2),1);

for i=1:size(landmarks.aal,1)-skipcount
    step=(perc2*25)*0.1*(1+((perc2<20 + perc2<10) && i==size(landmarks.aal,1)-skipcount));
    [landmarks.(sprintf('cal_%d',i)), leftpart, landmarks.(sprintf('car_%d',i)), rightpart]=slicesurf3(node,face,landmarks.aal(i,:), landmarks.sm(idxcz-i,:), landmarks.aar(i,:),step);
    if(showplot)
        plotmesh(leftpart,'k-','LineWidth',1);
        plotmesh(rightpart,'k-','LineWidth',1);

        plotmesh(landmarks.(sprintf('cal_%d',i)),'yo','LineWidth',2);
        plotmesh(landmarks.(sprintf('car_%d',i)),'co','LineWidth',2);
    end
end

%% Step 6: computing all posterior coronal cuts, moving away from the medial cut (cm) toward occipital

for i=1:size(landmarks.apl,1)-skipcount
    step=(perc2*25)*0.1*(1+((perc2<20 + perc2<10) && i==size(landmarks.apl,1)-skipcount));
    [landmarks.(sprintf('cpl_%d',i)), leftpart, landmarks.(sprintf('cpr_%d',i)), rightpart]=slicesurf3(node,face,landmarks.apl(i,:), landmarks.sm(idxcz+i,:), landmarks.apr(i,:),step);
    if(showplot)
        plotmesh(leftpart,'k-','LineWidth',1);
        plotmesh(rightpart,'k-','LineWidth',1);

        plotmesh(landmarks.(sprintf('cpl_%d',i)),'yo','LineWidth',2);
        plotmesh(landmarks.(sprintf('cpr_%d',i)),'co','LineWidth',2);
    end
end

%% Step 7: create the axial cuts across priciple ref. points: left: nz, lpa, iz, right: nz, rpa, iz

if(baseplane && perc2<=10)
    [landmarks.paal, curves.paal, landmarks.papl, curves.papl]=slicesurf3(node,face,landmarks.nz, landmarks.lpa, landmarks.iz, perc2*2);
    [landmarks.paar, curves.paar, landmarks.papr, curves.papr]=slicesurf3(node,face,landmarks.nz, landmarks.rpa, landmarks.iz, perc2*2);
    if(showplot)
        plotmesh(curves.paal,'k-','LineWidth',1);
        plotmesh(curves.paar,'k-','LineWidth',1);
        plotmesh(curves.papl,'k-','LineWidth',1);
        plotmesh(curves.papr,'k-','LineWidth',1);

        plotmesh(landmarks.paal,'yo','LineWidth',2);
        plotmesh(landmarks.papl,'co','LineWidth',2);
        plotmesh(landmarks.paar,'yo','LineWidth',2);
        plotmesh(landmarks.papr,'co','LineWidth',2);
    end
end

%%-------------------------------------------------------------------------------------
% helper functions
%--------------------------------------------------------------------------------------

% the respond function when a data-cursor tip to popup

function txt=myupdatefcn(empt,event_obj)
pt=get(gcf,'userdata');
pos = get(event_obj,'Position');
%idx=  get(event_obj,'DataIndex');
txt = {['x: ',num2str(pos(1))],...
       ['y: ',num2str(pos(2))],['z: ',num2str(pos(3))]};
if(~isempty(pt) && ismember(pos,pt,'rows'))
    return;
end
targetup=get(get(event_obj,'Target'),'parent');
idx=size(pt,1)+2;
landmarkname={'Nasion','Inion','Left-pre-auricular-point','Right-pre-auricular-point','Vertex/Cz','Done'};
title(sprintf('Rotate the mesh, select data cursor, click on P%d: %s',idx, landmarkname{idx}));
set(targetup,'userdata',struct('pos',pos));
pt=[pt;pos];
if(size(pt,1)<6)
     set(gcf,'userdata',pt);
end
hpt=findobj(gcf,'type','line');
if(~isempty(hpt))
    set(hpt,'xdata',pt(:,1),'ydata',pt(:,2),'zdata',pt(:,3));
else
    hold on;
    plotmesh([pt;pos],'gs', 'LineWidth',4);
end
disp(['Adding landmark ' landmarkname{idx-1} ':' txt]);
datacursormode(gcf,'off');
rotate3d('on');
