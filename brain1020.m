function [landmarks, initpoints, lines]=brain1020(node, elem, initpoints, perc1, perc2, varargin)
%
% landmarks=brain1020(node, elem)
%   or
% landmarks=brain1020(node, elem, initpoints)
% landmarks=brain1020(node, elem, initpoints, 10, 10)
%
% compute 10-20-like scalp landmarks with user-specified density on a head mesh
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    node: full head mesh node list
%    elem: full head mesh element list- a 3-column array defines face list
%          for the exterior (scalp) surface; a 4-column array defines the
%          tetrahedral mesh of the full head.
%    initpoints:(optional) one can provide the 3-D coordinates of the below
%          5 landmarks: nz, iz, lpa, rpa, cz
%          (Cz). initpoints can be a struct with the above landmark names
%          as subfield, or a 5x3 array definining these points in the above
%          mentioned order.
%    perc1:(optional) the percentage of geodesic distance towards the rim of 
%          the landmarks; this is the first number of the 10-20 or 10-10 or
%          10-5 systems, in this case, it is 10 (for 10%). default is 10.
%    perc2:(optional) the percentage of geodesic distance twoards the center 
%          of the landmarks; this is the 2nd number of the 10-20 or 10-10 or
%          10-5 systems, which are 20, 10, 5, respectively, default is 20
%
% output:
%    landmarks: a 3-column array defining the head surface landmark
%          positions at the requested spacing.
%
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

if(isempty(node) || isempty(elem) || size(elem,2)<=2 || size(node,2)<3)
    error('input node must have 3 columns, elem must have at least 3 columns');
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

opt=varargin2struct(varargin{:});

if(isstruct(initpoints))
    initpoints=orderfields(initpoints,{'nz', 'iz', 'lpa', 'rpa', 'cz'});
    landmarks=initpoints;
    initpoints=struct2array(initpoints);
    initpoints=reshape(initpoints(:),3,length(initpoints(:))/3)';
elseif(size(initpoints,1)>=5)
    landmarks=struct('nz', initpoints(1,:),'iz', initpoints(2,:),...
                     'lpa',initpoints(3,:),'rpa',initpoints(4,:),...
		             'cz', initpoints(5,:));
end

if(size(elem,2)>=4)
    elem=volface(elem(:,1:4));
end

if(isempty(initpoints) && size(initpoints,1)<5)
    hf=figure;
    plotmesh(node,elem,'facecolor','b');
    camlight; lighting phong
    datacursormode(hf,'on');
    set(datacursormode(hf),'UpdateFcn',@myupdatefcn);
end

if(exist('hf','var'))
    while(size(get(hf,'userdata'),1)<5)
        pause(0.1);
    end
    datacursormode(hf,'off');
    initpoints=get(hf,'userdata');
    close(hf);
end


% at this point, initpoints contains {nz, iz, lpa, rpa, cz0}

if(jsonopt('display',1,opt))
    figure;
    plotmesh(node,elem,'linestyle','none','facealpha',0.3,'facecolor','b');
    camlight; lighting phong
    hold on;
end
%% Step b: nz, iz and cz0 to determine saggital reference curve
nsagg=slicesurf(node, elem, initpoints([1,2,5],:));

%% Step c1: get cz1 as the mid-point between iz and nz
[slen, nsagg]=polylinelen(nsagg, initpoints(1,:), initpoints(2,:), initpoints(5,:));
[idx, weight, cz]=polylineinterp(slen, sum(slen)*0.5, nsagg);
initpoints(5,:)=cz(1,:);

%% Step c2: lpa, rpa and cz to determine coronal reference curve, get true cz
lines.ncoro=slicesurf(node, elem, initpoints([3,4,5],:));
[len, lines.ncoro]=polylinelen(lines.ncoro, initpoints(3,:), initpoints(4,:), initpoints(5,:));
[idx, weight, coro]=polylineinterp(len, sum(len)*[50 perc1:perc2:(100-perc1)]*0.01, lines.ncoro);

initpoints(5,:)=coro(1,:);
landmarks.cz=coro(1,:);      % using UI 10-10 approach
landmarks.c0=coro(2:end,:);  % t7, c3, cz, c4, t8

initpoints

%% Step d/e: subdivide saggital and coronal ref curves

lines.nsagg=slicesurf(node, elem, initpoints([1,2,5],:));
[slen, lines.nsagg]=polylinelen(lines.nsagg, initpoints(1,:), initpoints(2,:), initpoints(5,:));
[idx, weight, sagg]=polylineinterp(slen, sum(slen)*[perc1:perc2:(100-perc1)]*0.01, lines.nsagg);
landmarks.s0=sagg;           % fpz, fz, cz, pz, oz

%% Step f,h,i: fpz, t7 and oz to determine left 10% axial reference curve

[landmarks.nal, lines.nalaxis, landmarks.npl, lines.nplaxis]=slicesurf3(node,elem,landmarks.s0(1,:), landmarks.c0(1,:), landmarks.s0(end,:),perc2*2);

%% Step g: fpz, t8 and oz to determine right 10% axial reference curve

[landmarks.nar, lines.naraxis, landmarks.npr, lines.npraxis]=slicesurf3(node,elem,landmarks.s0(1,:), landmarks.c0(end,:),landmarks.s0(end,:), perc2*2);


%% debug
if(jsonopt('display',1,opt))
    plotmesh(lines.nsagg,'r-','LineWidth',1);
    plotmesh(lines.ncoro,'g-','LineWidth',1);
    plotmesh(lines.nalaxis,'k-','LineWidth',1);
    plotmesh(lines.nplaxis,'b-','LineWidth',1);
    plotmesh(lines.naraxis,'k-','LineWidth',1);
    plotmesh(lines.npraxis,'b-','LineWidth',1);

    plotmesh(landmarks.s0,'ro','MarkerSize',2);
    plotmesh(landmarks.c0,'go','MarkerSize',2);
    plotmesh(landmarks.nal,'ko','MarkerSize',2);
    plotmesh(landmarks.nar,'mo','MarkerSize',2);
    plotmesh(landmarks.npl,'ko','MarkerSize',2);
    plotmesh(landmarks.npr,'mo','MarkerSize',2);
end

%% Step j: f8, fz and f7 to determine front coronal cut

idxcz=closestnode(landmarks.s0,landmarks.cz);

for i=1:size(landmarks.nal,1)-1
    step=perc2/10*25*(1+(perc2<20 && i==size(landmarks.nal,1)-1));
    [landmarks.(sprintf('cal_%d',i)), leftpart, landmarks.(sprintf('car_%d',i)), rightpart]=slicesurf3(node,elem,landmarks.nal(i,:), landmarks.s0(idxcz-i,:), landmarks.nar(i,:),step);
    if(jsonopt('display',1,opt))
        plotmesh(leftpart,'k-','LineWidth',2);
        plotmesh(rightpart,'k-','LineWidth',2);

        plotmesh(landmarks.(sprintf('cal_%d',i)),'yo','MarkerSize',2);
        plotmesh(landmarks.(sprintf('car_%d',i)),'co','MarkerSize',2);
    end
end

for i=1:size(landmarks.npl,1)-1
    step=perc2/10*25*(1+(perc2<20 && i==size(landmarks.nal,1)-1));
    [landmarks.(sprintf('cpl_%d',i)), leftpart, landmarks.(sprintf('cpr_%d',i)), rightpart]=slicesurf3(node,elem,landmarks.npl(i,:), landmarks.s0(idxcz+i,:), landmarks.npr(i,:),step);
    if(jsonopt('display',1,opt))
        plotmesh(leftpart,'k-','LineWidth',1);
        plotmesh(rightpart,'k-','LineWidth',1);

        plotmesh(landmarks.(sprintf('cpl_%d',i)),'yo','MarkerSize',2);
        plotmesh(landmarks.(sprintf('cpr_%d',i)),'co','MarkerSize',2);
    end
end


%% helper functions
% the respond function when there is a data-tip to popup

function txt=myupdatefcn(empt,event_obj)
pos = get(event_obj,'Position');
%idx=  get(event_obj,'DataIndex');
txt = {['x: ',num2str(pos(1))],...
       ['y: ',num2str(pos(2))],['z: ',num2str(pos(3))]};
targetup=get(get(event_obj,'Target'),'parent');
set(targetup,'userdata',struct('pos',pos));
pt=get(gcf,'userdata');
if(size(pt,1)<5)
     set(gcf,'userdata',[pt;pos]);
end
