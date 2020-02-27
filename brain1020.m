function [landmarks, initpoints]=brain1020(node, elem, initpoints, perc1, perc2, varargin)
%
% landmarks=brain1020(node, elem)
%   or
% landmarks=brain1020(node, elem, initpoints)
% landmarks=brain1020(node, elem, initpoints, 10, 5)
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
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v2 or later, see LICENSE.txt for details
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

if(isstruct(initpoints))
    initpoints=orderfields(initpoints,{'nz', 'iz', 'lpa', 'rpa', 'cz'});
    landmarks=initpoints;
    initpoints=struct2array(initpoints);
    initpoints=reshape(initpoints(:),3,length(initpoints(:))/3)';
else
    landmarks=struct('nz', initpoints(1,:),'iz', initpoints(2,:),...
                     'lpa',initpoints(3,:),'rpa',initpoints(4,:),...
		     'cz', initpoints(5,:));
end

if(size(elem,2)>=4)
    elem=volface(elem(:,1:4));
end

if(isempty(initpoints) && size(initpoints,1)<5)
    hf=figure;
    plotmesh(node,elem);
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

initpoints

% at this point, initpoints contains {nz, iz, lpa, rpa, cz0}

%% Step 1: nz, iz and cz0 to determine saggital reference curve
[nsagg, isagg]=slicehead(node, elem, initpoints([1,2,5],:));
isagg(isnan(isagg))=[];

nsagg=nsagg(isagg,:);
ncoro=ncoro(icoro,:);

%% Step 2: get true cz as the mid-point between iz and nz
[len, nsagg]=polylinelen(nsagg, initpoints(1,:), initpoints(2,:));
[cz, index, weight]=polylineinterp(nsagg, len/2);
initpoints(5,:)=cz;
landmarks.cz=cz;

%% Step 3: lpa, rpa and cz to determine coronal reference curve
[ncoro, icoro]=slicehead(node, elem, initpoints([3,4,5],:));
icoro(isnan(icoro))=[];

%% debug
plotmesh(node,elem);
hold on;
plot3(nsagg(isagg,1),nsagg(isagg,2),nsagg(isagg,3),'r','LineWidth',4);
plot3(ncoro(icoro,1),ncoro(icoro,2),ncoro(icoro,3),'g','LineWidth',4);

%% helper functions
%---------------------------------------------------------------------------
% the respond function when there is a data-tip to popup
%---------------------------------------------------------------------------
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


