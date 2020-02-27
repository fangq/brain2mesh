function [len, node]=polylinelen(node, p0, p1)
%
% [len, node]=polylinelen(node, p0, p1)
%
% Calculate the polyline line segment length vector in sequential order
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    node: an N x 3 array defining each vertex of the polyline in
%          sequential order
%    p0:(optional) a given node to define the start of the polyline, if not
%         defined, start position is assumed to be 1st node
%    p1:(optional) a given node to define the end of the polyline, if not
%         defined, start position is assumed to be last node
%
% output:
%    len: the length of each segment between the start and the end points
%    node: the node list between the start and end points of the polyline
%
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v2 or later, see LICENSE.txt for details
%

if(nargin<3)
    p1=size(node,1);
    if(nargin<2)
        p0=1;
    end
end

if(size(p0,2)==3)
    p0=closestpt(node,p0);
end
if(size(p1,2)==3)
    p1=closestpt(node,p1);
end

if(p0>p1)
    node=node([p0:end 1:p1],:);
else
    node=node(p0:p1,:);
end

dist= node(1:end-1,:) - node(2:end,:);
len=sum(sqrt(sum(dist.*dist,2)));
