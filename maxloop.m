function newloop = maxloop(curveloop)
%
% newloop = maxloop(curveloop)
%
% return the curve segment that has the largest number of nodes
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    curveloop: curves defined by node indices, separated by nan
%
% output:
%    newloop: the node indices defining the longest segment
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

newloop = curveloop;

if (~isnan(curveloop(end)))
    curveloop(end + 1) = nan;
end

loopend = find(isnan(curveloop));
if (length(loopend) > 1)
    seglen = [loopend(1), diff(loopend)];
    [maxlen, maxloc] = max(seglen);
    loopend = [0 loopend];
    newloop = curveloop((loopend(maxloc) + 1):(loopend(maxloc + 1) - maxloc));
end
newloop(isnan(newloop)) = [];
