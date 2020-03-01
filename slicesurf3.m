function [leftpt,leftcurve,rightpt, rightcurve]=slicesurf3(node,elem,p1,p2,pmid, step)

fullcurve=slicesurf(node, elem, [p1;pmid;p2]);
[fulllen, fullcurve]=polylinelen(fullcurve, p1,pmid,p2);

[leftlen,  leftcurve]=polylinelen(fullcurve, p2, p1);
[idx, weight, leftpt]=polylineinterp(leftlen, sum(leftlen)*[step:step:(100-step*0.5)]*0.01, leftcurve);
if(nargout>2)
    [rightlen, rightcurve]=polylinelen(fullcurve, p2, pmid);
    [idx, weight, rightpt]=polylineinterp(rightlen, sum(rightlen)*[step:step:(100-step*0.5)]*0.01, rightcurve);
end
