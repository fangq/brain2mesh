function [leftpt,leftcurve,rightpt, rightcurve]=slicebetween(node,elem,p1,p2,p3, step)

fullcurve=slicehead(node, elem, [p1;p2;p3]);
[fulllen, fullcurve]=polylinelen(fullcurve, p1,p3,p2);

[leftlen,  leftcurve]=polylinelen(fullcurve, p2, p1);
[idx, weight, leftpt]=polylineinterp(leftlen, sum(leftlen)*[step:step:(100-step*0.5)]*0.01, leftcurve);
if(nargout>2)
    [rightlen, rightcurve]=polylinelen(fullcurve, p2, p3);
    [idx, weight, rightpt]=polylineinterp(rightlen, sum(rightlen)*[step:step:(100-step*0.5)]*0.01, rightcurve);
end