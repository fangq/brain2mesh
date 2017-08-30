function I=max_filter(I,n,f)
% I input 3D image
% n neighbourhood size. for example 3 for 3*3
% f=1 max filter
% f=0 min filter
if f==1 l=n; else l=1; end
m=size(I);
for ii=1:m(3)
I(:,:,ii) = ordfilt2(I(:,:,ii),l*l,ones(n,n));
end
 for ii=1:m(1)
I(ii,:,:) = ordfilt2(squeeze(I(ii,:,:)),l,ones(1,n));
end