function [vol, names]=tpm2label(seg, segorder)
%
% vol=tpm2label(seg)
%   or
% vol=tpm2label(seg, segorder)
% [vol, names]=tpm2label(seg, segorder)
%
% converting tissue probablistic maps (TPMs) to a multi-lable volume
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    seg: a struct, a cell array, or a 3D or 4D array; if seg is a cell
%         or a struct, their subfields must be 2D/3D arrays of the same
%         sizes;
%    segorder: if seg is a struct, segorder allows one to assign output
%         labels using customized order instead of the creation order
%
% output:
%    vol: a 2-D or 3-D array of the same type/size of the input arrays. The
%         label for each voxel is determined by the index to the highest
%         value in TPM of the same voxel. If a voxel is a background voxel
%         - i.e. zeros for all TPMs, it stays 0
%    names: a cell array storing the names of the labels (if input is a
%         struct), the first string is the name for label 1, and so on
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

mask=seg;
names={};
if(isstruct(seg))
    if(nargin>1)
        seg=orderfields(seg,segorder);
    end
    names=fieldnames(seg);
    mask=cellfun(@(x) seg.(x), names,'UniformOutput',false);
end
if(iscell(mask))
    mask=cat(ndims(mask{1})+1,mask{:});
end
if(~isnumeric(mask))
    error('input must be a cell/struct array with numeric elements of matching dimensions');
end

[newmask, vol]=max(mask,[],ndims(mask));
vol=vol .* (sum(mask,ndims(mask))>0);