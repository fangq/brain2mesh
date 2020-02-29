function seg=label2tpm(vol, names)
%
% seg=label2tpm(vol)
%   or
% seg=label2tpm(vol, names)
%
% converting a multi-lable volume to binary tissue probablistic maps (TPMs) 
%
% author: Qianqian Fang (q.fang at neu.edu)
%
% input:
%    vol: a 2-D or 3-D array of integer values.
%    names (optional): a cell array of strings defining the names of each label,
%         alternatively, a containers.Map object with label (integer) as
%         the key and name as the value. The names can be a subset of all
%         labels.
%
% output:
%    seg: a struct with subfields of 3D or 4D uint8 array; the subfield 
%         names are defined via the optional names input; the unnamed
%         labels will be named as 'label_#'. 
%
% -- this function is part of brain2mesh toolbox (http://mcx.space/brain2mesh)
%    License: GPL v3 or later, see LICENSE.txt for details
%

if(~isnumeric(vol))
    error('input must be a numerical array');
end

val=setdiff(sort(unique(vol(:))),0);
if(length(val)>numel(vol)*(0.1^ndims(vol)))
    error('please convert the input to labels first');
end

seg=struct;
for i=1:length(val)
    nm=sprintf('label_%d',val(i));
    if(nargin>1)
        if(iscell(names))
            if(i<=length(names))
                nm=names{i};
            end
        elseif(isa(names,'containers.Map') && isKey(names,val(i)))
            nm=names(val(i));
        end
    end
    seg.(nm)=uint8(vol==val(i));
end
