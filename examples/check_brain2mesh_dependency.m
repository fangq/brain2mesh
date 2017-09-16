if(~exist('brain2mesh','file'))
    error('Missing dependency. You must download and addpath to brain2mesh.m, URL: https://github.com/fangq/brain2mesh')
end

if(~exist('load_nii','file'))
    error('Missing dependency. You must download and addpath to load_nii.m, URL: https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image')
end

if(~exist('v2m','file'))
    error('Missing dependency. You must download and addpath to iso2mesh toolbox, URL: https://github.com/fangq/iso2mesh')
end
