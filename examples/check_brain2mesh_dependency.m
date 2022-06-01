if(~exist('brain2mesh','file'))
    error('Missing dependency. You must download and addpath to brain2mesh.m, URL: https://github.com/fangq/brain2mesh')
end

if(~exist('zmat','file'))
    error('Missing dependency. You must download and addpath to zmat.m, URL: https://github.com/fangq/zmat')
end

if(~exist('nii2jnii','file'))
    error('Missing dependency. You must download and addpath to iso2mesh toolbox, URL: https://github.com/fangq/iso2mesh')
end

disp('all dependencies are present. brain2mesh is ready to run');