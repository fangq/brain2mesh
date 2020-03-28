clear
%% Dependencies: brain2mesh, iso2mesh, and zmat (http://github.com/fangq/zmat)

%% Reading in data from c1-c5 SPM segmentations using a 4D array

% the SPM segmentation files are stored in the text-based JNIfTI format
% (.jnii) as defined in the JNIfTI file specification
% https://github.com/fangq/jnifti/ 

names={'gm','wm','csf','skull','scalp','air'};
for i = 1:5
    A = loadjnifti(sprintf('jnii/c%iANTS19-5Years_head.jnii',i));
    dim = size(A.NIFTIData);
    seg.(names{i}) = A.NIFTIData;
end

%% call brain2mesh to create multi-layered brain mesh; cfg.smooth performs 10-iter. surface smoothing 
cfg.dotruncate = 1;
%cfg.smooth=8; % apply smoothing to all surfaces, use a struct to smooth each
%cfg.radbound=struct('scalp',10,'skull',10,'csf',10,'gm',5,'wm',5); % set mesh density per layer
tic
[node,elem,face] = brain2mesh(seg, cfg);
toc

%% call brain1020 to create the 10-5 landmarks on the scalp
initpoints=[
   87.4123  188.5120   93.7087
   87.4360    5.9213  116.0523
   21.1737   93.9561   84.9446
  159.6440   89.8472   86.3139
   91.2203  104.7490  213.7747];

headsurf=volface(elem(:,1:4));
tic;
[landmarks, curves]=brain1020(node, headsurf, initpoints, 10,5,'cztol',1e-8, 'minangle', 0.75*pi);
toc
view([-0.6 1.5 0.6]);

%% Plotting of the result
plotmesh(node,elem(elem(:,5)==5,:),'FaceColor',[1 1 1]) %%wm
hold on;
plotmesh(node,elem(elem(:,5)==4,:),'x>78| y<110','FaceColor',[0.35 0.35 0.35]) %%pial
plotmesh(node,elem(elem(:,5)==3,:),'x>78 | z<135','FaceColor',[0.2 0.6 1]) %%csf
plotmesh(node,elem(elem(:,5)==2,:),'x>90 | z<130','FaceColor',[1 1 0.9]) %%bone
plotmesh(node,elem(elem(:,5)==1,:),'x>100 | z<125','FaceColor',[1 0.8 0.7]) %%scalp
view([-0.6 1.5 0.6])