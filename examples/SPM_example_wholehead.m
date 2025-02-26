clear;
%% Dependencies: brain2mesh, iso2mesh, and zmat (https://github.com/NeuroJSON/zmat)

%% Reading in data from c1-c5 SPM segmentations

% the SPM segmentation files are stored in the text-based JNIfTI format
% (.jnii) as defined in the JNIfTI file specification
% https://github.com/NeuroJSON/jnifti/

names = {'gm', 'wm', 'csf', 'skull', 'scalp', 'air'};
for i = 1:5
    A = loadjnifti(sprintf('jnii/c%iANTS40-44Years_head.jnii', i));
    dim = size(A.NIFTIData);
    seg.(names{i}) = A.NIFTIData;
end

%% call brain2mesh to create multi-layered brain mesh; dotruncate cuts the mesh below brain

cfg.smooth = 10;
tic;
[node, elem, face] = brain2mesh(seg, cfg);
toc;

%% call brain1020 to create the 10-5 landmarks on the scalp
initpoints = [
              86.4888  191.1470  100.8055
              82.8704    1.2961  114.1993
              15.1882   93.7385   71.0191
              158.4230   90.3180   77.2665
              83.7306  102.2434  207.2162];

headsurf = volface(elem(:, 1:4));
tic;
[landmarks, curves] = brain1020(node, headsurf, initpoints, 10, 10, 'cztol', 1e-8);
toc;
view([-0.6 1.5 0.6]);

%% Plotting of the result
figure;
plotmesh(node, elem(elem(:, 5) == 5, :), 'FaceColor', [1 1 1], 'EdgeAlpha', 0.6); %% wm
hold on;
plotmesh(node, elem(elem(:, 5) == 4, :), 'x>78 | y<110', 'FaceColor', [0.35 0.35 0.35], 'EdgeAlpha', 0.6); %% pial
plotmesh(node, elem(elem(:, 5) == 3, :), 'x>78 | z<135', 'FaceColor', [0.2 0.6 1], 'EdgeAlpha', 0.6); %% csf
plotmesh(node, elem(elem(:, 5) == 2, :), 'x>90 | z<130', 'FaceColor', [1 1 0.9], 'EdgeAlpha', 0.6); %% bone
plotmesh(node, elem(elem(:, 5) == 1, :), 'x>100 | z<125', 'FaceColor', [1 0.8 0.7], 'EdgeAlpha', 0.6); %% scalp
view([-0.6 1.5 0.6]);
