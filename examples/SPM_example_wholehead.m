clear
%% Dependencies: brain2mesh, iso2mesh, and zmat (http://github.com/fangq/zmat)

%% Reading in data from c1-c5 SPM segmentations

% the SPM segmentation files are stored in the text-based JNIfTI format
% (.jnii) as defined in the JNIfTI file specification
% https://github.com/fangq/jnifti/ 

names={'gm','wm','csf','skull','scalp','air'};
for i = 1:5
    A = loadjnifti(sprintf('jnii/c%iANTS40-44Years_head.jnii',i));
    dim = size(A.NIFTIData);
    seg.(names{i}) = A.NIFTIData;
end

%% Alternative loading option using a 4D array.
% tissue_order = [5 4 3 1 2]; %Reordering tissues with scalp first going inward to WM
% for i = 1:5
%    seg(:,:,:,i) = data(:,:,:,tissue_order(i));
% end


%% The one liner brain2mesh is called. The third parameter controls activates the whole head
%% pathway
cfg.dotruncate = 1;
[node,elem,face] = brain2mesh(seg,cfg);


%% Plotting of the result
plotmesh(node,elem(elem(:,5)==5,:),'FaceColor',[1 1 1],'EdgeAlpha',0.6) %%wm
hold on;
plotmesh(node,elem(elem(:,5)==4,:),'x>78| y<110','FaceColor',[0.35 0.35 0.35],'EdgeAlpha',0.6) %%pial
plotmesh(node,elem(elem(:,5)==3,:),'x>78 | z<135','FaceColor',[0.2 0.6 1],'EdgeAlpha',0.6) %%csf
plotmesh(node,elem(elem(:,5)==2,:),'x>90 | z<130','FaceColor',[1 1 0.9],'EdgeAlpha',0.6) %%bone
plotmesh(node,elem(elem(:,5)==1,:),'x>100 | z<125','FaceColor',[1 0.8 0.7],'EdgeAlpha',0.6) %%scalp
view([-0.6 1.5 0.6])