clear
%% Dependencies brain2mesh, iso2mesh, nii
%% This uses the nii viewer at: 
%% www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

%% Reading in data from c1-c5 SPM segmentations using a 4D array
for i = 1:5
    A = load_nii(sprintf('c%iANTS19-5Years_head.nii.gz',i));
    dim = size(A.img);
    data(:,:,:,i) = A.img;
end

%% Loading using a structure
seg.wm = data(:,:,:,2); %loading the white matter with the corresponding SPM file
seg.gm = data(:,:,:,1);
seg.csf = data(:,:,:,3); 
seg.skull = data(:,:,:,4); 
seg.scalp = data(:,:,:,5); 

%% Alternative loading option using a 4D array.
% tissue_order = [5 4 3 1 2]; %Reordering tissues with scalp first going inward to WM
% for i = 1:5
%    seg(:,:,:,i) = data(:,:,:,tissue_order(i));
% end

%% The one liner brain2mesh is called. The third parameter controls activates the whole head
%% pathway
[node,elem,face] = brain2mesh(seg,1,0);

%% Plotting of the result
plotmesh(node,elem(elem(:,5)==5,:),'FaceColor',[1 1 1],'EdgeAlpha',0.6) %%wm
hold on;
plotmesh(node,elem(elem(:,5)==4,:),'x>78| y<110','FaceColor',[0.35 0.35 0.35],'EdgeAlpha',0.6) %%pial
plotmesh(node,elem(elem(:,5)==3,:),'x>78 | z<135','FaceColor',[0.2 0.6 1],'EdgeAlpha',0.6) %%csf
plotmesh(node,elem(elem(:,5)==2,:),'x>90 | z<130','FaceColor',[1 1 0.9],'EdgeAlpha',0.6) %%bone
plotmesh(node,elem(elem(:,5)==1,:),'x>100 | z<125','FaceColor',[1 0.8 0.7],'EdgeAlpha',0.6) %%scalp
view([-0.6 1.5 0.6])