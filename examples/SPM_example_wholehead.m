clear
%% Dependencies brain2mesh, iso2mesh, nii
%% This uses the nii viewer at: 
%% www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

%% For brain2mesh, the WM is considered first. Thus the order of the tissues is
%% reversed with the GM.

check_brain2mesh_dependency;

tissues = [2 1 3 4 5];

%% The SPM files are read and stored into the 4D array 'seg'
for i = 1:5
    A = load_nii(sprintf('SPM/c%iANTS19-5Years_head.nii.gz',tissues(i)));
    dim = size(A.img);
    seg(:,:,:,i) = A.img;
end

%% The one liner brain2mesh is called. The third parameter controls activates the whole head
%% pathway
[node,elem,face] = brain2mesh(seg,1,1,[2 2 3 3.5 3.5]);


%% Plotting of the result
plotmesh(node,elem(elem(:,5)==1,:),'FaceColor',[1 1 1],'EdgeAlpha',0.6) %%wm
hold on;
plotmesh(node,elem(elem(:,5)==2,:),'x>78| y<110','FaceColor',[0.35 0.35 0.35],'EdgeAlpha',0.6) %%pial
plotmesh(node,elem(elem(:,5)==3,:),'x>78 | z<135','FaceColor',[0.2 0.6 1],'EdgeAlpha',0.6) %%csf
plotmesh(node,elem(elem(:,5)==4,:),'x>90 | z<130','FaceColor',[1 1 0.9],'EdgeAlpha',0.6) %%bone
plotmesh(node,elem(elem(:,5)==5,:),'x>100 | z<125','FaceColor',[1 0.8 0.7],'EdgeAlpha',0.6) %%scalp
view([-0.6 1.5 0.6])
