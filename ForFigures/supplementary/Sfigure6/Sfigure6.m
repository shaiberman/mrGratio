clear   , close all

%% mean FVF and MTV value along corpus callosum (mean over subjects)

currentPath=fileparts(which(mfilename));

% load data
load(fullfile(currentPath,'Sfigure6_fa_mtv.mat'))
% load cleaned segmentation
maskPath=fullfile(currentPath,'fiberGroups_fig6/CCseg_highRes_sub2.nii.gz');
CCseg=readFileNifti(maskPath);

figure,  colormap hot,

% Make the CC segmentation, along one slice
% Each voxel in the ROI is [1,8], according to the functional subregions
CC=rot90(squeeze(permute(CCseg.data,[2,3,1])));
slice=100;
CCslice=squeeze(CC(:,:,slice));

% Initialize variables
areaMap_tv=zeros(size(CCslice));
areaMap_fa=zeros(size(CCslice));

% Color coded according to mean g value colors.

for kk=1:8;
    Ind=find(CCslice==kk);
    % Each voxel which belongs to a subregion is assigned the mean FA/MTV value
    % within the subregion.
    areaMap_tv(Ind)=meanTV(kk);
    areaMap_fa(Ind)=meanFA(kk);
end

% plot MTV
subplot(2,1,1),imagesc(areaMap_tv(50:130,50:200)),
colorbar,caxis([0.275 0.295])
title('MTV')
set(gca,'xTick',[],'yTick',[])

% plot FA
subplot(2,1,2),imagesc(areaMap_fa(50:130,50:200)),
colorbar,caxis([0.73    0.83])
title('FA')
set(gca,'xTick',[],'yTick',[])