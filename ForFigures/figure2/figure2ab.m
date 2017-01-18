clear
close all,

currentPath=fileparts(which(mfilename));

%% 2a.plot the histogram  of all voxels - in all subjects
load(fullfile(currentPath,'figure2a_histData.mat'))
data = gTotal; 

exgMin=prctile(data,1);
exgMax=prctile(data,100);
median_g=prctile(data(data>exgMin & data<exgMax),50);
mad_g=mad(data(data>exgMin & data<exgMax),1);

figure,hold on,
histogram(data(data>exgMin & data<exgMax),200)
plot([median_g,median_g,],[5, 880],'r','linewidth',4)

xlabel('g ratio'), ylabel('frequency')

%% 2b. sub-region wise g-ratio

% load cleaned segmentation
maskPath=fullfile(currentPath,'fiberGroups_fig2/CCseg_highRes_sub2.nii.gz');
CCseg=readFileNifti(maskPath);
CC=rot90(squeeze(permute(CCseg.data,[2,3,1])));

% load the weights, a T1 map and  brainmask
R1=readFileNifti(fullfile(currentPath,'fiberGroups_fig2/R1_map_2DTI_ForFig.nii.gz'));
% create a slice of T1 midSagittal map
slice=100;
R1=rot90(squeeze(permute(R1.data,[2,3,1])));
R1Slice=squeeze(R1(:,:,slice));
CCslice=squeeze(CC(:,:,slice));

% make the CC segmentation, color coded according to mean g value colors.
load(fullfile(currentPath,'figure2_SubData.mat'))
meanG=nanmean(gSub.mean);

gColors={[1 1 0.4353],[1 0.9569, 0],[1 1 0.2471],[1 0.7059 0] ,[1 0 0],...
    [0.5412, 0 0],[0.9569 0 0],[0.6667 0 0]};
% these are the RGB values that correspond to the mean g values with a hot
% color map and the caxis([0.63 0.72])
% gMap=zeros(size(T1Slice));
gVals=zeros([size(R1Slice),3]);
areaMap=zeros(size(R1Slice));
for kk=1:8;
    
    c=gColors{kk};
    
    [Indx,Indy]=find(CCslice==kk);
    Ind=find(CCslice==kk);
    areaMap(Ind)=meanG(kk);
    
    for jj=1:length(Indx)
        gVals(Indx(jj),Indy(jj),:)=[c(1),c(2),c(3)];
    end
    
end

% create a transparent map to overlay CC and T1 map
gMask=mean(gVals,3);gMask=logical(gMask);gMask(gMask==1)=0.3;

%  project the segmentation on the T1 map
figure, imshow(R1Slice),caxis([0.2 1]),hold on
h=imshow(gVals);
set(h,'AlphaData',gMask)

 figure, imagesc(areaMap), 
 colormap hot, 
 caxis([0.63 0.72]),
 colorbar
