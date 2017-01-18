clear
close all
currentPath=fileparts(which(mfilename));
%% choose colors
colors = colormap('lines');
colors  = [colors(1:7,:);0,0.4,0.2];

%% plot gratio
valName='gr';
fgNum=1:8;
models={'linear'};
fgNames={'Occipital',  'Temporal',  'Post-Parietal',   'Sup-Parietal',  'Motor',  'Sup-Frontal',  'Ant-Frontal',  'Orb-Frontal'};
fitPath=fullfile(currentPath,'fitting_gRatio_allAges_VoxelWise_weighted.mat');
load(fitPath)
for j=1:length(models)
    coefstemp=coefs(:,j); Rstemp=Rs(:,j);
    f = nc_PlotModelFitsJ_fig2(coefstemp,valName,fgNames,fgNum,colors,Rstemp);
end
%% add the segmentation 

load(fullfile(currentPath, 'CCseg.mat'))

% load the weights, a T1 map and  brainmask
load(fullfile(currentPath,'T1Slice.mat'));

% create a slice of T1 midSagittal map
gVals=repmat(double((T1Slice==0)),1,1,3);

 for kk=1:8;
    c=colors(kk,:);
    [Indx,Indy]=find(areaMap==kk);
    for jj=1:length(Indx)
        gVals(Indx(jj),Indy(jj),:)=[c(1),c(2),c(3)];
    end
end

gMask=mean(gVals,3);gMask=logical(gMask);gMask(gMask==1)=0.3;
g = subplot(3,3,9);hold on
imshow(T1Slice),caxis([0 5])
h=imshow(gVals);
set(h,'AlphaData',gMask)
p=get(g,'position');
set(g,'Position',[p(1)-0.02,p(2)-0.07,p(3)+0.07,p(4)+0.07])

 