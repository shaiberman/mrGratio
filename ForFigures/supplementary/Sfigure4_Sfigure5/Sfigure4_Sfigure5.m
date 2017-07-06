clear   , close all

%% plot MTV and FVF as a functionb of age and Sex

currentPath=fileparts(which(mfilename));

% choose colors
colors = colormap('lines');close(gcf)
colors  = [colors(1:7,:);0,0.4,0.2];
colors2=colors;
colors1=repmat([0.45 0.45 0.45],8,1);

% define CC subregion
fgNames={'Occipital',  'Temporal',  'Post-Parietal',   'Sup-Parietal',  'Motor',  'Sup-Frontal',  'Ant-Frontal',  'Orb-Frontal'};
fgNum=1:8;
%% MTV
% load MTV data
models={'quadratic'}; 
load('fitting_TV_allAges_VoxelWise_gender_median.mat');

% plot MTV
valName='MTV';
f = nc_PlotModelFits_gender(coefsM,coefsF,valName,fgNames,fgNum,colors1,colors2,RsM,RsF);

%% FVF
% load FVF data fits
clear RsM ErrsM ghatsM coefsM RsF ErrsF ghatsF coefsF
models={'linear','quadratic','piecewise','exponent','lowess','poisson'}; % left out "piecewise", and exponent analyses
load('fitting_FVF_allAges_VoxelWise_gender_median.mat')

% plot FVF
valName='FVF'; j=4; % I tried different models and plotting the exoponential fit.
coefstempM=coefsM(:,j);     coefstempF=coefsF(:,j);
RstempM=RsM(:,j);    RstempF=RsF(:,j);
f = nc_PlotModelFits_gender(coefstempM,coefstempF,valName,fgNames,fgNum,colors1,colors2,RstempM,RstempF);



%%  load the segmentation file;
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
for ii=1:2
    figure(ii)
    g = subplot(3,3,9);hold on
    imshow(T1Slice),caxis([0 5])
    h=imshow(gVals);
    set(h,'AlphaData',gMask)
    p=get(g,'position');
    set(g,'Position',[p(1)-0.02,p(2)-0.07,p(3)+0.07,p(4)+0.07])
end