
clear
close all
currentPath=fileparts(which(mfilename));
%% choose colors
colors = colormap('lines');
colors  = [colors(1:7,:);0,0.4,0.2];
colors2=colors;
colors1=repmat([0.45 0.45 0.45],8,1);

%% load the data with the fits
load(fullfile(currentPath,'fitting_gRatio_allAges_VoxelWise_gender.mat'))

%% plot gratio
fgNum=1:8;
valName='gr';
fgNames={'Occipital',  'Temporal',  'Post-Parietal',   'Sup-Parietal',  'Motor',  'Sup-Frontal',  'Ant-Frontal',  'Orb-Frontal'};

f = nc_PlotModelFits_gender(coefsM(:),coefsF(:),valName,fgNames,fgNum,colors1,colors2,RsM,RsF);

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


%% create legend
figure, hold on
plot([1 2],[1 1],'LineWidth',20,'Color',[0.45 0.45 0.45])
xSegs=linspace(1,2,9);
for ii=1:8
    plot([xSegs(ii), xSegs(ii)+(xSegs(2)-xSegs(1))],[0.8 0.8],'LineWidth',20,'Color',colors2(ii,:))
end
% because of the line width, the segments over lap. remove the extra bit by adding a white section
plot([xSegs(ii+1), xSegs(ii+1)+(xSegs(2)-xSegs(1))],[0.8 0.8],'LineWidth',20,'Color',[1 1 1])
plot([xSegs(ii+1), xSegs(ii+1)+(xSegs(2)-xSegs(1))],[1 1],'LineWidth',20,'Color',[1 1 1])
axis([0 4 0 2])
text(1.95,1,'Males','fontSize',16,'fontName','times')
text(1.95,0.8,'Females','fontSize',16,'fontName','times')



 