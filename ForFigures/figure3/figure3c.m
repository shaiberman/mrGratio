
clear
close all,

currentPath=fileparts(which(mfilename));
%% choose colors
cl=[0.1,0.1,0.9 ; 0.1,0.4,0.4 ; 0.9,0.4,0.4];
%% bar plot comparing our data vs Macaque data in the CC

% Macaque data: 
% the following taken from  a table published by: 
% Stikov, N., Campbell, J. S., Stroh, T., Lavelée, M., Frey, S., Novek, J.,
% ... & Leppert, I. R. (2015). In vivo histology of the myelin g-ratio with
% magnetic resonance imaging. NeuroImage, 118, 397-405.‏
meanG3=[0.74, 0.6, 0.64, 0.72, 0.69, 0.64, 0.68, 0.69]; 
errG3=[0.06, 0.04, 0.01, 0.03, 0.04, 0.06, 0.04, 0.05];

% load our human data
load(fullfile(currentPath, 'figure2_SubData.mat'))
meanG1=nanmean(gSub.mean);
errG1=nanmean(gSub.std);

% take only the areas that match in terms of segmentation
Splenium = [meanG1(1),  meanG3(1)];
Body     = [ meanG1(5), meanG3(5)];
Genu     = [meanG1(8),  meanG3(8)];
bar_input=[Splenium; Body;Genu];

SpleniumEr= [errG1(1),errG3(1)];
BodyEr    = [ errG1(5),errG3(5)];  
GenuEr    = [errG1(8),errG3(8)];
errorbar_input=[SpleniumEr; BodyEr;GenuEr];

figure,hold on,
x=[0.8 1.2,1.8, 2.2, 2.8, 3.2];

% Splenium:
bar(x(1),meanG1(1),'lineWidth',2.5,'EdgeColor',cl(1,:),'FaceColor',[1 1 1],'barWidth',0.4);
bar(x(2), meanG3(1),'lineWidth',2.5,'EdgeColor',cl(1,:),'FaceColor',[1 1 1],'barWidth',0.4);

% Body
bar(x(3),meanG1(5),'lineWidth',2.5,'EdgeColor',cl(2,:),'FaceColor',[1 1 1],'barWidth',0.4);
bar(x(4),meanG3(5),'lineWidth',2.5,'EdgeColor',cl(2,:),'FaceColor',[1 1 1],'barWidth',0.4);

% Genu
bar(x(5),meanG1(8),'lineWidth',2.5,'EdgeColor',cl(3,:),'FaceColor',[1 1 1],'barWidth',0.4);
bar(x(6),meanG3(8),'lineWidth',2.5,'EdgeColor',cl(3,:),'FaceColor',[1 1 1],'barWidth',0.4);

set(gca,'xTick',[1 2 3],'xTickLabel',{'Splenium','midBody','Genu'},'FontSize',15)

errorbar(x(1:2), [meanG1(1),meanG3(1)], [errG1(1)/10,0],'lineWidth',1.8,'Color',cl(1,:),'Marker','.','lineStyle','none'),
errorbar(x(3:4),[meanG1(5),meanG3(5)],  [errG1(5)/10,0],'lineWidth',1.8,'Color',cl(2,:),'Marker','.','lineStyle','none')
errorbar(x(5:6),[meanG1(8),meanG3(8)], [errG1(8)/10,0],'lineWidth',1.8,'Color',cl(3,:),'Marker','.','lineStyle','none')

ylim([0.55 0.8]), xlim([0.3 3.6])
set(gca,'yTick',[0.55:0.1:0.75],'fontSize',15)
ylabel(' g-ratio ','fontSize',15)

%% visualize the fibers whose CC regions we are comparing 
load(fullfile(currentPath,'fiberGroups_fig2/CC_Occipital_clean_D5_L4.mat')), fg_a{1}=fg;
load(fullfile(currentPath,'fiberGroups_fig2/CC_Motor_clean_D5_L4.mat')),     fg_a{2}=fg;
load(fullfile(currentPath,'fiberGroups_fig2/CC_Orb_Frontal_clean_D5_L4.mat')),     fg_a{3}=fg;

AFQ_RenderFibers(fg_a{1},'numfibers',500,'color',cl(1,:),'camera','rightsag'); % splenium:
AFQ_RenderFibers(fg_a{2},'numfibers',500,'color',cl(2,:),'newfig',false,'camera','rightsag');% body:
AFQ_RenderFibers(fg_a{3},'numfibers',500,'color',cl(3,:),'newfig',false,'camera','rightsag');% genu            

R1file=fullfile(currentPath,'fiberGroups_fig2/R1_map_2DTI_ForFig.nii.gz');
R1=readFileNifti(R1file);
AFQ_AddImageTo3dPlot(R1,[-3,0,0]);
delete(findall(gcf,'type','light')); h = light; h.Position = [3 -50 20];
h2 = light; h2.Position = [3 -5 -20];

set(gca,'xtick',[],'ytick',[], 'ztick',[],'xlabel',[],'ylabel',[],'zlabel',[])
