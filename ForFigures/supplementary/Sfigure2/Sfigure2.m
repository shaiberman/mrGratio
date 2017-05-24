close all, clear
%% Bland-Altman plot of the same data 

figure(1),
currentPath=fileparts(which(mfilename));
colors = colormap('lines');

%% MTV-MRI vs histological MVF

% load data
load(fullfile(currentPath,'figureS2_West.mat'))

% plot
gr={'Control','Rictor CKO','Tsc2 CKO','Pten CKO'};
cl  = [colors([1,2,4],:);0,0.4,0.2];
subplot(2,1,1), hold on

% calculate the necessary parameters for the plot
diff = TV-MVFhist;
mn   = (TV+MVFhist)/2;
sd = nanstd(diff(:));
mndiff = nanmean(diff(:));

for kk = 1:size(TV,2);
            
    x=mn(:,kk);
    y=diff(:,kk);
    % plot the difference between each point pair, as a function of their mean
    h(kk)= scatter(x,y,70,cl(kk,:),'filled');
end   
hold on
plot([min(mn(:)), max(mn(:))], [mndiff mndiff],'k','linewidth',2);
plot([min(mn(:)), max(mn(:))], [mndiff+1.96*sd mndiff+1.96*sd],'--','color',[0.68 ,.68 ,.68],'linewidth',2);
plot([min(mn(:)), max(mn(:))], [mndiff-1.96*sd mndiff-1.96*sd],'--','color',[0.68 ,.68 ,.68],'linewidth',2);
%  title([ 'Mean Difference =' sprintf('%.4f',mndiff), '/n STD =' sprintf('%.4f',sd) ]);
 
xlabel('mean MVF_h_i_s_t  &  MTV_e_x_-_v_i_v_o','FontSize',13)
ylabel('MVF_h_i_s_t  -  MTV_e_x_-_v_i_v_o','FontSize',13)
legend(h,gr,'location','southwest', 'EdgeColor','k', 'location','northwest', 'FontSize',12)


%% Reproducibility of four subjects
% load data
load(fullfile(currentPath,'figureS2_reprodData.mat'))

% plot
subplot(2,1,2), hold on
cl  =  [colors([3,5:7],:)];
l=length(gSub.mean);k=0;
x=gSub.mean([1:2:l],:);y=gSub.mean([2:2:l],:);

% calculate the necessary parameters for the plot
diff = x-y;
mn   = (x+y)/2;
sd = nanstd(diff(:));
mndiff = nanmean(diff(:));

for kk = 1:size(x,1);
            
    m=mn(kk,:);
    d=diff(kk,:);
     % plot the difference between each point pair, as a function of their
     % mean
     h(kk)= scatter(m,d,70,cl(kk,:),'filled');
end   
hold on
plot([min(mn(:)), max(mn(:))], [mndiff mndiff],'k','linewidth',2);
plot([min(mn(:)), max(mn(:))], [mndiff+1.96*sd mndiff+1.96*sd],'--','color',[0.68 ,.68 ,.68],'linewidth',2);
plot([min(mn(:)), max(mn(:))], [mndiff-1.96*sd mndiff-1.96*sd],'--','color',[0.68 ,.68 ,.68],'linewidth',2);
% title([ 'Mean Difference =' sprintf('%.4f',mndiff), '/n STD =' sprintf('%.4f',sd) ]);
xlabel('mean scan1  &  scan2','FontSize',13)
ylabel('scan1  -  scan2','FontSize',13)
legend(h,gr,'location','southwest', 'EdgeColor','k', 'location','northwest', 'FontSize',12)
