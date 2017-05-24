 clear,  close all,
figure(1),
currentPath=fileparts(which(mfilename));
%% choose colors
colors = colormap('lines');
%%  1.1: plot West et al.'s data

% load data:
load(fullfile(currentPath,'figure1_West.mat'))

% plot:
gr={'Control','Rictor CKO','Tsc2 CKO','Pten CKO'};
subplot(1,2,1), hold on
cl  = [colors([1,2,4],:);0,0.4,0.2];

for kk=1:length(gr)
    
    %   plot the error bars
    errorbar(MFhist.mean(kk),   TV.mean(kk),TV.std(kk),'.','Color',[0.68 ,.68 ,.68]);
    temp=herrorbar(MFhist.mean(kk),   TV.mean(kk), MFhist.std(kk),'.');
    
    for j = 1:length(temp)
        set(temp(j),'color',[0.68 ,.68 ,.68]);
    end
end

% replot the data point to make them visible
for ii=1:length(gr)
    h(ii)= scatter(MFhist.mean(ii),TV.mean(ii),70,cl(ii,:),'filled');
end

% add the identity line
h(5)=identityLine(gca);

% add a regression line
x0=0:0.05:0.35;
p=polyfit(MFhist.mean(:),TV.mean(:),1);
y0=p(1)*x0+p(2);
h(6)=plot(x0,y0,'-k');

% calculate the coefficient of determination
r=calccod(MFhist.mean,TV.mean);  

% add information and design
set(gca,'FontWeight','normal','fontSize',14)
title(['R^2 = ',num2str(round(r,2))], 'FontSize',18,'FontWeight','normal')
xlabel('MVF_h_i_s_t','FontSize',18),  ylabel('MTV_e_x_-_v_i_v_o','FontSize',18)
set(gca,'xTick',[0:0.05:0.3],'yTick',[0:0.05:0.3]),
axis equal, axis([0 0.35 0 0.35])
grid on, box on
set(gca,'xTickLabel',{'0',' ','0.1',' ','0.2',' ','0.3'},'yTickLabel',{'0',' ','0.1',' ','0.2',' ','0.3'})
set(gca,'FontWeight','normal')
legend(h,{gr{:},'identity line','linear fit'},'location','northwest', 'EdgeColor','k', 'location','northwest', 'FontSize',13)

%%  2: plot the reproducibility of 4  subjects

% load the data of these subjects
load(fullfile(currentPath,'figure1_reprodData.mat'))

% this structure contains mean and std of the g ratio values n 4 subjects'
% 2 separate scan. from each scan we have the data within each of the 8
% subregions of the CC. te 2 scans for each subject are consecutive, such
% that gSub.mean(2,4) is the mean g ratio of the 4thj CC subregion, iun the
% secoind scan of the 1st subject.

cl  =  [colors([3,5:7],:)];
subplot(1,2,2), hold on

l=length(gSub.mean);k=0;
for ii=1:2:l
    k=k+1;
    %   plot the error bars
    errorbar(gSub.mean(ii,:),gSub.mean(ii+1,:),gSub.std(ii,:),'.','color',[0.68 ,.68 ,.68]);
    temp = herrorbar(gSub.mean(ii,:),gSub.mean(ii+1,:),gSub.std(ii+1,:),'.');
    for j = 1:length(temp)
        set(temp(j),'color',[0.68 ,.68 ,.68]);
    end
end

% replot the data point to make them visible
k=0;
for ii=1:2:l
    k=k+1;
    h(k)= scatter(gSub.mean(ii,:),gSub.mean(ii+1,:),80,cl(k,:),'filled');
end

% add the identity line
h(5)=identityLine(gca);


% calculate the coefficient if determination
x=gSub.mean([1:2:l],:); y=gSub.mean([2:2:l],:);
r=calccod(x(:),y(:));

% add a linear regression
x0=0.5:0.1:0.8;
p=polyfit(x(:),y(:),1);
y0=p(1)*x0+p(2);
h(6)=plot(x0,y0,'-k');

% add information and design
set(gca,'FontWeight','normal','fontSize',14)
title(['R^2=',num2str(round(r,2))], 'FontSize',18,'FontWeight','normal'),
xlabel('GR*: scan 1','FontSize',18), ylabel('GR*: scan 2','FontSize',18)
axis equal, axis([0.5 0.8 0.5 0.8])
set(gca,'xTick',[0.5:0.05:0.8],'yTick',[0.5:0.05:0.8]),
set(gca,'xTickLabel',{'0.5',' ','0.6',' ','0.7',' ','0.8'},'yTickLabel',{'0.5',' ','0.6',' ','0.7',' ','0.8'})
box on, grid on
legend(h,{'subject 1','subject 2','subject 3','subject 4','identitity line','linear fit'}, 'EdgeColor','k', 'location','northwest', 'FontSize',13)


%% individual correlation coefficient for the different subjects
for ii=1:4
    tmp=corrcoef(x(ii,:),y(ii,:));
    rs(ii)=tmp(2);
end
%  0.9522    0.9681    0.9806    0.9868
meanR=mean(rs);


