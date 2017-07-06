close all, clear
%% 2D histogram of FA values in different diffusion weightings

currentPath=fileparts(which(mfilename));

% load data
load(fullfile(currentPath,'Sfigure3_FVF.mat'))

% plot
figure('units','normalized','position',[.1 .1 .7 .7])
plot2Dhist(fvf_highb,fvf_lowb,155,[0 1],[0 1],'FVF, b=2000','FVF, b=1000',0:0.2:1,0:0.2:1), hold on

% add a regression line
mdl=fitlm(fvf_highb,fvf_lowb,'linear');
p=mdl.Coefficients.Estimate(end:-1:1)';
x=[0,1]; y_fit=polyval(p,x);
plot(x,y_fit,'-','lineWidth',2,'Color',[0.26 0.26 0.26])

% add legend
legend('identity line','linear regression','location','southeast')
text(0.7,0.12,['y~',num2str(round(p(2),2)),'+',num2str(round(p(1),2)),'x'],'fontSize',14)
grid on,
