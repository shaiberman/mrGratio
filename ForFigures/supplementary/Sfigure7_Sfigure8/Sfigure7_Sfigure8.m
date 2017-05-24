clear   , close all

%% comparing NODDI based estinated to MVF and FA based estimes
currentPath=fileparts(which(mfilename));
cols=[[220,170,50]/255 ; [50,170,210]/255 ;[170,210,50]/255 ; [170,50,210]/255 ];
MarkSiz=40;
a=0.5;

% load data
load(fullfile(currentPath,'Sfigure78_NODDIvsFA.mat'))


% plot
for     ii=1:length(subData)
    
    
    fvf_fa=subData(ii).fvf_fa;
    fvf_noddi=subData(ii).fvf_noddi;
    fic=subData(ii).fic_noddi;
    
    %% plot FVF : supplementary figure 7
    
    figure(1),
    g=subplot(2,2,ii);
    axis([0.2 1 0.2 1]), axis square,hold on
    
    % adjust the subplot location
    pos = get(g,'position');
    pos(4) = pos(4)*1.10; pos(3) = pos(3)*1.10; % Add 10 percent to height and width
    if (ii==2 | ii==4),  pos(1) = pos(1)*0.85;  end % put right column closer to left
    set(g, 'position', pos);
    
    % plot the data
    scatter(fvf_fa, fvf_noddi,MarkSiz,  'markerEdgeColor','k','markerFaceColor',cols(ii,:));
    xlabel('FVF_F_A','fontSize',16), ylabel('FVF_N_O_D_D_I','fontSize',16)
    set(gca,'xTick',[0.2:0.2:1],'yTick',[0.2:0.2:1])
    
    r=corrcoef(fvf_fa, fvf_noddi);r=r(2);
    title(sprintf(['subject ',num2str(ii),'s : r=',num2str(round(r,2))]))
        
    % add linear fit
    x0=0.2:0.1:1;
    p=polyfit(fvf_fa, fvf_noddi,1);
    y0=p(1)*x0+p(2);
    plot(x0,y0,'k')
    
    identityLine(gca);box on
    l=legend('callosum data','identity line','linear fit','location','southeast');
    set(l,'fontSize',14)
    
    
    
    %% plot Vic(1-Viso), supplementary figure 8
    figure(2)
    g=subplot(2,2,ii);
    axis([0.1 0.9 0.1 0.9]),axis square,hold on
    
    % adjust the subplot location
    pos = get(g,'position');
    pos(4) = pos(4)*1.10; pos(3) = pos(3)*1.10; % Add 10 percent to height and width
    if (ii==2 | ii==4),  pos(1) = pos(1)*0.85;   end % put right column closer to left
    set(g, 'position', pos);
    
    % plot
    scatter(fvf_fa,fic,MarkSiz,  'markerEdgeColor','k','markerFaceColor',cols(ii,:));
    xlabel('FVF_F_A','fontSize',16), ylabel('V_i_c(1-V_i_s_o)','fontSize',16)
    set(gca,'xTick',[0.1:0.2:0.9],'yTick',[0.1:0.2:0.9])
    
    r=corrcoef(fvf_fa,fic);r=r(2);
    title(sprintf(['subject ',num2str(round(ii,2)),'s : r=',num2str(round(r,2))]))
    
    % add linear fit
    x0=0.1:0.1:0.9;
    p=polyfitZero(fvf_fa,fic,1);
    y0=p(1)*x0+p(2);
    plot(x0,y0,'k')
    identityLine(gca);box on
    
    l=legend('callosum data','identity line','linear fit','location','southeast');
    set(l,'fontSize',14)
    %
end

%%
 figure(1),set(gcf,'units','normalized','position',[.1 .1 .8 .8])
 figure(2),set(gcf,'units','normalized','position',[.1 .1 .8 .8])

