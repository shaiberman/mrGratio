function f = nc_PlotModelFits_gender(coefs1,coefs2,valName,fgNames,fgnums,color1,color2,R1,R2,pVals)
% Make figures from coeficient estimates (nc_ModelSelection)
%
% f = nc_PlotModelFits(coefs,valName,fgNames,fgnums,color)
%
% Inputs:
%
% coefs   - Structure containing model coefficient estimates (see
%           nc_ModelSelection)
% valName - The name of the qMR parameter (e.g., R1)
% fgNames - Cell array of fiber group names (in the order they are in the
%           AFQ structure)
% fgnums  - Vector of fiber group numbers. The plots will be in the order
%           denoted here. For example if fgnums(1) = 3, than the third
%           fiber group (Left CST) will be plotted first.
% color   - rgb values for each plot
%
% Copyright Jason D. Yeatman, August 2014. Code released with:
% Yeatman JD, Wandell BA & Mezer AM (2014). Lifespan maturation
% and degeneration of human brain white matter. Nature Communications

% Open figure window
figure;

%% Plotting options

if notDefined('valName')
    valName = 'T1_map_lsq_2DTI';
end

% Which fiber tracts to plot
if notDefined('fgnums')
    fgnums = 1:20;
end

% Get the fiber group names
if notDefined('fgNames')
    fgNames = {'Left Thalamic Radiation','Right Thalamic Radiation','Left Corticospinal','Right Corticospinal', 'Left Cingulum Cingulate', 'Right Cingulum Cingulate'...
        'Left Cingulum Hippocampus','Right Cingulum Hippocampus', 'Callosum Forceps Major', 'Callosum Forceps Minor'...
        'Left IFOF','Right IFOF','Left ILF','Right ILF','Left SLF','Right SLF','Left Uncinate','Right Uncinate','Left Arcuate','Right Arcuate'...'
        'CC_Occipital',  'CC_Temporal',  'CC_Post_Parietal',   'CC_Sup_Parietal',  'CC_Motor',  'CC_Sup_Frontal',  'CC_Ant_Frontal',  'CC_Orb_Frontal'};
end
fgNames = fgNames(fgnums);

% axis scaling
if strcmp(valName,'MTV')
    yticks = [.22 .3  .38];
    ylims = [.2 .4];
    YLab = 'MTV'
    xtl = .01;
    mat = 'max';
elseif strcmp(valName, 'md')
    ylims = [.53 .85];
    yticks = [.6 .7 .8];
    YLab = 'ADC (\mum^2/ms)';
    xtl = .013;
    mat = 'min';
elseif strcmp(valName,'fa')
       ylims = [.58 .95];
    yticks = [.65 .75 .85];
       YLab = 'FA';
    xtl = .013;
    mat = 'max';
elseif strcmp(valName,'T1_map_lsq') || strcmp(valName,'T1_map_lsq_2DTI')
    ylims = [.8 1.4];
    yticks=[.8 1 1.2 1.4];
    YLab = 'T1 (seconds)';
    xtl = .03;
    mat = 'min';
elseif strcmp(valName,'R1_2DTI') || strcmp(valName,'R1')
    ylims = [.8 1.2];
    yticks=[.9 1 1.1];
    YLab = 'R1';
    xtl = .02;
    mat = 'max';
elseif strcmp(valName,'gr') || strcmp(valName,'gratio')
    ylims = [.38 .9];
    yticks=[.4 .6 .8];
    YLab = 'g-ratio';
    xtl = .02;
    mat = 'max';
elseif strcmp(valName,'FVF') || strcmp(valName,'fvf')
    ylims = [.28 .8];
    yticks=[.3 .5 .7];
    YLab = 'FVF';
    xtl = .02;
    mat = 'max';
end
xticks = 8:10:78;
x0 = min(coefs1(1).x):max(coefs2(1).x);
xlims = [min(x0)-1 max(x0)+1];
        
%% Loop over tracts and plot
pnum = 0;
for ii = fgnums
    % count
    pnum = pnum+1;
    % Open a plot
    axh(ii)=subplot(3,3,pnum);
    hold on;
    for jj=1:2
        clear bootCI
        eval([' coefs = coefs' num2str(jj) ';'])
        eval([' color = color' num2str(jj) ';'])
%         eval([' shape = shape' num2str(jj) ';'])
        eval([' R = R' num2str(jj) ';'])
        % Points to evaluate model
     
        % plot the data
        plot(coefs(ii).x, coefs(ii).y,'o','color',[0 0 0],'markerfacecolor',color(pnum,:),'markersize',6);
        
        % compute model prediction
        switch(coefs(ii).name)
            case {'lowess' 'lowess21' 'lowess22'}
                yhat = coefs(ii).full(:,2);
                x0   = coefs(ii).full(:,1)';
                bootCI = coefs(ii).boot;
                v95 = [nan nan];
            case {'piecewise' 'piecewise2' 'piecewisenoflat'}
                yhat = piecewiseEval(coefs(ii).full,x0);
                % Calculate confidence intervals for each bootstrap iteration
                for kk = 1:size(coefs(ii).boot,1)
                    bootCI(kk,:) = piecewiseEval(coefs(ii).boot(kk,:),x0);
                end
                v95 = [nan nan];
            case {'quadratic'}
                yhat = polyval(coefs(ii).full,x0);
                % Calculate confidence intervals for each bootstrap iteration
                for kk = 1:size(coefs(ii).boot,1)
                    bootCI(kk,:) = polyval(coefs(ii).boot(kk,:),x0);
                    % and confidence interval for vertex
                    vCI(kk) = -(coefs(ii).boot(kk,2)./(2*coefs(ii).boot(kk,1)));
                end
                v95 =  [nan nan];
            case {'poisson'}
                yhat = evalPoissonCurve(coefs(ii).full,x0);
                % Calculate confidence intervals for each bootstrap iteration
                for kk = 1:size(coefs(ii).boot,1)
                    bootCI(kk,:) = evalPoissonCurve(coefs(ii).boot(kk,:),x0);
                    % and confidence interval for vertex
                    vCI(kk) = 1./coefs(ii).boot(kk,2);
                end
                v95 =  [nan nan];
            case {'linear'}
                yhat  = coefs(ii).full(1) .* x0 + coefs(ii).full(2);
                % Calculate confidence intervals for each bootstrap iteration
                for kk = 1:size(coefs(ii).boot,1)
                    bootCI(kk,:) = coefs(ii).boot(kk,1) .* x0 + coefs(ii).boot(kk,2);
                end
                v95 = [nan nan];
            case {'exponent'}
                expfun = @(p,x) p(1).*x.^p(2) + p(3);
                yhat =   feval(expfun, coefs(ii).full, x0);
                % Calculate confidence intervals for each bootstrap iteration
                for kk = 1:size(coefs(ii).boot,1)
                    bootCI(kk,:) = feval(expfun, coefs(ii).boot(kk,:), x0);
                end
                v95 = [nan nan];
        end
        
        % Compute the 95% confidence intervals on the fits
        p95{ii} = prctile(bootCI,[2.5 97.5]);
        % Plot the 95% confidence interval
        fill([x0 fliplr(x0)],[p95{ii}(1,:) fliplr(p95{ii}(2,:))],color(pnum,:),'facealpha',.5,'edgealpha',0);
        % Plot the model fit
        plot(x0,yhat, 'color',color(pnum,:),'linewidth',2);
        
        %% Redraw the axes because matlab has a bug with opengl
    plot([xlims(1)+.01 xlims(1)+.01 xlims(2)],[ylims(2) ylims(1)+.0001 ylims(1)+.0001],'-k');

        %% Draw the confidence interval around the vertex
        plot(v95,repmat(min(ylims).*1.05,1,2),'linewidth',3,'color',color(pnum,:));
        plot(repmat(v95(1),1,2),[min(ylims).*1.03 min(ylims).*1.07],'linewidth',2,'color',color(pnum,:));
        plot(repmat(v95(2),1,2),[min(ylims).*1.03 min(ylims).*1.07],'linewidth',2,'color',color(pnum,:));
        
    end
    
    if ismember(pnum,6:8)
        xlabel('Age','fontname','times','fontsize',17);
        set(gca,'xtick',xticks,'fontsize',15);
    else
        set(gca,'xtick',xticks,'xticklabel',[]);
    end
    if ismember(pnum,1:3:24)
        ylabel(YLab,'fontname','times','fontsize',16);
        set(gca,'ytick',yticks,'fontsize',15)
    else
        set(gca,'ytick',yticks,'yticklabel',[]);
    end
    axis([xlims ylims]);
    
    % add a title
    if ~notDefined('pVals')
    if strcmp(valName,'gr') || strcmp(valName,'gratio')
               tit=['P = ',num2str(round(pVals(ii),4))];
        if strcmp(mat,'max')
            th = text((xlims(2)-xlims(1))./2 + xlims(1),min(ylims).*1.2,tit,...
                'fontsize',15,'fontname','times','HorizontalAlignment','center');
        elseif strcmp(mat,'min')
            th = text((xlims(2)-xlims(1))./2 + xlims(1),max(ylims).*.97,tit,...
                'fontsize',15,'fontname','times','HorizontalAlignment','center');
        end
    end
    end
    h=title(fgNames{pnum},'fontname','times','fontsize',14);
    if strcmp(valName,'gr') || strcmp(valName,'gratio') || strcmp(valName,'FVF')
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-0.1 P(3)])
    elseif strcmp(valName,'MTV')
        P = get(h,'Position');
        set(h,'Position',[P(1) P(2)-0.05 P(3)])
    end
end

for ii = fgnums
    % Get the position of each subplot
    p = get(axh(ii),'position');
    % Add this much
    pnew = [0 0 .03 .05];
    p = p+pnew;
    set(axh(ii),'position',p);
    
end

%% Set figure properties
f(1)=gcf;
set(gcf,'inverthardcopy','off','color',[1 1 1]);

