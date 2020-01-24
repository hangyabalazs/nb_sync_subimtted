function nbtonicphasic_pChAT(cb)
%NBTONICPHASIC   Group cells based on auto-correlation.
%   Cholinergic and putative cholinergic cells are grouped into tonic,
%   phasic - bursting and poisson-like groups based on refractory
%   and burst index calculated from their auto-correlograms (see NBACG for
%   details). ACGs are smoothed and normaized before averaging.
%
%   See also NBACG.

%   Balazs Hangya, Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu

% Load auto-correlations
%   Select corresponding ACG matrix

global RESDIR;
global PATH;
global cCode;
global plotWN;


% Select pooled CB mode
if strcmp(cb, 'POOLED')
    poolCBs = true;
else
    poolCBs = false;
    choosecb(cb);
end

% Load ACG file
if poolCBs
    fnm_acg = [RESDIR 'POOLED\acgPOOLED\ACG_matrices_POOLED.mat'];
    load(fnm_acg);   % Load ACG matrix
    % TP group indexes
    phasicB = Refractory < 40 & BurstIndex > 0.2;   % phasic, bursting
    poissonL = Refractory < 40 & BurstIndex <= 0.2;   % poisson-like
    tonic = Refractory >= 40;   % refractory above 40 ms
    
    % Define results directory
    resdir = [RESDIR 'POOLED\tonicphasicPOOLED\chat_pchat\'];
    
else
    
    fnm_acg = [RESDIR cb '\acg' cb PATH '\ACG_matrices_' cb '.mat'];
    load(fnm_acg);   % Load ACG matrix
    % Define results directory
    resdir = [RESDIR cb '\tonicphasic' cb PATH '\chat_pchat\'];
    
    % Tonic and phasic cells
    [phasicB, poissonL, tonic] = groupsTP(cellids);
end

win = plotWN;
allChAT = cellids;
NumChAT = length(allChAT);   % number of cells

if ~poolCBs
    mouse = getvalue('RatID',allChAT);  % animal ID
    nms = length(unique(mouse));   % number of mice (n=9)
    mouse2 = nan(size(mouse));
    for k = 1:nms
        mouse2(mouse==min(mouse)) = k;   % label for indiv. mice
        mouse(mouse==min(mouse)) = Inf;  % temp. variable
    end
    
    % Color code
    clr = colormap(jet(nms));   % unique color for each mouse
    close;
    depth = getvalue('DVpos', allChAT(1:45));
    Baseline = getvalue('Baseline', cellids);
end

if poolCBs
    nbNum = [];
    hdbNum = [];
    phdbNum = [];
    cbUsed = [0 0 0];
    for iCC = 1:NumChAT
        cbSwitcher(allChAT(iCC))
        currCB = whichcb;
        
        % Baseline
        VS = loadcb(allChAT{iCC},'Spikes');   % load prealigned spikes
        FR(iCC) = (length(VS)-1)/(VS(end)-VS(1));
        maxACG(iCC) = max(CCR(iCC,:));
        
        
        switch currCB
            case 'NB'
                ChAT(iCC,1) = logical(getvalue('ChAT+',allChAT(iCC)));
                pChAT(iCC,1) = logical(getvalue('pChAT+',allChAT(iCC)));
                nbNum = iCC;
                cbUsed(1) = 1;
            case 'HDB'
                ChAT(iCC,1) = true;
                pChAT(iCC,1) = false;
                hdbNum = iCC;
                cbUsed(2) = 1;
            case 'PannaHDB'
                ChAT(iCC,1) = true;
                pChAT(iCC,1) = false;
                phdbNum = iCC;
                cbUsed(3) = 1;
        end
    end
else
    FR = Baseline;
    ChAT = getvalue('ChAT+',allChAT);
    pChAT = getvalue('pChAT+',allChAT);
end

if isempty(pChAT)
    pChAT = zeros(NumChAT,1);
end
% coloring according to identified (red) and putative (green)ChATs
clr2 = [ChAT>0 pChAT>0 zeros(NumChAT,1)];   

if ~poolCBs
    % Plot ACG properties
    H1 = figure;
    hold on
    xlabel('Refractory (ms)')
    ylabel('Burst index')
    H2 = figure;
    hold on
    xlabel('Refractory (ms)')
    ylabel('Burst index')
    for k = 1:NumChAT
        figure(H1)
        plot(Refractory(k),BurstIndex(k),'o','MarkerSize',12,...
            'MarkerFaceColor',clr(mouse2(k),:),'MarkerEdgeColor','k')  % scatter plot
        title('Unique color for each mouse')
        setmyplot_tamas;
        figure(H2)
        plot(Refractory(k),BurstIndex(k),'o','MarkerSize',12,...
            'MarkerFaceColor',clr2(k,:),'MarkerEdgeColor','k')  % scatter plot
        title('ChAT (red) vs pChAT (green)')
        setmyplot_tamas;
    end
    fnm = ['chat_pchat_scatter_' cb '.fig'];
    saveas(H1,fullfile(resdir,fnm))   % save plot
    close;
    fnm = ['chat_pchat_scatter_' cb '.jpeg'];
    saveas(H1,fullfile(resdir,fnm))   % save plot
    close;
end

% Normalize ACG
CCRnorm = nan(size(CCR));
for k = 1:NumChAT
    sccr = smooth(CCR(k,:),'linear',5);   % smooth
    b = mean(sccr(lags>-500&lags<500));
    CCRnorm(k,:) = sccr / b;
end

% Tonic and phasic cells
tonicC = Refractory >= 40 & ChAT == 1;   % refractory above 40 ms
phasicBC = Refractory < 40 & BurstIndex > 0.2 & ChAT == 1;   % phasic, bursting
poissonLC = Refractory < 40 & BurstIndex <= 0.2 & ChAT == 1;   % poisson-like
tonicP = Refractory >= 40 & pChAT == 1;   % refractory above 40 ms
phasicBP = Refractory < 40 & BurstIndex > 0.2 & pChAT == 1;   % phasic, bursting
poissonLP = Refractory < 40 & BurstIndex <= 0.2 & pChAT == 1;   % poisson-like

tonicCCR_C = CCRnorm(tonicC,:);
phasicBCCR_C = CCRnorm(phasicBC,:);
poissonLCCR_C = CCRnorm(poissonLC,:);
tonicCCR_P = CCRnorm(tonicP,:);
phasicBCCR_P = CCRnorm(phasicBP,:);
poissonLCCR_P = CCRnorm(poissonLP,:);


% Plot tonic cells
%ChAT
[~, inx] = sort(Refractory(tonicC),'descend');   % sort according to refractory
figure
imagesc(lags(lags>win(1)&lags<win(2)),1:size(tonicCCR_C,1),...
    tonicCCR_C(inx,lags>win(1)&lags<win(2)))
set(gca,'CLim',[0 2.5])
title({'TONIC cholinergic cells'; 'sorted by Refractory'})
xlabel('Time (ms)')
ylabel('Cell number')
colormap jet
xticks([win(1), 0, win(2)]);
xlim(win);
setmyplot_tamas;
H2 = gcf;
fnm = ['tonic_ACG_ChAT_' cb '.fig'];
saveas(H2,fullfile(resdir,fnm))   % save plot
fnm = ['tonic_ACG_ChAT_' cb '.jpeg'];
saveas(H2,fullfile(resdir,fnm))   % save plot
close;
%pChAT
[~, inx] = sort(Refractory(tonicP),'descend');   % sort according to refractory
figure
imagesc(lags(lags>win(1)&lags<win(2)),1:size(tonicCCR_P,1),...
    tonicCCR_P(inx,lags>win(1)&lags<win(2)))
set(gca,'CLim',[0 2.5])
title({'TONIC cholinergic cells'; 'sorted by Refractory'})
xlabel('Time (ms)')
ylabel('Cell number')
colormap jet
xticks([win(1), 0, win(2)]);
xlim(win);
setmyplot_tamas;
H2 = gcf;
fnm = ['tonic_ACG_pChAT_' cb '.fig'];
saveas(H2,fullfile(resdir,fnm))   % save plot
fnm = ['tonic_ACG_pChAT_' cb '.jpeg'];
saveas(H2,fullfile(resdir,fnm))   % save plot
close;
% Plot phasicB cells
[~, inx] = sort(BurstIndex(phasicBC),'descend');   % sort according to burst index
figure
imagesc(lags(lags>win(1)&lags<win(2)),1:size(phasicBCCR_C,1),...
    phasicBCCR_C(inx,lags>win(1)&lags<win(2)))
set(gca,'CLim',[0 2.5])
title({'PHASIC BURSTING cholinergic cells'; 'sorted by BurstIndex'})
xlabel('Time (ms)')
ylabel('Cell number')
colormap jet
xticks([win(1), 0, win(2)]);
xlim(win);
setmyplot_tamas;
H3 = gcf;
fnm = ['phasicB_ACG_ChAT_' cb '.fig'];
saveas(H3,fullfile(resdir,fnm))   % save plot
fnm = ['phasicB_ACG_ChAT_' cb '.jpeg'];
saveas(H3,fullfile(resdir,fnm))   % save plot
close;

[~, inx] = sort(BurstIndex(phasicBP),'descend');   % sort according to burst index
figure
imagesc(lags(lags>win(1)&lags<win(2)),1:size(phasicBCCR_P,1),...
    phasicBCCR_P(inx,lags>win(1)&lags<win(2)))
set(gca,'CLim',[0 2.5])
title({'PHASIC BURSTING cholinergic cells'; 'sorted by BurstIndex'})
xlabel('Time (ms)')
ylabel('Cell number')
colormap jet
xticks([win(1), 0, win(2)]);
xlim(win);
setmyplot_tamas;
H3 = gcf;
fnm = ['phasicB_ACG_pChAT_' cb '.fig'];
saveas(H3,fullfile(resdir,fnm))   % save plot
fnm = ['phasicB_ACG_pChAT_' cb '.jpeg'];
saveas(H3,fullfile(resdir,fnm))   % save plot
close;

[~, inx] = sort(BurstIndex(poissonLC),'descend');   % sort according to burst index
figure
imagesc(lags(lags>win(1)&lags<win(2)),1:size(poissonLCCR_C,1),...
    poissonLCCR_C(inx,lags>win(1)&lags<win(2)))
set(gca,'CLim',[0 2.5])
title({'POISSON-LIKE cholinergic cells'; 'sorted by BurstIndex'})
xlabel('Time (ms)')
ylabel('Cell number')
colormap jet
xticks([win(1), 0, win(2)]);
xlim(win);
setmyplot_tamas;
H4 = gcf;
fnm = ['poissonL_ACG_ChAT_' cb '.fig'];
saveas(H4,fullfile(resdir,fnm))   % save plot
fnm = ['poissonL_ACG_ChAT_' cb '.jpeg'];
saveas(H4,fullfile(resdir,fnm))   % save plot
close;

[~, inx] = sort(BurstIndex(poissonLP),'descend');   % sort according to burst index
figure
imagesc(lags(lags>win(1)&lags<win(2)),1:size(poissonLCCR_P,1),...
    poissonLCCR_P(inx,lags>win(1)&lags<win(2)))
set(gca,'CLim',[0 2.5])
title({'POISSON-LIKE cholinergic cells'; 'sorted by BurstIndex'})
xlabel('Time (ms)')
ylabel('Cell number')
colormap jet
xticks([win(1), 0, win(2)]);
xlim(win);
setmyplot_tamas;
H4 = gcf;
fnm = ['poissonL_ACG_pChAT_' cb '.fig'];
saveas(H4,fullfile(resdir,fnm))   % save plot
fnm = ['poissonL_ACG_pChAT_' cb '.jpeg'];
saveas(H4,fullfile(resdir,fnm))   % save plot
close;

% Plot averages for ChAT, pChAT
%ChAT
% average ACG, tonic cells
mn_tonic_C = mean(tonicCCR_C(:,lags>win(1)&lags<win(2)));  
% SE, tonic cells
se_tonic_C = std(tonicCCR_C(:,lags>win(1)&lags<win(2))) /...
    sqrt(size(tonicCCR_C,1));   
% average ACG, phasic, bursting cells
mn_phasicB_C = mean(phasicBCCR_C(:,lags>win(1)&lags<win(2)));  
% SE, phasic, bursting cells
se_phasicB_C = std(phasicBCCR_C(:,lags>win(1)&lags<win(2))) /...
    sqrt(size(phasicBCCR_C,1));   
% average ACG, poisson-like cells
mn_poissonL_C = mean(poissonLCCR_C(:,lags>win(1)&lags<win(2)));   
% SE, poisson-like cells
se_poissonL_C = std(poissonLCCR_C(:,lags>win(1)&lags<win(2))) /...
    sqrt(size(poissonLCCR_C,1));   

% SUPPLEMENTARY FIGURE 1 PANEL A
% ChAT
H4 = figure;
hold on
errorshade(lags(lags>win(1)&lags<win(2)),mn_phasicB_C,se_phasicB_C,...
    'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',2)
errorshade(lags(lags>win(1)&lags<win(2)),mn_poissonL_C,se_poissonL_C,...
    'LineColor',cCode(2,:),'ShadeColor',cCode(2,:),'LineWidth',2)
errorshade(lags(lags>win(1)&lags<win(2)),mn_tonic_C,se_tonic_C,...
    'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',2)
xlim(win);
xticks([win(1), 0, win(2)]);
xlabel('Lag (ms)')
ylabel('Normalized ACG')
axis square;
title('ChAT+');
legend({['n=' num2str(sum(phasicBC))],['n=' num2str(sum(poissonLC))],...
    ['n=' num2str(sum(tonicC))]}, 'FontSize', 12, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
fnm = ['AVG_ACG_groups_ChAT_' cb '.fig'];
saveas(H4,fullfile(resdir,fnm))   % save plot
fnm = ['AVG_ACG_groups_ChAT_' cb '.jpeg'];
saveas(H4,fullfile(resdir,fnm))   % save plot
close(H4);

% pChAT
% average ACG, tonic cells
mn_tonic_P = mean(tonicCCR_P(:,lags>win(1)&lags<win(2)));  
% SE, tonic cells
se_tonic_P = std(tonicCCR_P(:,lags>win(1)&lags<win(2))) /...
    sqrt(size(tonicCCR_P,1));   
% average ACG, phasic, bursting cells
mn_phasicB_P = mean(phasicBCCR_P(:,lags>win(1)&lags<win(2))); 
% SE, phasic, bursting cells
se_phasicB_P = std(phasicBCCR_P(:,lags>win(1)&lags<win(2))) /...
    sqrt(size(phasicBCCR_P,1));   
% average ACG, poisson-like cells
mn_poissonL_P = mean(poissonLCCR_P(:,lags>win(1)&lags<win(2))); 
% SE, poisson-like cells
se_poissonL_P = std(poissonLCCR_P(:,lags>win(1)&lags<win(2))) /...
    sqrt(size(poissonLCCR_P,1));   

H4 = figure;
hold on
errorshade(lags(lags>win(1)&lags<win(2)),mn_phasicB_P,se_phasicB_P,...
    'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',2)
errorshade(lags(lags>win(1)&lags<win(2)),mn_poissonL_P,se_poissonL_P,...
    'LineColor',cCode(2,:),'ShadeColor',cCode(2,:),'LineWidth',2)
errorshade(lags(lags>win(1)&lags<win(2)),mn_tonic_P,se_tonic_P,...
    'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',2)
xlim(win);
xticks([win(1), 0, win(2)]);
xlabel('Lag (ms)')
ylabel('Normalized ACG')
axis square;
title('pChAT+');
legend({['n=' num2str(sum(phasicBP))],['n=' num2str(sum(poissonLP))],...
    ['n=' num2str(sum(tonicP))]}, 'FontSize', 12, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
fnm = ['AVG_ACG_groups_pChAT_' cb '.fig'];
saveas(H4,fullfile(resdir,fnm))   % save plot
fnm = ['AVG_ACG_groups_pChAT_' cb '.jpeg'];
saveas(H4,fullfile(resdir,fnm))   % save plot
close;

% SUPPLEMENTARY FIGURE1 PANEL D
% Plot scatter for Tonic - Phasic groups BI vs Ref
%ChAT vs pChAT

%Figure
H9 = figure;
hold on
xlabel('Refractory (ms)')
ylabel('Burst Index')
plot(Refractory(tonicC)',BurstIndex(tonicC)','o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:))
plot(Refractory(tonicP)',BurstIndex(tonicP)','^','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:))
plot(Refractory(phasicBC),BurstIndex(phasicBC),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(Refractory(poissonLC),BurstIndex(poissonLC),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
plot(Refractory(phasicBP),BurstIndex(phasicBP),'^','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(Refractory(poissonLP),BurstIndex(poissonLP),'^','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
title('   ')
x_lim = xlim;
y_lim = ylim;
y_lim(2) = 1;
ylim(y_lim);
xticks([0 x_lim(2)/2 x_lim(2)])
yticks([-1 0 1]);
axis square
setmyplot_tamas;
fnm = ['BurstIndex_vs_Refractory' cb '.fig'];
saveas(H9,fullfile(resdir,fnm))   % save plot
fnm = ['BurstIndex_vs_Refractory' cb '.jpeg'];
saveas(H9,fullfile(resdir,fnm))   % save plot
close;

% SUPPLEMENTARY FIGURE1 PANEL E
% Plot scatter for Tonic - Phasic groups TI vs Ref
%ChAT vs pChAT
% Regression
y = ThetaIndex_IN(tonic);
x = Refractory(tonic)';
[b,~,~,~,stats] = regress(y',[ones(length(x),1),x']);
p = stats(3);           % F-test significance

%Figure
H9 = figure;
hold on
xlabel('Refractory')
ylabel('ThetaIndex')
plot(Refractory(tonicC)',ThetaIndex_IN(tonicC)','o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:))
plot(Refractory(tonicP)',ThetaIndex_IN(tonicP)','^','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:))
axis square
icp = b(1);   % intercept
gr = b(2);   % gradient
xx = 0:0.01:max(x);
yy = xx .* gr + icp;
hold on
plot(xx,yy, '--', 'Color',cCode(3,:),'LineWidth',2)   % overlay regression line
plot(Refractory(phasicBC),ThetaIndex_IN(phasicBC),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(Refractory(poissonLC),ThetaIndex_IN(poissonLC),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
plot(Refractory(phasicBP),ThetaIndex_IN(phasicBP),'^','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(Refractory(poissonLP),ThetaIndex_IN(poissonLP),'^','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
title('   ')
x_lim = xlim;
y_lim = ylim;
y_lim(2) = 1;
ylim(y_lim);
xticks([0 x_lim(2)/2 x_lim(2)])
text((x_lim(2)*0.45), y_lim(2)*0.95, ['p=' num2str(round(p,2))], 'Color',...
    cCode(3,:));
yticks([0 0.5 1]);
setmyplot_tamas;
fnm = ['ThetaIndex_vs_Refractory' cb '.fig'];
saveas(H9,fullfile(resdir,fnm))   % save plot
fnm = ['ThetaIndex_vs_Refractory' cb '.jpeg'];
saveas(H9,fullfile(resdir,fnm))   % save plot
close;

% SUPPLEMENTARY FIGURE1 PANEL F
% Baseline for ChAT and pChAT
Baseline = FR;
if ~isempty(Baseline)
    H10 = figure;
    hold on;
    bar([1, 2, 3], [nanmedian(Baseline(phasicB)) nanmedian(Baseline(poissonL))...
        nanmedian(Baseline(tonic))], 'FaceColor', 'none', 'LineWidth',2)
    plot(0.75, Baseline(phasicBC), 'o', 'MarkerFaceColor', 'none',...
        'MarkerEdgeColor',cCode(1,:), 'LineWidth',1)
    plot(1.75, Baseline(poissonLC), 'o', 'MarkerFaceColor', 'none',...
        'MarkerEdgeColor',cCode(2,:), 'LineWidth',1)
    plot(2.75, Baseline(tonicC), 'o', 'MarkerFaceColor', 'none',...
        'MarkerEdgeColor',cCode(3,:), 'LineWidth',1)
    plot(1.25, Baseline(phasicBP), '^', 'MarkerFaceColor', 'none',...
        'MarkerEdgeColor',cCode(1,:), 'LineWidth',1)
    plot(2.25, Baseline(poissonLP), '^', 'MarkerFaceColor', 'none',...
        'MarkerEdgeColor',cCode(2,:), 'LineWidth',1)
    plot(3.25, Baseline(tonicP), '^', 'MarkerFaceColor', 'none',...
        'MarkerEdgeColor',cCode(3,:), 'LineWidth',1)
    ylabel('Firing rate (Hz)');
    xlabel('Groups')
    xticks([1 2 3])
    title({'ChAT - circle'; 'pChAT - triangle'});
    y_lim = ylim;
    y_lim = ceil(y_lim/10)*10;
    ylim([y_lim])
    yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
    xticks([])
    axis square;
    setmyplot_tamas;
    saveas(H10, [resdir 'ChAT_BaselineFR_median.fig']);
    saveas(H10, [resdir 'ChAT_BaselineFR_median.jpeg']);
    close(H10);
end


