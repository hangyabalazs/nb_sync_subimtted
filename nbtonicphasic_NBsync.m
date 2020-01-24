function nbtonicphasic_NBsync(cb, varargin)
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
% Select corresponding ACG matrix
global RESDIR;
global PATH;
global cCode;
global plotWN;


% Select pooled CB mode
if strcmp(cb, 'POOLED')
    poolCBs = true;
else
    poolCBs = false;
end

if nargin>1
    fnm_acg = [RESDIR cb '\acg' cb PATH '\longSeg\ACG_matrices_' cb '.mat'];
else
    fnm_acg = [RESDIR cb '\acg' cb PATH '\ACG_matrices_' cb '.mat'];
end

% Load ACG file
load(fnm_acg);   % Load ACG matrix
% Define results directory
resdir = [RESDIR cb '\tonicphasic' cb PATH '\'];

% TP group indexes
phasicB = Refractory < 40 & BurstIndex > 0.2;   % phasic, bursting
poissonL = Refractory < 40 & BurstIndex <= 0.2;   % poisson-like
tonic = Refractory >= 40;   % refractory above 40 ms

plotStr= {['BURST-SB n=',num2str(sum(phasicB))];
    ['BURST-PL n=',num2str(sum(poissonL))];
    ['REG n=', num2str(sum(tonic))]
    };
win = plotWN;

% Example ACG-s selected for presentation
[Bindex, Pindex, Tindex] = exampleACGs(cb, cellids, Refractory,...
    BurstIndex);
allChAT = cellids;
NumChAT = length(allChAT);   % number of cells

for iC = 1:NumChAT
    cbSwitcher(allChAT{iC})
    VS = loadcb(allChAT{iC},'Spikes');   % load prealigned spikes
    FR(iC) = (length(VS)-1)/(VS(end)-VS(1)); % Firing rate
    maxACG(iC) = max(CCR(iC,:));
    posAP(iC) = getvalue('APpos', allChAT{iC}); % AP position
end

% Normalize ACG
CCRnorm = nan(size(CCR));
for k = 1:NumChAT
    sccr = smooth(CCR(k,:),'linear',5);   % smooth
    b = mean(sccr(lags>-500&lags<500));
    CCRnorm(k,:) = sccr / b;
end

%% FIGURE1 PANEL G2: Plot tonic cells
if sum(tonic)>0 % No tonic for HDB
    tonicCCR = CCRnorm(tonic,:);  % ACG for tonic cells
    [srt, inx] = sort(Refractory(tonic),'descend');   % sort according to refractory
    H2 = figure;
    imagesc(lags(lags>win(1)&lags<win(2)),1:size(tonicCCR,1),...
        tonicCCR(inx,lags>win(1)&lags<win(2)))
    set(gca,'CLim',[0 2.5])
    title({'Regular'})
    xlabel('Lag (ms)')
    ylabel('Neuron#')
    colormap jet
    xticks([win(1), 0, win(2)]);
    if sum(Tindex) > 0
        yticks([1, find(srt==Refractory(Tindex)), length(srt)]);
        yticklabels({'1', 'T', num2str(length(srt))});
    else
        yticks([1 length(srt)]);
        yticklabels({num2str(1), num2str(length(srt))});
    end
    xlim(win);
%     axis square;
    setmyplot_tamas;
    
    fnm = ['tonic_ACG' cb '.fig'];
    saveas(H2,fullfile(resdir,fnm))   % save plot
    fnm = ['tonic_ACG' cb '.jpeg'];
    saveas(H2,fullfile(resdir,fnm))   % save plot
    close(H2);
end

%% FIGURE1 PANEL F3: Plot phasicB cells
if sum(phasicB)>0
    phasicBCCR = CCRnorm(phasicB,:);  % ACG for phasic bursting cells
    [srt, inx] = sort(BurstIndex(phasicB),'descend');   % sort according to burst index
    H3 = figure;
    imagesc(lags(lags>win(1)&lags<win(2)),1:size(phasicBCCR,1),...
        phasicBCCR(inx,lags>win(1)&lags<win(2)))
    set(gca,'CLim',[0 2.5])
    title({'Strongly bursting'})
    xlabel('Lag (ms)')
    ylabel('Neuron#')
    colormap jet
    xticks([win(1), 0, win(2)]);
    if sum(Bindex) > 0
        yticks([1, find(srt==BurstIndex(Bindex)), length(srt)]);
        yticklabels({'1', 'B', num2str(length(srt))});
    else
        yticks([1 length(srt)]);
        yticklabels({num2str(1), num2str(length(srt))});
    end
    xlim(win);
%     axis square;
    setmyplot_tamas;
    fnm = ['phasicB_ACG' cb '.fig'];
    saveas(H3,fullfile(resdir,fnm))   % save plot
    fnm = ['phasicB_ACG' cb '.jpeg'];
    saveas(H3,fullfile(resdir,fnm))   % save plot
    close(H3);
end
%% FIGURE1 PANEL F4: Plot poisson cells
if sum(poissonL)>0
    poissonLCCR = CCRnorm(poissonL,:);  % ACG for poisson-like cells
    [srt, inx] = sort(BurstIndex(poissonL),'descend');   % sort according to burst index
    H4 = figure;
    imagesc(lags(lags>win(1)&lags<win(2)),1:size(poissonLCCR,1),...
        poissonLCCR(inx,lags>win(1)&lags<win(2)))
    set(gca,'CLim',[0 2.5])
    title({'Poisson-like'})
    xlabel('Lag (ms)')
    ylabel('Neuron#')
    colormap jet
    xticks([win(1), 0, win(2)]);
    if sum(Pindex) > 0
        y_ticks = unique([1, find(srt==BurstIndex(Pindex)), length(srt)]);
        yticks(y_ticks);
        yticklabels({'1', 'P', num2str(length(srt))});
    else
        yticks([1 length(srt)]);
        yticklabels({num2str(1), num2str(length(srt))});
    end
    xlim(win);
%     axis square;
    setmyplot_tamas;
    fnm = ['poissonL_ACG' cb '.fig'];
    saveas(H4,fullfile(resdir,fnm))   % save plot
    fnm = ['poissonL_ACG' cb '.jpeg'];
    saveas(H4,fullfile(resdir,fnm))   % save plot
    close(H4);
end

%% FIGURE1 PANEL H: Plot averages
if sum(tonic)>0
    mn_tonic = mean(tonicCCR(:,lags>win(1)&lags<win(2)));   % average ACG, tonic cells
    se_tonic = std(tonicCCR(:,lags>win(1)&lags<win(2))) /...
        sqrt(size(tonicCCR,1));   % SE, tonic cells
end
if sum(phasicB)>0
    mn_phasicB = mean(phasicBCCR(:,lags>win(1)&lags<win(2)));   % average ACG, phasic, bursting cells
    se_phasicB = std(phasicBCCR(:,lags>win(1)&lags<win(2))) /...
        sqrt(size(phasicB,1));   % SE, phasic, bursting cells
end
if sum(poissonL)>0
    mn_poissonL = mean(poissonLCCR(:,lags>win(1)&lags<win(2)));   % average ACG, poisson-like cells
    se_poissonL = std(poissonLCCR(:,lags>win(1)&lags<win(2))) /...
        sqrt(size(poissonL,1));   % SE, poisson-like cells
end

if poolCBs
    % Load untagged NB cells
    [mn_tonic_UT, se_tonic_UT, numT] = nbUntagged(lags, win);
    plotStr{4,1}= ['NB REG n=' num2str(numT)];
end

H5 = figure;
hold on

if poolCBs
    errorshade(lags(lags>win(1)&lags<win(2)),mn_tonic_UT,se_tonic_UT,...
        'LineColor',[0.6 0.6 0.6],'ShadeColor',[0.6 0.6 0.6],'LineWidth',2)
end

if sum(tonic)>0
    errorshade(lags(lags>win(1)&lags<win(2)),mn_tonic,se_tonic,...
        'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',2)
end
if sum(phasicB)>0
    errorshade(lags(lags>win(1)&lags<win(2)),mn_phasicB,se_phasicB,...
        'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',2)
end
if sum(poissonL)>0
    errorshade(lags(lags>win(1)&lags<win(2)),mn_poissonL,se_poissonL,...
        'LineColor',cCode(2,:),'ShadeColor',cCode(2,:),'LineWidth',2)
end
xlim(win);
xticks([win(1), 0, win(2)]);
xlabel('Lag (ms)')
ylabel('Normalized count')
axis square;
title('  ');
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
fnm = ['AVG_ACG_groups' cb '.fig'];
saveas(H5,fullfile(resdir,fnm))   % save plot
fnm = ['AVG_ACG_groups' cb '.jpeg'];
saveas(H5,fullfile(resdir,fnm))   % save plot
close(H5);

plotStr = plotStr(1:3,1); % Revert plotStr to original
%% FIGURE1 PANEL I: Plot scatter for Tonic - Phasic groups
% Also FIGURE7 PANEL B
H6 = figure;
hold on
if ~strcmp(cb, 'HDB')
    plot(Refractory(tonic),BurstIndex(tonic),'o','MarkerSize',12,...
        'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:))  % scatter plot
    plot(Refractory(Tindex),BurstIndex(Tindex),'o','MarkerSize',12,...
        'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:), 'LineWidth', 2)  % scatter plot
end
plot(Refractory(phasicB),BurstIndex(phasicB),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(Refractory(Bindex),BurstIndex(Bindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:), 'LineWidth', 2)  % scatter plot
plot(Refractory(poissonL),BurstIndex(poissonL),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
plot(Refractory(Pindex),BurstIndex(Pindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:), 'LineWidth', 2)  % scatter plot
title('   ')
xlabel('Refractory (ms)')
ylabel('Burst Index')
xlim([0 150]);
ylim([-1 1])
axis square;
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
fnm = ['Tonic_phasic_groups' cb '.fig'];
saveas(H6,fullfile(resdir,fnm))   % save plot
fnm = ['Tonic_phasic_groups' cb '.jpeg'];
saveas(H6,fullfile(resdir,fnm))   % save plot
close(H6);

keyboard; % Continue in case of POOLED data

%% FIGURE7 PANEL C1: Plot scatter for Tonic - Phasic groups BI vs AP

% SB AP
regAP = posAP(tonic);
burstAP = posAP(phasicB | poissonL);
% SB BI
burstBI = BurstIndex(phasicB | poissonL);
regBI = BurstIndex(tonic);
% SB REF
burstREF = Refractory(phasicB | poissonL);
regREF = Refractory(tonic);

% Calculate means grouped by AP positions
distS = unique(posAP);
for iD = 1:numel(distS)
    burstISum(iD) = mean(BurstIndex(posAP==distS(iD)));
    refSum(iD) = mean(Refractory(posAP==distS(iD)));
    burstSum(iD) = mean(burstBI(burstAP==distS(iD)));
    regSum(iD) = mean(regBI(regAP==distS(iD)));
    bRSum(iD) = mean(burstREF(burstAP==distS(iD)));
    rRSum(iD) = mean(regREF(regAP==distS(iD)));
end

% replace NaN-s
regSum(isnan(regSum)) = 0;
burstSum(isnan(burstSum)) = 0;
bRSum(isnan(bRSum)) = 0;
rRSum(isnan(rRSum)) = 0;

% 3-point moving average
ratioBI = movmean(burstISum,3);
ratioRef = movmean(refSum,3);
ratioRB = (regSum+burstSum)/2;
smoothRB = movmean(ratioRB,3);

% Line colors
lineCol = [255 127 127]./255;

H6 = figure;
hold on
if ~strcmp(cb, 'HDB')
    plot(posAP(tonic),BurstIndex(tonic),'o','MarkerSize',12,...
        'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:))  % scatter plot
    plot(posAP(Tindex),BurstIndex(Tindex),'o','MarkerSize',12,...
        'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:), 'LineWidth', 2)  % scatter plot
end
plot(posAP(phasicB),BurstIndex(phasicB),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(posAP(Bindex),BurstIndex(Bindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:), 'LineWidth', 2)  % scatter plot
plot(posAP(poissonL),BurstIndex(poissonL),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
plot(posAP(Pindex),BurstIndex(Pindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:), 'LineWidth', 2)  % scatter plot
plot(unique(posAP), ratioBI, 'Color', lineCol, 'LineWidth', 2)
title('   ')
yticks([-1 0 1])
xticks([-1000 0 1000])
xticklabels({'-1', '0', '1'});
xlabel('AP dist. from bregma (mm)')
ylabel('Burst Index')
xlim([-1500 1500]);
ylim([-1 1])
axis square;
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
fnm = ['Tonic_phasic_groups_BI_AP' cb '.fig'];
saveas(H6,fullfile(resdir,fnm))   % save plot
fnm = ['Tonic_phasic_groups_BI_AP' cb '.jpeg'];
saveas(H6,fullfile(resdir,fnm))   % save plot
close(H6);


%% FIGURE7 PANEL C2: Plot scatter for Tonic - Phasic groups BI vs AP
H6 = figure;
hold on
if ~strcmp(cb, 'HDB')
    plot(posAP(tonic),Refractory(tonic),'o','MarkerSize',12,...
        'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:))  % scatter plot
    plot(posAP(Tindex),Refractory(Tindex),'o','MarkerSize',12,...
        'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:), 'LineWidth', 2)  % scatter plot
end
plot(posAP(phasicB),Refractory(phasicB),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(posAP(Bindex),Refractory(Bindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:), 'LineWidth', 2)  % scatter plot
plot(posAP(poissonL),Refractory(poissonL),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
plot(posAP(Pindex),Refractory(Pindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:), 'LineWidth', 2)  % scatter plot
plot(unique(posAP), ratioRef, 'Color', lineCol, 'LineWidth', 2)
title('   ')
xlabel('posAP (um)')
ylabel('Refractory')
xticks([-1000 0 1000])
xticklabels({'-1', '0', '1'});
xlabel('AP dist. from bregma (mm)')
xlim([-1500 1500]);
ylim([0 150])
axis square;
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
fnm = ['Tonic_phasic_groups_Ref_AP' cb '.fig'];
saveas(H6,fullfile(resdir,fnm))   % save plot
fnm = ['Tonic_phasic_groups_Ref_AP' cb '.jpeg'];
saveas(H6,fullfile(resdir,fnm))   % save plot
close(H6);


%% FIGURE1 PANEL J: NEW CALCULATION TYPE LONG
% Regression
if ~strcmp(cb, 'HDB')
    y = ThetaIndex_IN(tonic);
    x = Refractory(tonic)';
    [b,~,~,~,stats] = regress(y',[ones(length(x),1),x']);
    p = round(stats(3),4);           % F-test significance
end
H10 = figure;
hold on
if ~strcmp(cb, 'HDB')
    plot(x,y,'o','MarkerSize',12,'MarkerFaceColor','none',...
        'MarkerEdgeColor',cCode(3,:))
    plot(x(Tindex),y(Tindex),'o','MarkerSize',12,'MarkerFaceColor',...
        'none','MarkerEdgeColor',cCode(3,:), 'LineWidth', 2)  % scatter plot
    icp = b(1);   % intercept
    gr = b(2);   % gradient
    xx = 0:0.01:max(x);
    yy = xx .* gr + icp;
end
axis square
hold on
if ~strcmp(cb, 'HDB')
    plot(xx,yy, '--', 'Color',cCode(3,:),'LineWidth',2)   % overlay regression line
end
plot(Refractory(phasicB),ThetaIndex_IN(phasicB),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(Refractory(poissonL),ThetaIndex_IN(poissonL),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
plot(Refractory(Pindex),ThetaIndex_IN(Pindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:), 'LineWidth', 2)  % scatter plot
plot(Refractory(Bindex),ThetaIndex_IN(Bindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:), 'LineWidth', 2)  % scatter plot
title('   ')
xlabel('Refractory (ms)')
ylabel('Theta Index')
x_lim = xlim;
y_lim = ylim;
if ~strcmp(cb, 'HDB')
    text((x_lim(2)*0.45), y_lim(2)*0.95, ['p: ' num2str(p)], 'Color',...
        cCode(3,:));
end
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
% ylim([-0.4 0.8])
% yticks([-0.4 0 0.4 0.8])
axis square;
fnm = ['ThetaIndex_IN_vs_Refractory' cb '.fig'];
saveas(H10,fullfile(resdir,fnm))   % save plot
fnm = ['ThetaIndex_IN_vs_Refractory' cb '.jpeg'];
saveas(H10,fullfile(resdir,fnm))   % save plot
close(H10);


%% REVIEWER FIGURE 2 PANEL E: ThetaIndex_iN 
% MEDIAN

coords = [0.5, 2, 4];

% Significance test
[groups, pAll] = sigTest(ThetaIndex_IN, phasicB, poissonL, tonic, coords);

% x coordinate generation to scatter the plot
x1 = linspace(coords(1)-0.4, coords(1)+0.4, sum(phasicB));
x2 = linspace(coords(2)-0.4, coords(2)+0.4, sum(poissonL));
x3 = linspace(coords(3)-0.4, coords(3)+0.4, sum(tonic));

H4 = figure;
hold on;
if sum(phasicB)>0
    plot(x1, ThetaIndex_IN(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
end
plot(x2, ThetaIndex_IN(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
if sum(tonic)>0
    plot(x3, ThetaIndex_IN(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
end
bar(coords, [nanmedian(ThetaIndex_IN(phasicB)) nanmedian(ThetaIndex_IN(poissonL))...
    nanmedian(ThetaIndex_IN(tonic))],0.8, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index');
xlabel('Groups');
xticks([])
title('   ')
% y_lim = [-0.5 1];
% ylim(y_lim);
% yticks([-0.4 0 0.4 0.8])
sigstar(groups,pAll);
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
axis square;
saveas(H4, [RESDIR cb '\tonicphasic' cb PATH '\ThetaIndex_IN_median.fig']);
saveas(H4, [RESDIR cb '\tonicphasic' cb PATH '\ThetaIndex_IN_median.jpeg']);
close(H4);

%% FIGURE PANEL K: Plot scatter for Tonic - Phasic groups TI vs Ref
% Regression

if ~strcmp(cb, 'HDB')
    y = ThetaIndex(tonic)';
    x = Refractory(tonic)';
    [b,~,~,~,stats] = regress(y',[ones(length(x),1),x']);
    p = round(stats(3),4);           % F-test significance
end
H10 = figure;
hold on
if ~strcmp(cb, 'HDB')
    plot(x,y,'o','MarkerSize',12,'MarkerFaceColor','none',...
        'MarkerEdgeColor',cCode(3,:))
    plot(x(Tindex),y(Tindex),'o','MarkerSize',12,'MarkerFaceColor',...
        'none','MarkerEdgeColor',cCode(3,:), 'LineWidth', 2)  % scatter plot
    icp = b(1);   % intercept
    gr = b(2);   % gradient
    xx = 0:0.01:max(x);
    yy = xx .* gr + icp;
end
axis square
hold on
if ~strcmp(cb, 'HDB')
    plot(xx,yy, '--', 'Color',cCode(3,:),'LineWidth',2)   % overlay regression line
end
plot(Refractory(phasicB),ThetaIndex(phasicB),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:))  % scatter plot
plot(Refractory(poissonL),ThetaIndex(poissonL),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:))  % scatter plot
plot(Refractory(Pindex),ThetaIndex(Pindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:), 'LineWidth', 2)  % scatter plot
plot(Refractory(Bindex),ThetaIndex(Bindex),'o','MarkerSize',12,...
    'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:), 'LineWidth', 2)  % scatter plot
title('   ')
xlabel('Refractory (ms)')
ylabel('Theta Index')
x_lim = xlim;
y_lim = ylim;
if ~strcmp(cb, 'HDB')
    text((x_lim(2)*0.45), y_lim(2)*0.95, ['p: ' num2str(p)], 'Color',...
        cCode(3,:));
end
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
ylim([-0.4 0.8])
yticks([-0.4 0 0.4 0.8])
axis square;
fnm = ['ThetaIndex_vs_Refractory' cb '.fig'];
saveas(H10,fullfile(resdir,fnm))   % save plot
fnm = ['ThetaIndex_vs_Refractory' cb '.jpeg'];
saveas(H10,fullfile(resdir,fnm))   % save plot
close(H10);

%% REVIEWER FIGURE2 PANEL B: ThetaIndex barplot with overlaying individual points
% MEDIAN

coords = [0.5, 2, 4];

% Significance test
[groups, pAll] = sigTest(ThetaIndex, phasicB, poissonL, tonic, coords);

% x coordinate generation to scatter the plot
x1 = linspace(coords(1)-0.4, coords(1)+0.4, sum(phasicB));
x2 = linspace(coords(2)-0.4, coords(2)+0.4, sum(poissonL));
x3 = linspace(coords(3)-0.4, coords(3)+0.4, sum(tonic));

H4 = figure;
hold on;
if sum(phasicB)>0
    plot(x1, ThetaIndex(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
end
plot(x2, ThetaIndex(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
if sum(tonic)>0
    plot(x3, ThetaIndex(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
end
bar(coords, [nanmedian(ThetaIndex(phasicB)) nanmedian(ThetaIndex(poissonL))...
    nanmedian(ThetaIndex(tonic))],0.8, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index');
xlabel('Groups');
xticks([])
title('   ')
y_lim = [-0.4 0.8];
ylim(y_lim);
yticks([-0.4 0 0.4 0.8])
sigstar(groups,pAll);
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
axis square;
saveas(H4, [RESDIR cb '\tonicphasic' cb PATH '\ThetaIndex_median.fig']);
saveas(H4, [RESDIR cb '\tonicphasic' cb PATH '\ThetaIndex_median.jpeg']);
close(H4);


%% FIGURE 3 PANEL F: FiringRate

% Significance test
[groups, pAll] = sigTest(FR, phasicB, poissonL, tonic, coords);

phasicB_se = nanse(FR(phasicB));
poisson_se = nanse(FR(poissonL));
tonic_se = nanse(FR(tonic));

phasicB_sem = se_of_median(FR(phasicB));
poisson_sem = se_of_median(FR(poissonL));
tonic_sem = se_of_median(FR(tonic));


% MEDIAN
H1 = figure;
hold on;
bar(coords, [nanmedian(FR(phasicB)) nanmedian(FR(poissonL))...
    nanmedian(FR(tonic))], 'FaceColor', 'none', 'LineWidth',3)
plot(x1, FR(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
plot(x2, FR(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
plot(x3, FR(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
sigstar(groups,pAll);
ylabel('Firing rate (Hz)');
xlabel('Groups')
xticks([])
y_lim = ylim;
y_lim = ceil(y_lim/10)*10;
ylim([y_lim])
yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
axis square;
annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
setmyplot_tamas;
saveas(H1, [RESDIR cb '\tonicphasic' cb PATH '\BaselineFR_median.fig']);
saveas(H1, [RESDIR cb '\tonicphasic' cb PATH '\BaselineFR_median.jpeg']);
close(H1);





%--------------------------------------------------------------------------
function [Bindex, Pindex, Tindex] = exampleACGs(cb, cellids,...
    Refractory, BurstIndex)
% Example cellids

if strcmp(cb, 'NB') | strcmp(cb, 'POOLED')
    choosecb('NB');
    Bexp='n028_120211a_8.1';
    Pexp='n046_130104a_6.2';
    Texp='n046_130101a_6.1';
    Tindex = contains(cellids,Texp);
    Bindex = contains(cellids,Bexp);
    Pindex= contains(cellids,Pexp);
    
elseif strcmp(cb, 'HDB')
    choosecb('HDB');
    Bexp='n078_150104a_1.1';
    Pexp='n067_141009b_2.2';
    Bcells = strfind(cellids,Bexp);
    BPos = find(cellfun(@isempty,Bcells));
    Bindex = ones(1, length(Bcells));
    Bindex(BPos)=0;
    Bindex = logical(Bindex);
    Pcells = strfind(cellids,Pexp);
    PPos = find(cellfun(@isempty,Pcells));
    Pindex = ones(1, length(Pcells));
    Pindex(PPos)=0;
    Pindex = logical(Pindex);
    Tindex = zeros(1,length(cellids));
    Tindex = logical(Tindex);
else
    [Bindex, Pindex, Tindex] = deal(logical(zeros(1, length(cellids))));
end

%--------------------------------------------------------------------------
function [groups, pAll] = sigTest(currData, group1, group2, group3, coords)

% Significance test
% Boxplots
g1Num = sum(group1);
g2Num = sum(group2);
g3Num = sum(group3);

p1 = [];
p2 = [];
p3 = [];
if g1Num>0 & g2Num>0
    [H1, p1] = boxstat(currData(group1),currData(group2),'','',0.05,'nonpaired');
    close(H1);
    p(1) = round(p1,3);
end
if g1Num>0 & g3Num>0
    [H2, p2] = boxstat(currData(group1),currData(group3),'','',0.05,'nonpaired');
    close(H2);
    p(2) = round(p2,3);
end
if g2Num>0 & g3Num>0
    [H3, p3] = boxstat(currData(group2),currData(group3),'','',0.05,'nonpaired');
    close(H3);
    p(3) = round(p3,3);
end

groupInc = [p1<=0.05, p2<=0.05, p3<0.05];

groups={[coords(1),coords(2)],[coords(1),coords(3)],...
    [coords(2),coords(3)]};
groups = groups(groupInc);
pAll = p(groupInc);

%--------------------------------------------------------------------------
function [mn_tonic, se_tonic, numT] = nbUntagged(lags, win)

% Load untagged acg matrix
load('D:\_MATLAB_DATA\NB\acgNB\other\ACG_matrices_NB.mat')

tonic = Refractory >= 40;   % refractory above 40 ms

tCCR = CCR(tonic,:);
% Normalize ACG
CCRnorm = nan(size(tCCR));
numT = size(tCCR,1);
for k = 1:numT
    sccr = smooth(tCCR(k,:),'linear',5);   % smooth
    b = mean(sccr(lags>-500&lags<500));
    CCRnorm(k,:) = sccr / b;
end

tonicCCR = CCRnorm;  % ACG for tonic cells
mn_tonic = mean(tonicCCR(:,lags>win(1)&lags<win(2)));   % average ACG, tonic cells
se_tonic = std(tonicCCR(:,lags>win(1)&lags<win(2))) /...
    sqrt(size(tonicCCR,1));   % SE, tonic cells
