function parameterPlot(cellbase)
%PARAMETERPLOT  Generates baseline, spikeshape, and thetaindex
% figures for the cholinergic cells in the current cellbase.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

dbstop if error;
global RESDIR;
global PATH;
global cCode;

% ACG matrix loader
acgDir = [RESDIR cellbase '\acg' cellbase PATH '\'];
acgPath = [acgDir 'ACG_matrices_' cellbase '.mat'];
load(acgPath);

% ChAT, pChAT properties
for iC = 1:numel(cellids)
    cbSwitcher(cellids{iC});
    currCB = whichcb;
    if ~strcmp(currCB, 'PannaHDB')
        ChAT(iC) = getvalue('ChAT+',cellids(iC));
    else
        ChAT(iC) = 1;
    end
    pChAT(iC) = ~ChAT(iC);
end
ch=ChAT==1; % tagged
pch=pChAT==1; % putative

% TP groups
phasicB = groupID == 1;
poissonL = groupID == 2;
tonic = groupID == 3;
BFCNburst = groupID==1 | groupID==2;
BFCNreg = groupID==3;

wnSize = num2str(lags(end));

%% Baseline
Baseline = getvalue('Baseline', cellids);

%Baseline bar plot
if ~isnan(Baseline)
    if ~isempty(Baseline)
        % MEAN
        H1 = figure;
        hold on;
        bar([1,2,3.5], [nanmean(Baseline(phasicB)) nanmean(Baseline(poissonL)) nanmean(Baseline(tonic))], 'FaceColor', [0.7 0.7 0.7], 'LineWidth',2)
        plot(1, Baseline(phasicB), '*', 'Color', cCode(1,:), 'LineWidth',3)
        plot(2, Baseline(poissonL), '*', 'Color', cCode(2,:), 'LineWidth',3)
        plot(3.5, Baseline(tonic), '*', 'Color', cCode(3,:), 'LineWidth',3)
        ylabel('Baseline firing rate (Hz)');
        xlabel('Groups')
        xticks([])
        y_lim = ylim;
        yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
        axis square;
        setmyplot_tamas;
        saveas(H1, [RESDIR cellbase '\tonicphasic' cellbase PATH '\BaselineFR_mean.fig']);
        saveas(H1, [RESDIR cellbase '\tonicphasic' cellbase PATH '\BaselineFR_mean.jpeg']);
        close(H1);
        % MEDIAN
        H1 = figure;
        hold on;
        bar([1,2,3.5], [nanmedian(Baseline(phasicB)) nanmedian(Baseline(poissonL)) nanmedian(Baseline(tonic))], 'FaceColor', [0.7 0.7 0.7], 'LineWidth',2)
        plot(1, Baseline(phasicB), '*', 'Color', cCode(1,:), 'LineWidth',3)
        plot(2, Baseline(poissonL), '*', 'Color', cCode(2,:), 'LineWidth',3)
        plot(3.5, Baseline(tonic), '*', 'Color', cCode(3,:), 'LineWidth',3)
        ylabel('Baseline firing rate (Hz)');
        xlabel('Groups')
        xticks([])
        y_lim = ylim;
        yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
        axis square;
        setmyplot_tamas;
        saveas(H1, [RESDIR cellbase '\tonicphasic' cellbase PATH '\BaselineFR_median.fig']);
        saveas(H1, [RESDIR cellbase '\tonicphasic' cellbase PATH '\BaselineFR_median.jpeg']);
        close(H1);
        % Stats, boxplots
        [H1 Wp1] = boxstat(Baseline(phasicB),Baseline(poissonL),'Bursting','PoissonL',0.05,'nonpaired');
        [H2 Wp2] = boxstat(Baseline(phasicB),Baseline(tonic),'Bursting','Tonic',0.05,'nonpaired');
        [H3 Wp3] = boxstat(Baseline(tonic),Baseline(poissonL),'Tonic','PoissonL',0.05,'nonpaired');
        close all;
    end
end


%% SpikeShape
if ~strcmp(cellbase, 'POOLED')
    for i = 1:length(cellids)
        Spike = getvalue('SpikeShape', cellids(i));
        spikeShapes(i,:) = Spike{1}.Spike;
    end
    
    
    H2 = figure;
    hold on;
    maxDiff = [max(mean(spikeShapes(phasicB,:),1)) max(mean(spikeShapes(poissonL,:),1)) max(mean(spikeShapes(tonic,:),1))];
    shapeB = mean(spikeShapes(phasicB,:),1)/maxDiff(1);
    shapeP = mean(spikeShapes(poissonL,:),1)/maxDiff(2);
    shapeT = mean(spikeShapes(tonic,:),1)/maxDiff(3);
    plot(shapeB,'Color', cCode(1,:), 'LineWidth', 3)
    plot(shapeP,'Color', cCode(2,:), 'LineWidth', 3)
    plot(shapeT,'Color', cCode(3,:), 'LineWidth', 3)
    ylabel('Baseline firing rate (Hz)');
    axis off;
    setmyplot_tamas;
    saveas(H2, [RESDIR cellbase '\tonicphasic' cellbase PATH '\AverageSpikeShape.fig']);
    saveas(H2, [RESDIR cellbase '\tonicphasic' cellbase PATH '\AverageSpikeShape.jpeg']);
    close(H2);
    
    % SpikeShape for Chat vs pChat
    H22 = figure;
    hold on;
    maxDiff = [max(mean(spikeShapes(ch,:),1)) max(mean(spikeShapes(pch,:),1))];
    shapeCh = mean(spikeShapes(ch,:),1)/maxDiff(1);
    shapePCh = mean(spikeShapes(pch,:),1)/maxDiff(2);
    plot(shapeCh,'Color', cCode(1,:), 'LineWidth', 3)
    plot(shapePCh,'Color', cCode(2,:), 'LineWidth', 3)
    ylabel('Baseline firing rate (Hz)');
    axis off;
    setmyplot_tamas;
    saveas(H22, [RESDIR cellbase '\tonicphasic' cellbase PATH '\AverageSpikeShape_chatPchat.fig']);
    saveas(H22, [RESDIR cellbase '\tonicphasic' cellbase PATH '\AverageSpikeShape_chatPchat.jpeg']);
    close(H22);
    
    % Burst vs regular firing
    H23 = figure;
    hold on;
    maxDiff = [max(mean(spikeShapes(BFCNburst,:),1)) max(mean(spikeShapes(BFCNreg,:),1))];
    shapeBurst = mean(spikeShapes(BFCNburst,:),1)/maxDiff(1);
    shapeReg = mean(spikeShapes(BFCNreg,:),1)/maxDiff(2);
    plot(shapeBurst,'Color', cCode(1,:), 'LineWidth', 3)
    plot(shapeReg,'Color', cCode(3,:), 'LineWidth', 3)
    ylabel('Baseline firing rate (Hz)');
    axis off;
    setmyplot_tamas;
    saveas(H23, [RESDIR cellbase '\tonicphasic' cellbase PATH '\AverageSpikeShape_ef_lf.fig']);
    saveas(H23, [RESDIR cellbase '\tonicphasic' cellbase PATH '\AverageSpikeShape_ef_lf.jpeg']);
    close(H23);
end
%% ThetaIndex
thetaIndex = ThetaIndex;
thetaIndexN = ThetaIndexN;

coords = [1,2,3.5];

% Significance test
[groups, pAll] = sigTest(thetaIndex, phasicB, poissonL, tonic, coords);


% x coordinate generation to scatter the plot
x1 = rand(sum(phasicB),1)+0.5;
x2 = rand(sum(poissonL),1)+1.5;
x3 = rand(sum(tonic),1)+3;

% ThetaIndex barplot with scattered individual points (MEAN and MEDIAN)
% MEDIAN
H3 = figure;
hold on;
plot(x1, thetaIndex(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
plot(x2, thetaIndex(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
plot(x3, thetaIndex(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
bar(coords, [nanmedian(thetaIndex(phasicB)) nanmedian(thetaIndex(poissonL)) nanmedian(thetaIndex(tonic))],1, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index');
xlabel('Groups');
title('ThetaIndex median')
xticks([])
y_lim = [-0.2 1];
ylim(y_lim);
yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
sigstar(groups,pAll);
setmyplot_tamas;
saveas(H3, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_scatter_median' wnSize '.fig']);
saveas(H3, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_scatter_median' wnSize '.jpeg']);
close(H3);

% MEAN
H3 = figure;
hold on;
plot(x1, thetaIndex(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
plot(x2, thetaIndex(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
plot(x3, thetaIndex(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
bar(coords, [nanmean(thetaIndex(phasicB)) nanmean(thetaIndex(poissonL)) nanmean(thetaIndex(tonic))],1, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index');
xlabel('Groups');
title('ThetaIndex mean')
xticks([])
y_lim = [-0.1 1];
ylim(y_lim);
yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
sigstar(groups,pAll);
setmyplot_tamas;
saveas(H3, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_scatter_mean' wnSize '.fig']);
saveas(H3, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_scatter_mean' wnSize '.jpeg']);
close(H3);

% ThetaIndex barplot with overlaying individual points
% MEDIAN
H4 = figure;
hold on;
if sum(phasicB)>0
    plot(coords(1), thetaIndex(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
end
plot(coords(2), thetaIndex(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
if sum(tonic)>0
    plot(coords(3), thetaIndex(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
end
bar(coords, [nanmedian(thetaIndex(phasicB)) nanmedian(thetaIndex(poissonL)) nanmedian(thetaIndex(tonic))],1, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index');
xlabel('Groups');
xticks([])
title('ThetaIndex median')
y_lim = [-0.1 1];
ylim(y_lim);
yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
sigstar(groups,pAll);
setmyplot_tamas;
saveas(H4, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_median' wnSize '.fig']);
saveas(H4, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_median' wnSize '.jpeg']);
close(H4);

% MEAN
H4 = figure;
hold on;
if sum(phasicB)>0
    plot(coords(1), thetaIndex(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
end
plot(coords(2), thetaIndex(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
if sum(tonic)>0
    plot(coords(3), thetaIndex(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
end
bar(coords, [nanmean(thetaIndex(phasicB)) nanmean(thetaIndex(poissonL)) nanmean(thetaIndex(tonic))],1, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index');
xlabel('Groups');
xticks([])
title('ThetaIndex mean')
y_lim = [-0.1 1];
ylim(y_lim);
yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
sigstar(groups,pAll);
setmyplot_tamas;
saveas(H4, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_mean' wnSize '.fig']);
saveas(H4, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndex_mean' wnSize '.jpeg']);
close(H4);

%% ThetaIndex normalized with ACGmax/FR
% MEDIAN
[groups, pAll] = sigTest(thetaIndexN, phasicB, poissonL, tonic, coords);

% ThetaIndex barplot with scattered individual points (MEAN and MEDIAN)
% MEDIAN
H3 = figure;
hold on;
plot(x1, thetaIndexN(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
plot(x2, thetaIndexN(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
plot(x3, thetaIndexN(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
bar(coords, [nanmedian(thetaIndexN(phasicB)) nanmedian(thetaIndexN(poissonL)) nanmedian(thetaIndexN(tonic))],1, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index normalized');
xlabel('Groups');
title('ThetaIndex median normalized')
xticks([])
y_lim = [-1 1];
ylim(y_lim);
yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
sigstar(groups,pAll);
setmyplot_tamas;
saveas(H3, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndexN_scatter_median' wnSize '.fig']);
saveas(H3, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndexN_scatter_median' wnSize '.jpeg']);
close(H3);


H4 = figure;
hold on;
if sum(phasicB)>0
    plot(coords(1), thetaIndexN(phasicB), 'o', 'Color', cCode(1,:), 'LineWidth',3)
end
plot(coords(2), thetaIndexN(poissonL), 'o', 'Color', cCode(2,:), 'LineWidth',3)
if sum(tonic)>0
    plot(coords(3), thetaIndexN(tonic), 'o', 'Color', cCode(3,:), 'LineWidth',3)
end
bar(coords, [nanmedian(thetaIndexN(phasicB)) nanmedian(thetaIndexN(poissonL)) nanmedian(thetaIndexN(tonic))],1, 'FaceColor', 'none', 'LineWidth',3)
ylabel('Theta Index Normalized');
xlabel('Groups');
xticks([])
title('ThetaIndex median normalized')
y_lim = [-1 1];
ylim(y_lim);
yticks([y_lim(1) y_lim(2)/2 y_lim(2)])
sigstar(groups,pAll);
setmyplot_tamas;
saveas(H4, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndexN_median' wnSize '.fig']);
saveas(H4, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndexN_median' wnSize '.jpeg']);
close(H4);


%% Plot scatter for Tonic - Phasic groups TI vs Ref
% Regression
if ~strcmp(cellbase, 'HDB')
    y = ThetaIndexN(tonic);
    x = Refractory(tonic)';
    [b,~,~,~,stats] = regress(y',[ones(length(x),1),x']);
    p = round(stats(3),4);           % F-test significance
end
H10 = figure;
hold on
if ~strcmp(cellbase, 'HDB')
    plot(x,y,'o','MarkerSize',12,'MarkerFaceColor',cCode(3,:),'MarkerEdgeColor','k')
    icp = b(1);   % intercept
    gr = b(2);   % gradient
    xx = 0:0.01:max(x);
    yy = xx .* gr + icp;
end
axis square
hold on
if ~strcmp(cellbase, 'HDB')
    plot(xx,yy, '--', 'Color',cCode(3,:),'LineWidth',2)   % overlay regression line
end
plot(Refractory(phasicB),ThetaIndexN(phasicB),'o','MarkerSize',12,'MarkerFaceColor',cCode(1,:),'MarkerEdgeColor','k')  % scatter plot
plot(Refractory(poissonL),ThetaIndexN(poissonL),'o','MarkerSize',12,'MarkerFaceColor',cCode(2,:),'MarkerEdgeColor','k')  % scatter plot
title('ThetaIndex Normalized vs Refractory')
xlabel('Refractory')
ylabel('ThetaIndex normalized')
x_lim = xlim;
y_lim = ylim;
if ~strcmp(cellbase, 'HDB')
    text((x_lim(2)*0.45), y_lim(2)*0.95, ['p: ' num2str(p)], 'Color', cCode(3,:));
end
axis square;
setmyplot_tamas;
saveas(H10, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndexN_scatterNew' wnSize '.fig']);
saveas(H10, [RESDIR cellbase '\tonicphasic' cellbase PATH '\ThetaIndexN_scatterNew' wnSize '.jpeg']);
close(H10);


%%
%---------------------------------------------------------------------------------------------------------
function [groups, pAll] = sigTest(currData, group1, group2, group3, coords)


% Significance test
% Boxplots
[H1, p1] = boxstat(currData(group1),currData(group2),'','',0.05,'nonpaired');
close(H1);
[H2, p2] = boxstat(currData(group1),currData(group3),'','',0.05,'nonpaired');
close(H2);
[H3, p3] = boxstat(currData(group2),currData(group3),'','',0.05,'nonpaired');
close(H3);

p(1) = round(p1,3);
p(2) = round(p2,3);
p(3) = round(p3,3);

groupInc = [p1<=0.05, p2<=0.05, p3<0.05];

groups={[coords(1),coords(2)],[coords(1),coords(3)],...
    [coords(2),coords(3)]};
groups = groups(groupInc);
pAll = p(groupInc);
