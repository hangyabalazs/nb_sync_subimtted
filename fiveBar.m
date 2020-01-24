function fiveBar(varargin)
%FIVEBAR   Stacked barplot for phasicB cells in theta and gamma.
%   FIVEBAR(BARDATA) Creates barplots for input data from groupStaErs,
%   or for barData matrix.
%
%   See also GROUPSTAERS

%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

choosecb('NB');

% Input arguments
prs = inputParser;
addParameter(prs,'barData',[],@(s)isempty(s)|isstruct(s)) % barData struct
parse(prs,varargin{:})
g = prs.Results;
dbstop if error;
cb = whichcb;
fs = filesep;
numVer = '20';

global RESDIR;
global cCode;
resdir = [RESDIR cb fs 'sta' cb fs '_ALLDATA'];
if isempty(g.barData)
    load([resdir fs 'barData20.mat'])
else
    barData = g.barData;
end

sessType = fieldnames(barData);

% THETA
allB_T = barData.(sessType{1}).none.original.Bursting(2,:); % Raw
burstB_T = barData.(sessType{1}).burst1.original.Bursting(2,:); % Burst1
singleB_T = barData.(sessType{1}).single.original.Bursting(2,:); % Single
syncB_T = barData.(sessType{1}).burst1Single.sync.Bursting(2,:); % Sync
nonsyncB_T = barData.(sessType{1}).burst1Single.nonsync.Bursting(2,:); % Nonsync

% FiveBar for Theta: All, burst1, single, sync, nonsync
burstT = [mean(allB_T), mean(burstB_T), mean(singleB_T),...
    mean(syncB_T), mean(nonsyncB_T)];

% SE: Theta for phasicB
SE_all = std(allB_T) / sqrt(numel(allB_T));
SE_BT = std(burstB_T) / sqrt(numel(burstB_T));
SE_ST = std(singleB_T) / sqrt(numel(singleB_T));
SE_syT = std(syncB_T) / sqrt(numel(syncB_T));
SE_asyT = std(nonsyncB_T) / sqrt(numel(nonsyncB_T));
SE_T = [SE_all SE_BT SE_ST SE_syT SE_asyT];

% Boxplots
[H1, p1] = boxstat(burstB_T,singleB_T,'Bursts','Single spikes',0.05,'paired');
close(H1);
[H2, p2] = boxstat(syncB_T,nonsyncB_T,'Sync','Async',0.05,'paired');
close(H2);

pT(1) = round(p1,3);
pT(2) = round(p2,3);

groupInc = [p1<=0.05, p2<=0.05];

groups={[2,3],[4.5,5.5]};
groups = groups(groupInc);
pTAll = pT(groupInc);

% Theta fiveBar plot
tickPos = [1 2 3 4.5 5.5]; % Bar x positions

xmin=-0.25;
xmax=0.25;
% x coordinate generation to scatter the plot
x1 = xmin+rand(numel(allB_T),1)*(xmax-xmin)+1;
x2 = xmin+rand(numel(burstB_T),1)*(xmax-xmin)+2;
x3 = xmin+rand(numel(singleB_T),1)*(xmax-xmin)+3;
x4 = xmin+rand(numel(syncB_T),1)*(xmax-xmin)+4.5;
x5 = xmin+rand(numel(nonsyncB_T),1)*(xmax-xmin)+5.5;


H3 = figure;
hold on;
plot(x1, allB_T, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth',3)
plot(x2, burstB_T, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 0 1], 'LineWidth',3)
plot(x3, singleB_T, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(3,:), 'LineWidth',3)
plot(x4, syncB_T, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(1,:), 'LineWidth',3)
plot(x5, nonsyncB_T, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0.6], 'LineWidth',3)
TG = bar(tickPos, burstT, 'FaceColor','none', 'LineWidth',2);
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 0 0];
TG.CData(2,:) = [1 0 1];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
TG.CData(5,:) = [0 0 0.6];
sigstar(groups,pTAll);
axis square;
ylabel('Theta Power (dB)');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 6])
ylim([-0.2 0.8])
yticks([0 0.4 0.8])
title('PhasicB Theta')
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fiveBar_new' num2str(numVer) '.fig']);
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fiveBar_new' num2str(numVer) '.jpeg']);
close(H3);

% GAMMA
allB_G = barData.(sessType{1}).none.original.Bursting(3,:); % Raw
burstB_G = barData.(sessType{1}).burst1.original.Bursting(3,:); % Burst1
singleB_G = barData.(sessType{1}).single.original.Bursting(3,:); % Single
syncB_G = barData.(sessType{1}).burst1Single.sync.Bursting(3,:); % Sync
nonsyncB_G = barData.(sessType{1}).burst1Single.nonsync.Bursting(3,:); % Nonsync

% FiveBar for Gamma: All, burst1, single, sync, nonsync
burstG = [mean(allB_G), mean(burstB_G), mean(singleB_G),...
    mean(syncB_G), mean(nonsyncB_G)];

% SE: Gamma for phasicB
SE_allG = std(allB_G) / sqrt(numel(allB_G));
SE_BG = std(burstB_G) / sqrt(numel(burstB_G));
SE_SG = std(singleB_G) / sqrt(numel(singleB_G));
SE_syG = std(syncB_G) / sqrt(numel(syncB_G));
SE_asyG = std(nonsyncB_G) / sqrt(numel(nonsyncB_G));
SE_G = [SE_allG SE_BG SE_SG SE_syG SE_asyG];

% Boxplots
[H4, p3] = boxstat(burstB_G,singleB_G,'Bursts','Single spikes',0.05,'paired');
close(H4);
[H5, p4] = boxstat(syncB_G,nonsyncB_G,'Sync','Async',0.05,'paired');
close(H5);

pG(1) = round(p3,3);
pG(2) = round(p4,3);

groupInc = [p3<=0.05, p4<=0.05];

groups={[2,3],[4.5,5.5]};
groups = groups(groupInc);
pGAll = pG(groupInc);


% Gamma fiveBar plot
H6 = figure;
hold on;
plot(x1, allB_G, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth',3)
plot(x2, burstB_G, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [1 0 1], 'LineWidth',3)
plot(x3, singleB_G, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(3,:), 'LineWidth',3)
plot(x4, syncB_G, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(1,:), 'LineWidth',3)
plot(x5, nonsyncB_G, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0.6], 'LineWidth',3)
TG = bar(tickPos, burstG, 'FaceColor','none', 'LineWidth',2);
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 0 0];
TG.CData(2,:) = [1 0 1];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
TG.CData(5,:) = [0 0 0.6];
sigstar(groups,pGAll);
axis square;
ylabel('Gamma Power (dB)');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 6])
ylim([-0.2 0.8])
yticks([0 0.4 0.8])
title('PhasicB Gamma')
saveas(H6, [resdir '\ExtraPlots\PhasicB_Gamma_fiveBar_new' num2str(numVer) '.fig']);
saveas(H6, [resdir '\ExtraPlots\PhasicB_Gamma_fiveBar_new' num2str(numVer) '.jpeg']);
close(H6);