function chatProbabilityDist

% Load auto-correlations
% Select corresponding ACG matrix
cellbase = 'NB';
global RESDIR;
global PATH;
global cCode;
global plotWN;

%% Prepare
% Global variables
cb = 'NB'; % select cellbase for nucleus basalis
initGlobals(cb); % initialize global variables
pathSelect(cb); % select path for saving

% Select cholinergic cells
if strcmp(cellbase, 'POOLED')
    [ChAT, pChAT, allChAT] = prepareScript;
else
    [ChAT, pChAT, allChAT] = selectChAT(cb);
end
NumChAT = length(allChAT);

binsize = 40;

% Load ACG file
fnm_acg = [RESDIR cellbase '\acg' cellbase PATH '\ACG_matrices_' cellbase '.mat'];
load(fnm_acg);   % Load ACG matrix

[cInx, pInx] = deal(zeros(numel(Refractory),1));

regInx = Refractory >= 40;
regCells = cellids(regInx);

ChATInx = find(contains(cellids, ChAT));
pChATInx = find(contains(cellids, pChAT));

cInx(ChATInx) = 1;
pInx(pChATInx) = 1;

cRef = (Refractory >= 40);

c1 = cRef + cInx;
c2 = cRef + pInx;

ChatRef = sort(Refractory(c1==2));
pChatRef = sort(Refractory(c2==2));
allChATRef = sort(Refractory(cRef));

% UNTAGGED
fnm_acg = [RESDIR  'NB\acgNB' PATH '\other\ACG_matrices_NB.mat'];
load(fnm_acg);   % Load ACG matrix

regUntaggedInx = groupID ==3;
regUntCells = cellids(regUntaggedInx);
untRef = sort(Refractory(regUntaggedInx));

pooledCells = sort([allChATRef; untRef]);
switch binsize
    case 10
        bins = 40:20:160;

    case 20
        bins = 40:20:160;

    case 40
        bins = [40, 80; 60, 100; 80, 120; 100, 140; 120, 180];

end
% bins = 40:20:160;
for iB = 1:size(bins,1)
        inRangeChAT = sum((allChATRef>= bins(iB,1)) & (allChATRef < bins(iB,2)));
        inRangeAll = sum((pooledCells>= bins(iB,1)) & (pooledCells < bins(iB,2)));
        binNumChAT(iB) = inRangeChAT;
        binNumAll(iB) = inRangeAll;
end
Ratio = binNumChAT ./ binNumAll;
allRatio = round(mean(Ratio),2);
% H1 = figure;
% hold on;
% bar([binNumChAT; binNumAll], 'LineWidth',2)
% bar(1:length(binNumChAT), [binNumAll], 'FaceColor', 'none','EdgeColor',[0.8 0 0], 'LineWidth',2)

%%
H1 = figure;
hold on;
title([cellbase ' ChAT probability for Regular cells'])
yyaxis left % Left axis for log(Refractory) hist values
bar(binNumAll','FaceColor',[0.6 0.6 0.6],'EdgeColor',cCode(3,:),'LineWidth',2)
bar(binNumChAT','FaceColor',cCode(3,:),'EdgeColor',cCode(3,:),'LineWidth',2)
ylabel('Number of cells')
ylim([0 max(binNumAll)+2])
xticks([1:5])
switch binsize
    case 10
        xticklabels({'40-50', '50-60', '60-70', '70-80', '80-90', '90-100', '100-110',...
            '110-120', '120-130', '130-140', '140-150'})
    case 20
        xticklabels({'40-60', '60-80', '80-100', '100-120', '120-140', '140-160'})
    case 40
        xticklabels({'40-80', '60-100', '80-120', '100-140', '120-180'})
end
xlabel('Refractory bins (ms)')



yyaxis right
plot(Ratio, 'o-', 'Color', [0,0,0], 'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[0,0,0],'LineWidth',2)
ylabel('Ratio')
ylim([0 1])
legend({'All (Untagged+ChAT+pChAT', 'ChAT+pChAT',...
    ['ChAT probability, Average: ', num2str(allRatio)]})
legend('boxoff');
% setmyplot_tamas;
axis square;
saveas(H1, ['D:\_MATLAB_DATA\NB\ChATprobabilityReg\' cellbase '_chatProbStacked.fig']);
saveas(H1, ['D:\_MATLAB_DATA\NB\ChATprobabilityReg\' cellbase '_chatProbStacked.jpeg']);
close(H1);



%% HDB

choosecb('HDB'); % Change to HDB cellbase
cb = whichcb;
[~, ~, allChAT] = selectChAT(cb);

[allChAT, ~] = selectPCells(cb, allChAT); % selectting putative cholinergic cells


%% UNTAGGED
choosecb('NB'); % Change to NB cellbase
[~, ~, allChAT] = selectChAT('other');
pathSelect('other');
currCB = whichcb;
nbtonicphasic_NBsync(currCB);   % ACG Tonic-Phasic properties
parameterPlot;
