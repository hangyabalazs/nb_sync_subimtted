function ccgCorrBars(varargin)
%CCGCORRBARS   Bar plots for CCG ratios
%   CCGCORRBARS Loads CCG matrices for Burst1 and Single. Calculates and plots  
%   average cross-correlations for TP groups.
%
%   See also CCG, CELLID2TPGROUP.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

dbstop if error;
if nargin>0
    cb = varargin{1};
else
    cb = whichcb;
end
global RESDIR;
global PATH;
fs = filesep;

resdir = [RESDIR cb fs 'ccgCorrBars' cb fs];

sessType = 'behavior';
Burst1 = load([RESDIR cb fs 'ccg' cb fs sessType fs 'burst1' cb PATH fs 'CCG_matrices_'...
    sessType '_burst1_' cb '.mat']);
Single = load([RESDIR cb fs 'ccg' cb fs sessType fs 'single' cb PATH fs 'CCG_matrices_'...
    sessType '_single_' cb '.mat']);
mode = {'Burst1', 'Single'};

plotWin = [-250 250]; % +/- 250ms lag
lags = [plotWin(1):1:-1, 1:1:plotWin(2)]; % +/- 1 s window, 1 ms sampling rate

% Loading Data into struct, define pair identities
for iM = 1:length(mode)
    currmode = mode{iM}; % Burst1 or Single
    currData = eval(currmode);
    PairOfCells = currData.PairOfCells;
    group1 = currData.group1;
    group2 = currData.group2;
    
    groupID = zeros(1,size(PairOfCells,1));
    yLabel=cell(1, size(PairOfCells,1));
    
    % Categorize Pairs to TP groups
    for numPair = 1:size(PairOfCells,1)
        currgroup1 = group1{numPair};
        currgroup2 = group2{numPair};
        if strcmp(currgroup1, currgroup2) % In case of the pairs group ID is the same
            switch currgroup1{1}
                case 'phasicB'
                    groupID(numPair) = 1;
                    yLabel{1,numPair} = 'B';
                case 'poissonL'
                    groupID(numPair) = 2;
                    yLabel{1,numPair} = 'PL';
                case 'tonic'
                    groupID(numPair) = 3;
                    yLabel{1,numPair} = 'T';
            end
        else % Mixed pair
            if strcmp(currgroup1, 'phasicB') && strcmp(currgroup2, 'poissonL') || strcmp(currgroup1, 'poissonL') && strcmp(currgroup2, 'phasicB')
                groupID(numPair) = 4;
                yLabel{1,numPair} = 'B PL';
            elseif strcmp(currgroup1, 'phasicB') && strcmp(currgroup2, 'tonic') || strcmp(currgroup1, 'tonic') && strcmp(currgroup2, 'phasicB')
                groupID(numPair) = 5;
                yLabel{1,numPair} = 'B T';
            elseif strcmp(currgroup1, 'tonic') && strcmp(currgroup2, 'poissonL') || strcmp(currgroup1, 'poissonL') && strcmp(currgroup2, 'tonic')
                groupID(numPair) = 6;
                yLabel{1,numPair} = 'PL T';
            end
            
        end
    end
    
    pairB.(currmode) = groupID == 1; % Bursting
    CCR = currData.CCR;
    
    % DELETING AND AVERAGING CENTRAL DATABINS: central bin can be
    % confounded by clustering problems
    centerP = ((length(CCR)-1)/2)+1;
    CCR(:,centerP) = mean([CCR(:,centerP-1), CCR(:,centerP+2)], 2); % Average the first and third databin
    CCR(:,centerP+1) = []; % deletes 2nd databin
    allData.(currmode) = CCR;
end

ccrB = allData.Burst1(pairB.Burst1,:);
ccrS = allData.Single(pairB.Single,:);

% Mean in the activation and baseline windows
actRatioB = mean(ccrB(:, (centerP-30):(centerP+30)),2);
baseRatioB = mean(ccrB(:, (centerP+100):end),2);
actRatioS = mean(ccrS(:, (centerP-30):(centerP+30)),2);
baseRatioS = mean(ccrS(:, (centerP+100):end),2);

% Activation index
ratioB = (actRatioB-baseRatioB)./(actRatioB+baseRatioB);
ratioS = (actRatioS-baseRatioS)./(actRatioS+baseRatioS);

% FIGURES

% CCG for Bursts vs Single spikes
uiopen([RESDIR cb fs 'groupCCG' cb '\behavior\single' cb '\avgCCG_Bursting_single_' cb '.fig'],1);
H_single = gcf;
uiopen([RESDIR cb fs 'groupCCG' cb '\behavior\burst1' cb '\avgCCG_Bursting_burst1_' cb '.fig'],1);
H_burst1 = gcf;
H_merge = figure;
A_merge = axes;
hold on
ln = findobj(H_single,'type','line');
lnn = copyobj(ln,A_merge);
set(lnn,'Color','k')
pt = findobj(H_single,'type','patch');
ptn = copyobj(pt,A_merge);
set(ptn,'FaceColor','k')
ln = findobj(H_burst1,'type','line');
lnn = copyobj(ln,A_merge);
set(lnn,'Color','m')
pt = findobj(H_burst1,'type','patch');
ptn = copyobj(pt,A_merge);
set(ptn,'FaceColor','m')
fnm = [resdir 'CCG_merge.fig'];
fnmJ = [resdir 'CCG_merge.jpeg'];
x_lim = xlim;
ylim([-0.8 5])
xticks([x_lim(1) 0 x_lim(2)]);
yticks([0 2 4])
xlabel('Lag (ms)')
ylabel('Normalized CCG')
legend({'Single spikes', 'Bursts'})
legend('boxoff')
setmyplot_tamas;
axis square;
saveas(H_merge, fnm);
saveas(H_merge, fnmJ);
close([H_burst1 H_single H_merge])

% Boxplot for synchrony index
[H1, Wp] = boxstat(ratioB,ratioS,'Burst1','Single',0.05,'paired');
setmyplot_tamas;
storeP = round(Wp,4);
close(H1);


% Bar graph for synchrony index - MEDIAN
H3 = figure;
hold on;
bar(1, median(ratioB), 'FaceColor', 'w', 'EdgeColor', 'm', 'LineWidth',2)
bar(2, median(ratioS), 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth',2)
line([1, 2], [ratioB, ratioS], 'Color', [0.6 0.6 0.6], 'LineWidth',3)
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)*0.5), y_lim(2)*0.95, ['p=' num2str(storeP)], 'Color', 'k');
xticks([]);
setmyplot_tamas;
axis square;
title('  ')
ylabel('Synchrony Index')
fnm = [resdir 'synchrony_index_burst1_single_median.fig'];
fnmJ = [resdir 'synchrony_index_burst1_single_median.jpeg'];
saveas(H3, fnm);
saveas(H3, fnmJ);
close(H3);