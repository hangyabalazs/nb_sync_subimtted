function untaggedCCGPlotter

choosecb('NB');

UT = load('CCG_matrices_behavior_none_NB.mat');
chatType1 = cellfun(@(x) getvalue('ChAT+', x), UT.PairOfCells(:,1));
chatType2 = cellfun(@(x) getvalue('ChAT+', x), UT.PairOfCells(:,2));

pairType = chatType1 | chatType2;
untagged = find(pairType == 0);
numPair = numel(untagged);
CCR = UT.CCR(untagged,:);
PairOfCells = UT.PairOfCells(untagged,:);
for iC = 1:size(PairOfCells,1)
   tetrodeP(iC) = strcmp(PairOfCells{iC,1}(end-2), PairOfCells{iC,2}(end-2)); 
end




% DELETING AND AVERAGING CENTRAL DATABINS
plotWin = [round(((size(CCR,2)-1)/2)*-1) round(((size(CCR,2)-1)/2))]; % +/- 250ms lag
lags = [plotWin(1):1:-1, 1:1:plotWin(2)]; % +/- 1 s window, 1 ms sampling rate
centerP = round(((size(CCR,2)-1)/2))+1; % Center point of the current data
CCR(:,centerP) = mean([CCR(:,centerP-1), CCR(:,centerP+2)], 2); % Average the first and third databin
CCR(:,centerP+1) = []; % deletes 2nd databin

CCRnorm = [];
for k = 1:numPair
    sccr = smooth(CCR(k,:),'linear',15);   % smooth
    b = mean(sccr);
    normValue(k) = b;
    stdValue(k) = std(sccr);
    CCRnorm(k,:) = (sccr - normValue(k)) / stdValue(k);  % mean-normalize
end

mn_all = median(CCRnorm, 1);   % average CCG, tonic cells
mn_all = normalizem(mn_all,centerP);
se_all = se_of_median(CCRnorm,1);   % SE, tonic cells
se_all = normalizem(se_all,centerP);


H_ut = figure;
errorshade(lags,mn_all,se_all,'LineColor',[0.6 0.6 0.6],'ShadeColor',[0.6 0.6 0.6],'LineWidth',3)
% xlim([plotWin(1) plotWin(2)])
% xlabel('Lag (ms)')
% ylabel('Normalized count')

ln = findobj(H_ut,'type','line');

uiopen('D:\_MATLAB_DATA\POOLED\groupCCGPOOLED\behavior\nonePOOLED\B_PC_TC_none_POOLED_all.fig',1);
H_chat = gcf;
A_chat = gca;
hold on;
lnn = copyobj(ln,A_chat);
pt = findobj(H_ut,'type','patch');
ptn = copyobj(pt,A_chat);


plotStr= {['BURST-SB n=',num2str(numB)];
    ['BURST-PL n=',num2str(numPL)];
    ['REG n=',num2str(numT)];
    ['B-PL n=', num2str(numB_PL)];
    ['B-T n=',num2str(numB_T)];
    ['PL-T n=',num2str(numPL_T)]
    };
legend(plotStr, 'FontSize', 12, 'Location','northeast')




%% Synchrony index

actRatioTP = mean(CCRnorm(tetrodeP, (centerP-30):(centerP+30)),2);
baseRatioTP = mean(CCRnorm(tetrodeP, (centerP+100):end),2);
actRatioNTP = mean(CCRnorm(~tetrodeP, (centerP-30):(centerP+30)),2);
baseRatioNTP = mean(CCRnorm(~tetrodeP, (centerP+100):end),2);

ratioTP = (actRatioTP-baseRatioTP)./(actRatioTP+baseRatioTP);
ratioNTP = (actRatioNTP-baseRatioNTP)./(actRatioNTP+baseRatioNTP);







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
saveas(H_merge, fnm);
saveas(H_merge, fnmJ);
close([H_burst1 H_single H_merge])




keyboard;


% -------------------------------------------------------------------------
function nM = normalizem(M,cp)

sM = smooth(M,'linear',15);
bsl = mean(sM(cp+100:end));
nM = sM - bsl;
