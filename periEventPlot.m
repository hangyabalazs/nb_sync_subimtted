function periEventPlot


g.sesstype = 'behavior';
g.psthFilt = 'spsth';
% Initilize
dbstop if error;
cb = whichcb;
global RESDIR;
global cCode;
global PATH;
fs = filesep;
resdir = [RESDIR cb fs];   % results directory
psthDir = [resdir 'psth' cb fs g.sesstype cb fs 'figures' PATH fs];

% Load all psth data from psthExtract.m
load([RESDIR cb fs 'PSTH_newdata' PATH fs 'EventData_review_190829_2.mat']);

% Parameters
centerPMC = round(size(EventData.MISS.Miss.fa.NotAdaptive.none.spsth,2)/2);
centerPFH = round(size(EventData.FB.FA.NotAdaptive.none.spsth,2)/2);
time = -500:1:500; % time vector for plotting
plotWindowMC = (centerPMC+time(1)):(centerPMC+time(end));
plotWindowFH = (centerPFH+time(1)):(centerPFH+time(end));
% actWindow = centerP:(centerP+50); % activation win for burstratio
% baseWindow = (centerP+100):plotWindow(end); % base win for burstratio
LWidth = 2;
plotTime = [-500 500];
xtick_labels = {'-500', '0', '500'};
x_ticks = [-500 0 500];
tType = {'FB', 'TONE', 'MISS', 'CR'};
tType = {'MISS-CR'};


% groups = EventData.FB.FA.NotAdaptive.none.stats;
groups = EventData.MISS.Miss.fa.NotAdaptive.none.stats;

for i = 1:length(groups)
    switch groups{i}.group
        case 1
            groupID(i) = 1;
        case 2
            groupID(i) = 2;
        case 3
            groupID(i) = 3;
    end
end

phasicB = groupID==1;
poissonL = groupID==2;
tonic = groupID==3;

faTrials = EventData.TONE.FA.NotAdaptive;
hitTrials = EventData.TONE.HIT.NotAdaptive;
missTrials = EventData.MISS.Miss.fa.NotAdaptive;
crTrials = EventData.CR.CorrectRejection.hit.NotAdaptive;

currDataN_Hit = select_normalization(hitTrials);
currDataN_FA = select_normalization(faTrials);
noneData_Hit = currDataN_Hit.none.(g.psthFilt); % raw data
noneDataFA = currDataN_FA.none.(g.psthFilt); % raw data

currDataN_Miss = select_normalization(missTrials);
currDataN_CR = select_normalization(crTrials);
noneData_Miss = currDataN_Miss.none.(g.psthFilt); % raw data
noneData_CR = currDataN_CR.none.(g.psthFilt); % raw data


%     burst1Data = currDataN.burst1.(g.psthFilt); % burst1
%     burstallData = currDataN.burstall.(g.psthFilt); % burstall
%     singleData = currDataN.single.(g.psthFilt); % single
%     burst1SingleData = currDataN.burst1Single.(g.psthFilt); % burst1+single
%     burstallSingleData = currDataN.burstall.(g.psthFilt)+currData.single.(g.psthFilt); % burstall+single

% AVG and SE for Raw Data HIT
mn_phasicB_HIT = mean(noneData_Hit(phasicB, plotWindowFH),1);   % average ACG, phasic, bursting cells
se_phasicB_HIT = std(noneData_Hit(phasicB, plotWindowFH)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
mn_poissonL_HIT = mean(noneData_Hit(poissonL, plotWindowFH),1);   % average ACG, poisson-like cells
se_poissonL_HIT = std(noneData_Hit(poissonL, plotWindowFH)) / sqrt(sum(poissonL));   % SE, poisson-like cells
mn_tonic_HIT = mean(noneData_Hit(tonic, plotWindowFH),1);   % average ACG, tonic cells
se_tonic_HIT = std(noneData_Hit(tonic, plotWindowFH)) / sqrt(sum(tonic));   % SE, tonic cells

% AVG and SE for Raw Data FA
mn_phasicB_FA = mean(noneDataFA(phasicB, plotWindowFH),1);   % average ACG, phasic, bursting cells
se_phasicB_FA = std(noneDataFA(phasicB, plotWindowFH)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
mn_poissonL_FA = mean(noneDataFA(poissonL, plotWindowFH),1);   % average ACG, poisson-like cells
se_poissonL_FA = std(noneDataFA(poissonL, plotWindowFH)) / sqrt(sum(poissonL));   % SE, poisson-like cells
mn_tonic_FA = mean(noneDataFA(tonic, plotWindowFH),1);   % average ACG, tonic cells
se_tonic_FA = std(noneDataFA(tonic, plotWindowFH)) / sqrt(sum(tonic));   % SE, tonic cells

% AVG and SE for Raw Data MISS
mn_phasicB_MISS = mean(noneData_Miss(phasicB, plotWindowMC),1);   % average ACG, phasic, bursting cells
se_phasicB_MISS = std(noneData_Miss(phasicB, plotWindowMC)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
mn_poissonL_MISS = mean(noneData_Miss(poissonL, plotWindowMC),1);   % average ACG, poisson-like cells
se_poissonL_MISS = std(noneData_Miss(poissonL, plotWindowMC)) / sqrt(sum(poissonL));   % SE, poisson-like cells
mn_tonic_MISS = mean(noneData_Miss(tonic, plotWindowMC),1);   % average ACG, tonic cells
se_tonic_MISS = std(noneData_Miss(tonic, plotWindowMC)) / sqrt(sum(tonic));   % SE, tonic cells

% AVG and SE for Raw Data CR
mn_phasicB_CR = mean(noneData_CR(phasicB, plotWindowMC),1);   % average ACG, phasic, bursting cells
se_phasicB_CR = std(noneData_CR(phasicB, plotWindowMC)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
mn_poissonL_CR = mean(noneData_CR(poissonL, plotWindowMC),1);   % average ACG, poisson-like cells
se_poissonL_CR = std(noneData_CR(poissonL, plotWindowMC)) / sqrt(sum(poissonL));   % SE, poisson-like cells
mn_tonic_CR = mean(noneData_CR(tonic, plotWindowMC),1);   % average ACG, tonic cells
se_tonic_CR = std(noneData_CR(tonic, plotWindowMC)) / sqrt(sum(tonic));   % SE, tonic cells

%% FIGURES
% FIGURE 1 BFCN-burstSB
% Bursting
H1 = figure;
hold on
set(gcf, 'Renderer', 'painters');
errorshade(time, mn_phasicB_HIT-min(mn_phasicB_HIT), se_phasicB_HIT, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
errorshade(time, mn_phasicB_FA-min(mn_phasicB_FA), se_phasicB_FA, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_phasicB_CR-min(mn_phasicB_CR), se_phasicB_CR, 'LineColor', [0 0 0.75], 'ShadeColor', [0 0 0.75], 'LineWidth', LWidth);
errorshade(time, mn_phasicB_MISS-min(mn_phasicB_MISS), se_phasicB_MISS, 'LineColor', [0 0 0], 'ShadeColor', [0 0 0], 'LineWidth', LWidth);
title(['Burst-SB'])
legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from cue (ms)']);
ylabel('Firing rate (Hz)');
axis square;
ylim([-2 20]);
fName = [psthDir 'reviewPlots' fs '_BurstSB.fig'];   % save PSTH figure
fNameJ = [psthDir 'reviewPlots' fs '_BurstSB.jpeg'];   % save PSTH figure
saveas(H1,fName);
saveas(H1,fNameJ);
% close(H1);


% FIGURE 2 BFCN-burstPL
% Bursting
H1 = figure;
hold on
set(gcf, 'Renderer', 'painters');
errorshade(time, mn_poissonL_HIT-min(mn_poissonL_HIT), se_poissonL_HIT, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL_FA-min(mn_poissonL_FA), se_poissonL_FA, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL_CR-min(mn_poissonL_CR), se_poissonL_CR, 'LineColor', [0 0 0.75], 'ShadeColor', [0 0 0.75], 'LineWidth', LWidth);
errorshade(time, mn_poissonL_MISS-min(mn_poissonL_MISS), se_poissonL_MISS, 'LineColor', [0 0 0], 'ShadeColor', [0 0 0], 'LineWidth', LWidth);
title(['Burst-PL'])
legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from cue (ms)']);
ylabel('Firing rate (Hz)');
axis square;
ylim([-2 20]);
fName = [psthDir 'reviewPlots' fs '_BurstPL.fig'];   % save PSTH figure
fNameJ = [psthDir 'reviewPlots' fs '_BurstPL.jpeg'];   % save PSTH figure
saveas(H1,fName);
saveas(H1,fNameJ);
% close(H1);

% FIGURE 3 BFCN-REG
% Bursting
H1 = figure;
hold on
set(gcf, 'Renderer', 'painters');
errorshade(time, mn_tonic_HIT-min(mn_tonic_HIT), se_tonic_HIT, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic_FA-min(mn_tonic_FA), se_tonic_FA, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic_CR-min(mn_tonic_CR), se_tonic_CR, 'LineColor', [0 0 0.75], 'ShadeColor', [0 0 0.75], 'LineWidth', LWidth);
errorshade(time, mn_tonic_MISS-min(mn_tonic_MISS), se_tonic_MISS, 'LineColor', [0 0 0], 'ShadeColor', [0 0 0], 'LineWidth', LWidth);
title(['REG'])
legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from cue (ms)']);
ylabel('Firing rate (Hz)');
axis square;
ylim([-2 20]);
fName = [psthDir 'reviewPlots' fs '_Regular.fig'];   % save PSTH figure
fNameJ = [psthDir 'reviewPlots' fs '_Regular.jpeg'];   % save PSTH figure
saveas(H1,fName);
saveas(H1,fNameJ);
% close(H1);
close all;



%%
% GROUPED BY EVENTTYPES

% FIGURE 2 BFCN-burstPL
% Bursting
H4 = figure;
hold on
set(gcf, 'Renderer', 'painters');

subplot(2,2,1)
errorshade(time, mn_phasicB_HIT-min(mn_phasicB_HIT), se_phasicB_HIT, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL_HIT-min(mn_poissonL_HIT), se_poissonL_HIT, 'LineColor', cCode(2,:), 'ShadeColor', cCode(2,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic_HIT-min(mn_tonic_HIT), se_tonic_HIT, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
title(['HIT'])
legend({'Burst-SB', 'Burst-PL', 'Reg'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from cue (ms)']);
ylabel('Firing rate (Hz)');
axis square;
ylim([-2 20]);

subplot(2,2,2)
errorshade(time, mn_phasicB_FA-min(mn_phasicB_FA), se_phasicB_FA, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL_FA-min(mn_poissonL_FA), se_poissonL_FA, 'LineColor', cCode(2,:), 'ShadeColor', cCode(2,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic_FA-min(mn_tonic_FA), se_tonic_FA, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
title(['FA'])
legend({'Burst-SB', 'Burst-PL', 'Reg'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from cue (ms)']);
ylabel('Firing rate (Hz)');
axis square;
ylim([-2 20]);

subplot(2,2,3)
errorshade(time, mn_phasicB_CR-min(mn_phasicB_CR), se_phasicB_CR, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL_CR-min(mn_poissonL_CR), se_poissonL_CR, 'LineColor', cCode(2,:), 'ShadeColor', cCode(2,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic_CR-min(mn_tonic_CR), se_tonic_CR, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
title(['CR'])
legend({'Burst-SB', 'Burst-PL', 'Reg'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from cue (ms)']);
ylabel('Firing rate (Hz)');
axis square;
ylim([-2 20]);

subplot(2,2,4)
errorshade(time, mn_phasicB_MISS-min(mn_phasicB_MISS), se_phasicB_MISS, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL_MISS-min(mn_poissonL_MISS), se_poissonL_MISS, 'LineColor', cCode(2,:), 'ShadeColor', cCode(2,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic_MISS-min(mn_tonic_MISS), se_tonic_MISS, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
title(['MISS'])
legend({'Burst-SB', 'Burst-PL', 'Reg'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from cue (ms)']);
ylabel('Firing rate (Hz)');
axis square;
ylim([-2 20]);

maximize_figure;






% -------------------------------------------------------------------------
function currDataNorm = select_normalization(currData)

fld = fieldnames(currData);
NumFields = length(fld);
currDataNorm = currData;
for iF = 1:NumFields
    nfcn = @normfun0;  % no normalization
    if ismember(iF,[2 3])
        nfcn = @normfun3;   % divide by baseline (relative increase)
    end
    fld2 = fieldnames(currData.(fld{iF}));
    for iF2 = 1:3
        lpsth = currData.(fld{iF}).(fld2{iF2});
        if ~ismember(iF,[2 3])
            nlpsth = feval(nfcn,lpsth);  % call normalization function
        else
            nlpsth = feval(nfcn,lpsth);  % normalization: divide by baseline (relative increase) for both burst1 and single for proper comparison
        end
        currDataNorm.(fld{iF}).(fld2{iF2}) = nlpsth;
    end
end

% -------------------------------------------------------------------------
function N = normfun0(M)

N = M;  % no normalization - note that Z-scoring will also abolish all vs burs1+single differences

% -------------------------------------------------------------------------
function [N, mn, sd] = normfun1(M)

s = size(M);
mn = repmat(mean(M,2),1,s(2));
sd = repmat(std(M,[],2),1,s(2));
N = (M - mn) ./ sd;   % Z-score
N = nan2zero(N);

% -------------------------------------------------------------------------
function [N, mn, sd] = normfun2(M,mn,sd)

N = (M - mn) ./ sd;   % Z-score with input mean and SD
N = nan2zero(N);

% -------------------------------------------------------------------------
function [N, mn] = normfun3(M)

s = size(M);
baseWindow = 600:750;
mn = repmat(mean(M(baseWindow),2),1,s(2));
N = M ./ (1 + mn);   % devide by baseline mean
N = nan2zero(N);