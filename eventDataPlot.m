function eventDataPlot(varargin)
%EVENTDATAPLOT Reads EventData.mat file from psthExtract containing all psth data.
%   EVENTDATAPLOT(VARARGIN) Input arguments can be:
%   Trialtype: Feedback, Tone, Trialstart, ITI; 
%   Align: Hit, FA;
%   sessType: behavior, burstOn, burstOff; 
%   psthFilt: psth, spsth, spsth_se
%   adapFilt: NotAdaptive, Adaptive, DoublyAdaptive
%   chatTest: selection for all or putative cholinergic cells - true, false
%   Plots psth for all spikes, burst1+single spikes, bar graph for event
%   selectivity, psth for bursts vs single for bursting cells, individual
%   psth and raster for all cells.
%
%   See also PSTHEXTRACT.

%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

prs = inputParser;
addParameter(prs,'sesstype','behavior',@(s)ischar(s)&ismember(s,{'behavior',...
    'burstOn', 'burstOff'}))   % sessiontype filter
addParameter(prs,'trialtype','Feedback',@(s)ischar(s)&ismember(s,{'Feedback',...
    'Tone', 'Trialstart', 'ITI', 'Miss', 'CorrectRejection'}))   % trialtype filter
addParameter(prs,'psthFilt','spsth',@(s)ischar(s)&ismember(s,{'spsth',...
    'psth', 'spsth_se'}))   % psth filter
addParameter(prs,'adapFilt','NotAdaptive',@(s)ischar(s)&ismember(s,{'NotAdaptive',...
    'Adaptive', 'DoublyAdaptive'}))   % Adaptivity filter
addParameter(prs,'align','FA',@(s)ischar(s)&ismember(s,{'FA',...
    'HIT'}))   % Rewarded of punished trials
addParameter(prs,'chatTest',false,@(s)islogical(s)|ismember(s,[0 1 2]))  % Control for ChAT/pChAT test
parse(prs,varargin{:})
g = prs.Results;

% Initilize
dbstop if error;
cb = whichcb;
global RESDIR;
global cCode;
global PATH;
fs = filesep;
resdir = [RESDIR cb fs];   % results directory
psthDir = [resdir 'psth' cb fs g.sesstype cb fs 'figures' PATH fs];

% load([RESDIR cb fs 'PSTH_newdata' PATH fs 'EventData.mat']);

% if g.chatTest
%     load([RESDIR cb fs 'PSTH_newdata' PATH fs 'EventData.mat']);
% else
%     % Load all psth data from psthExtract.m
load([RESDIR cb fs 'PSTH_newdata' cb PATH fs 'EventData_review_08-Nov-2019.mat']);
% end
% Parameters
centerP = round(size(EventData.FB.FA.NotAdaptive.none.spsth,2)/2);
time = -250:1:250; % time vector for plotting
plotWindow = (centerP+time(1)):(centerP+time(end));
actWindow = centerP:(centerP+50); % activation win for burstratio
baseWindow = (centerP+100):plotWindow(end); % base win for burstratio
LWidth = 3;
plotTime = [-20 80];
xtick_labels = {'0', '40', '80'};
x_ticks = [0 40 80];

% Corresponding struct fields to current trialtype
switch g.trialtype
    case 'Feedback'
        tT = 'FB';
    case 'Tone'
        tT = 'TONE';
    case 'Trialstart'
        tT = 'TRIALSTART';
    case 'ITI'
        tT = 'ITI';
end

currData = EventData.(tT).(g.align).(g.adapFilt); % stores data selected by filters
[currData, ~] = clearEmpty(currData); % Removes empty cells
if isfield(currData.none, 'cellid') % store cellids from data file
    allChAT =currData.none.cellid; % can differ from the result of selectChat
else
    allChAT = EventData.CELLID;
end
groups = currData.none.stats;

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

switch g.chatTest
    case 0 % AllChAT
        phasicB = groupID==1;
        poissonL = groupID==2;
        tonic = groupID==3;
        PCtitle = '';
    case 1 % Tagged cholinergic cells
        ChAT = getvalue('ChAT+',allChAT)';
        phasicB = groupID==1 & ChAT ==1;
        poissonL = groupID==2 & ChAT ==1;
        tonic = groupID==3 & ChAT ==1;
        PCtitle = '_tagged';
    case 2 % Putative cholinergic cells
        pChAT = getvalue('pChAT+',allChAT)';
        phasicB = groupID==1 & pChAT ==1;
        poissonL = groupID==2 & pChAT ==1;
        tonic = groupID==3 & pChAT ==1;
        PCtitle = '_putative';
end

currDataN = select_normalization(currData);
noneData = currDataN.none.(g.psthFilt); % raw data
burst1Data = currDataN.burst1.(g.psthFilt); % burst1
burstallData = currDataN.burstall.(g.psthFilt); % burstall
singleData = currDataN.single.(g.psthFilt); % single
burst1SingleData = currDataN.burst1Single.(g.psthFilt); % burst1+single
burstallSingleData = currDataN.burstall.(g.psthFilt)+currData.single.(g.psthFilt); % burstall+single

% AVG and SE for Raw Data
mn_phasicB = mean(noneData(phasicB, plotWindow),1);   % average ACG, phasic, bursting cells
se_phasicB = std(noneData(phasicB, plotWindow)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
mn_poissonL = mean(noneData(poissonL, plotWindow),1);   % average ACG, poisson-like cells
se_poissonL = std(noneData(poissonL, plotWindow)) / sqrt(sum(poissonL));   % SE, poisson-like cells
mn_tonic = mean(noneData(tonic, plotWindow),1);   % average ACG, tonic cells
se_tonic = std(noneData(tonic, plotWindow)) / sqrt(sum(tonic));   % SE, tonic cells

% AVG and SE for Burstall+Single Data
BSmn_phasicB = mean(burstallSingleData(phasicB, plotWindow),1);   % average ACG, phasic, bursting cells
BSse_phasicB = std(burstallSingleData(phasicB, plotWindow)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
BSmn_poissonL = mean(burstallSingleData(poissonL, plotWindow),1);   % average ACG, poisson-like cells
BSse_poissonL = std(burstallSingleData(poissonL, plotWindow)) / sqrt(sum(poissonL));   % SE, poisson-like cells
BSmn_tonic = mean(burstallSingleData(tonic, plotWindow),1);   % average ACG, tonic cells
BSse_tonic = std(burstallSingleData(tonic, plotWindow)) / sqrt(sum(tonic));   % SE, tonic cells

% AVG and SE for Single spikes
Smn_phasicB = mean(singleData(phasicB, plotWindow),1);   % average ACG, phasic, bursting cells
Sse_phasicB = std(singleData(phasicB, plotWindow)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
Smn_poissonL = mean(singleData(poissonL, plotWindow),1);   % average ACG, poisson-like cells
Sse_poissonL = std(singleData(poissonL, plotWindow)) / sqrt(sum(poissonL));   % SE, poisson-like cells
Smn_tonic = mean(singleData(tonic, plotWindow),1);   % average ACG, tonic cells
Sse_tonic = std(singleData(tonic, plotWindow)) / sqrt(sum(tonic));   % SE, tonic cells

% AVG and SE for Burst1 Data
Bmn_phasicB = mean(burst1Data(phasicB, plotWindow),1);   % average ACG, phasic, bursting cells
Bse_phasicB = std(burst1Data(phasicB, plotWindow)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
Bmn_poissonL = mean(burst1Data(poissonL, plotWindow),1);   % average ACG, poisson-like cells
Bse_poissonL = std(burst1Data(poissonL, plotWindow)) / sqrt(sum(poissonL));   % SE, poisson-like cells
Bmn_tonic = mean(burst1Data(tonic, plotWindow),1);   % average ACG, tonic cells
Bse_tonic = std(burst1Data(tonic, plotWindow)) / sqrt(sum(tonic));   % SE, tonic cells

% AVG and SE for Burst1+Single Data
B1Smn_phasicB = mean(burst1SingleData(phasicB, plotWindow),1);   % average ACG, phasic, bursting cells
B1Sse_phasicB = std(burst1SingleData(phasicB, plotWindow)) / sqrt(sum(phasicB));   % SE, phasic, bursting cells
B1Smn_poissonL = mean(burst1SingleData(poissonL, plotWindow),1);   % average ACG, poisson-like cells
B1Sse_poissonL = std(burst1SingleData(poissonL, plotWindow)) / sqrt(sum(poissonL));   % SE, poisson-like cells
B1Smn_tonic = mean(burst1SingleData(tonic, plotWindow),1);   % average ACG, tonic cells
B1Sse_tonic = std(burst1SingleData(tonic, plotWindow)) / sqrt(sum(tonic));   % SE, tonic cells

%--------------------------------------------------------------------------
% FIGURES
% Figure for Raw
H1 = figure;
set(gcf, 'Renderer', 'painters');
hold on;
errorshade(time, mn_phasicB-min(mn_phasicB), se_phasicB, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL-min(mn_poissonL), se_poissonL, 'LineColor', cCode(2,:), 'ShadeColor', cCode(2,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic-min(mn_tonic), se_tonic, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
title(['All action potentials ' PCtitle(2:end)])
legend({'Phasic Bursting','Poisson-like','Tonically active'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
setmyplot_tamas;
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel('Time from punishment (ms)');
ylabel('Firing rate (Hz)');
axis square;
fName = [psthDir 'raw_psth_' g.trialtype '_' g.align PCtitle '.fig'];   % save PSTH figure
fNameJ = [psthDir 'raw_psth_' g.trialtype '_' g.align PCtitle '.jpeg'];   % save STA jpeg
saveas(H1,fName);
saveas(H1,fNameJ);
close(H1);

% Figure for burst1+single
H2 = figure;
set(gcf, 'Renderer', 'painters');
hold on;
errorshade(time, B1Smn_phasicB-min(B1Smn_phasicB), B1Sse_phasicB, 'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:), 'LineWidth', LWidth);
errorshade(time, B1Smn_poissonL-min(B1Smn_poissonL), B1Sse_poissonL, 'LineColor', cCode(2,:), 'ShadeColor', cCode(2,:), 'LineWidth', LWidth);
errorshade(time, B1Smn_tonic-min(B1Smn_tonic), B1Sse_tonic, 'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth', LWidth);
title(['Bursts counted as single events ' PCtitle(2:end)])
setmyplot_tamas;
legend({'Phasic Bursting','Poisson-like','Tonically active'}, 'FontSize', 12, 'Location','northeast');
xlim(plotTime);
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel('Time from punishment (ms)');
ylabel('Firing rate (Hz)');
axis square;
fName = [psthDir 'burst1Single_psth_' g.trialtype '_' g.align PCtitle '.fig'];   % save PSTH figure
fNameJ = [psthDir 'burst1Single_psth_' g.trialtype '_' g.align PCtitle '.jpeg'];   % save STA jpeg
saveas(H2,fName);
saveas(H2,fNameJ);
close(H2);

% Burst selectivity
H3 = figure;
set(gcf, 'Renderer', 'painters');
hold on;
errorshade(time, Bmn_phasicB-min(Bmn_phasicB), Bse_phasicB, 'LineColor', 'm', 'ShadeColor', 'm', 'LineWidth', LWidth);
errorshade(time, Smn_phasicB-min(Smn_phasicB), Sse_phasicB, 'LineStyle', ':', 'LineColor', [0 0 0], 'ShadeColor', [0 0 0], 'LineWidth', LWidth);
title(['Burst selectivity ' PCtitle(2:end)])
setmyplot_tamas;
legend({'Bursts', 'Single spikes'}, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
xticks(x_ticks);
xticklabels(xtick_labels);
xlim(plotTime);
xlabel('Time from punishment (ms)');
ylabel('Firing rate (Hz)');
axis square;
fName = [psthDir 'burstSelectivity_' g.trialtype '_' g.align PCtitle '.fig'];   % save PSTH figure
fNameJ = [psthDir 'burstSelectivity_' g.trialtype '_' g.align PCtitle '.jpeg'];   % save STA jpeg
saveas(H3,fName);
saveas(H3,fNameJ);
close(H3);

%BurstIndex (Act. window mean - base win mean / act+base win means)
burstInd = (mean(burst1Data(phasicB,actWindow),2)-mean(burst1Data(phasicB,baseWindow),2))./(mean(burst1Data(phasicB,actWindow),2)+mean(burst1Data(phasicB,baseWindow),2));
singleInd = (mean(singleData(phasicB,actWindow),2)-mean(singleData(phasicB,baseWindow),2))./(mean(singleData(phasicB,actWindow),2)+mean(singleData(phasicB,baseWindow),2));
burstInd(isnan(burstInd)) = 0;
singleInd(isnan(singleInd)) = 0;

% Boxplot for burstIndexes
[H4, Wp] = boxstat(burstInd,singleInd,'Burst1','Single AP',0.05,'paired');
title(['Burst1 vs Single ' PCtitle(2:end)])
setmyplot_tamas;
storeP = round(Wp,4);
close(H4);

% Barplot for burstIndexes median
H5 = figure;
hold on;
line([1, 2], [burstInd, singleInd], 'Color', [0.6 0.6 0.6], 'LineWidth',3)
bar(1, median(burstInd), 'FaceColor','none','EdgeColor','m','LineWidth',3)
bar(2, median(singleInd), 'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',3)
xticks([0.9 2.1]);
xticklabels({'Bursts', 'Single spikes'});
ylabel('Selectivity Index')
axis square;
setmyplot_tamas;
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)*0.7), y_lim(2)*0.95, ['p: ' num2str(storeP)]);
setmyplot_tamas;
fName = [psthDir 'burstIndex_' g.trialtype '_' g.align PCtitle '_median.fig'];   % save PSTH figure
fNameJ = [psthDir 'burstIndex_' g.trialtype '_' g.align PCtitle '_median.jpeg'];   % save STA jpeg
saveas(H5,fName);
saveas(H5,fNameJ);
close(H5);

% Barplot for burstIndexes mean
H5 = figure;
hold on;
line([1, 2], [burstInd, singleInd], 'Color', [0.6 0.6 0.6], 'LineWidth',3)
bar(1, mean(burstInd), 'FaceColor','none','EdgeColor','m','LineWidth',3)
bar(2, mean(singleInd), 'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',3)
xticks([0.9 2.1]);
xticklabels({'Bursts', 'Single spikes'});
ylabel('Selectivity Index')
axis square;
setmyplot_tamas;
x_lim = xlim;
y_lim = ylim;
text((x_lim(2)*0.7), y_lim(2)*0.95, ['p: ' num2str(storeP)]);
setmyplot_tamas;
fName = [psthDir 'burstIndex_' g.trialtype '_' g.align PCtitle '_mean.fig'];   % save PSTH figure
fNameJ = [psthDir 'burstIndex_' g.trialtype '_' g.align PCtitle '_mean.jpeg'];   % save STA jpeg
saveas(H5,fName);
saveas(H5,fNameJ);
close(H5);

% Individual raster and psth plots for all cells
for aI = 1:length(allChAT)
    %Store raster data
    rasterNone = currData.none.raster{aI};
    rasterBurst1 = currData.burst1.raster{aI};
    cellidt = regexprep(allChAT{aI},'_',' ');
    
    switch groupID(aI)
        case 1
            group = 'PhasicB';
        case 2
            group = 'PoissonL';
        case 3
            group = 'Tonic';
    end
    
    fnmS1 = [resdir 'raster' cb '\individuals\' group '_' cellidt '_' g.align '_' g.trialtype PCtitle '_raster.fig'];
    fnmS2 = [resdir 'raster' cb '\individuals\' group '_' cellidt '_' g.align '_' g.trialtype PCtitle '_psth.fig'];
    fnmSE1 = [resdir 'raster' cb '\individuals\' group '_' cellidt '_' g.align '_' g.trialtype PCtitle '_raster.jpeg'];
    fnmSE2 = [resdir 'raster' cb '\individuals\' group '_' cellidt '_' g.align '_' g.trialtype PCtitle '_psth.jpeg'];
    
    % Creates plot matrix (1 - Single spikes, 2 - Burst1)
    plotImage = rasterNone(:,plotWindow)+rasterBurst1(:,plotWindow);
    middle = round(size(plotImage,2)/2);
    
    % Individual raster plot
    H6 = figure;
    hold on;
    imagesc(plotImage);
    if isempty(find(plotImage==2))
        map = [1 1 1; 0 0 0];
    else
        map = [1 1 1; 0 0 0; 1 0 0];
    end
    colormap(map);
    setmyplot_tamas;
    line([middle; middle], ylim, 'LineStyle', ':', 'Color', [0.8 0 0], 'LineWidth', LWidth);
    yticks(ylim)
    yticklabels({num2str(1), num2str(size(rasterNone, 1))})
    xticks(x_ticks);
    xticklabels(xtick_labels);
    ylabel('#Trials');
    xlabel('Time from punishment (ms)')
    title([cellidt ' ' group ' ' PCtitle(2:end)])
    ylabel('#Trials');
    axis square;
    xlim([middle middle]+plotTime);
    saveas(H6, fnmS1);
    saveas(H6, fnmSE1);
    close(H6);
    
    % Individual PSTH plot
    H7 = figure;
    hold on;
    plot(time, noneData(aI,plotWindow), 'Color', 'k', 'LineWidth', LWidth);
    xlabel('Time from punishment (ms)');
    ylabel('Firing rate (Hz)');
    setmyplot_tamas;
    xticks(x_ticks);
    xticklabels(xtick_labels);
    title([cellidt ' Group: ' group ' ' PCtitle(2:end)])
    axis square;
    line([0; 0], ylim, 'LineStyle', ':', 'Color', [0.8 0 0], 'LineWidth', LWidth);
    xlim(plotTime);
    saveas(H7,fnmS2);
    saveas(H7,fnmSE2);
    close(H7); 
    
    % Merged-Plot
    H8 = figure;
    hold on;
    imagesc(plotImage);
    if isempty(find(plotImage==2))
        map = [1 1 1; 0 0 0];
    else
        map = [1 1 1; 0 0 0; 1 0 0];
    end
    colormap(map);
    setmyplot_tamas;
    plot(time, noneData(aI,plotWindow), 'Color', 'k', 'LineWidth', LWidth);
    line([middle; middle], ylim, 'Color', [0 0 0.8], 'LineWidth', LWidth);
    yticks(ylim)
    yticklabels({num2str(1), num2str(size(rasterNone, 1))})
    xticks(x_ticks);
    xticklabels(xtick_labels);
    ylabel('#Trials');
    xlabel('Time from punishment (ms)')
    title([cellidt ' ' group ' ' PCtitle(2:end)])
    ylabel('#Trials');
    axis square;
end


% -------------------------------------------------------------------------
function [currData, indexes] = clearEmpty(currData)

indexes = [];
filters = fieldnames(currData);

for i=1:length(filters)
    indexes = find(cellfun(@isempty,currData.(filters{i}).stats));
    currData.(filters{i}).stats(indexes)=[];
    currData.(filters{i}).psth(indexes,:) = [];
    currData.(filters{i}).spsth(indexes,:) = [];
    currData.(filters{i}).spsth_se(indexes,:) = [];
    currData.(filters{i}).raster(indexes) = [];
    if isfield(currData.(filters{i}), 'cellid')
        currData.(filters{i}).cellid(indexes) = [];
    end
end

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