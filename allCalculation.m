function allCalculation
%ALLCALCULATION Wrapper for scripts necessary to run
% for BF_SYNC_MASTER_REVIEW plotter functions.

%% NB CALCULATIONS
cells = 'NB';
choosecb(cells)
[nbChAT, nbpChAT, allChAT] = selectChAT(cells);
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false); % RAW ACG

pathSelect('other')
[nbChAT, nbpChAT, allChAT] = selectChAT('other');
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false); % RAW ACG



nbtonicphasic_NBsync(cells);  % ACG Tonic-Phasic properties
burstRatio(cells); % Pie chart of burst/single spikes between groups

issave = true;
actBurst = true;
psthExtract(allChAT, issave, actBurst, cells); % Runs psth for all events
doraster = true;
trialType = 'FB';
chatType = 'allC';
feedbackType = 'FA';
eventDataPlot_All(cells, actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

[allChAT, ~] = selectPCells(cb, allChAT);   % adding alternative sessions
sessiontypes = {'behavior'};   % also possible: 'ITI', 'sleep', 'freely moving', 'quiet wakefulness'
% SET 'RUNTIME' TO 'FIRST' IN CASE OF 1ST EXECUTION, AFTER THAT YOU CAN USE THE
% VALUE 'REPEATED' WHICH LOADS THE ALREADY DETECTED PAIRS
% PANEL A, LEFT: EXAMPLE PLOTS

ccgWrapperCB;
poolCCG = load('D:\_MATLAB_DATA\POOLED\ccgPOOLED\CCG_matrices_POOLED.mat');
allChAT = [poolCCG.PairOfCells(:,1); poolCCG.PairOfCells(:,2)];
allChAT = unique(allChAT);

for i = 1:length(sessiontypes)
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'none', 'runtime', 'repeated', 'sessiontype', sessiontypes{i});% RAW CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burst1', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURST1 CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'single', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % SINGLE CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burstall', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURSTALL CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burst1Single', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURSTALL CCG
end

ccgPATH = 'D:\_MATLAB_DATA\POOLED\ccgPOOLED\CCG_matrices_POOLED.mat';


% PANEL B, C
% Makes imagesc, and CCG's 
groupCCG('ccgPath', ccgPATH);
groupCCG('burstfilter', 'burst1');
groupCCG('burstfilter', 'burstall');
groupCCG('burstfilter', 'single');
groupCCG('burstfilter', 'burst1Single');

% PANEL C, RIGHT: BAR GRAPH
% Bar plots for CCG ratios
ccgCorrBars('POOLED');


%% HDB CALCULATIONS
cells = 'HDB';
choosecb(cells)
[nbChAT, nbpChAT, allChAT] = selectChAT(cells);
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false); % RAW ACG
nbtonicphasic_NBsync(cells);  % ACG Tonic-Phasic properties
burstRatio(cells); % Pie chart of burst/single spikes between groups

issave = true;
actBurst = true;
psthExtract(allChAT, issave, actBurst, cells); % Runs psth for all events
doraster = true;
trialType = 'FB';
chatType = 'allC';
feedbackType = 'FA';
eventDataPlot_All(cells, actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

%% PANNAHDB CALCULATIONS
cells = 'PannaHDB';
choosecb(cells)
[nbChAT, nbpChAT, allChAT] = selectChAT(cells);
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false); % RAW ACG
nbtonicphasic_NBsync(cells);  % ACG Tonic-Phasic properties
burstRatio(cells); % Pie chart of burst/single spikes between groups

issave = true;
actBurst = true;
psthExtract(allChAT, issave, actBurst, cells); % Runs psth for all events
doraster = true;
trialType = 'FB';
chatType = 'allC';
feedbackType = 'FA';
eventDataPlot_All(cells, actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event


%% PANNAHDB CALCULATIONS
cells = 'bothHDB';
choosecb('HDB')
[nbChAT, nbpChAT, allChAT] = selectChAT(cells);
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false); % RAW ACG
nbtonicphasic_NBsync(cells);  % ACG Tonic-Phasic properties
burstRatio(cells); % Pie chart of burst/single spikes between groups

issave = true;
actBurst = false;
psthExtract(allChAT, issave, actBurst, cells); % Runs psth for all events


%% CPU
actBurst = false;
issave = true;
% Striatal TANs
cells = 'NB';
choosecb(cells);
[~, ~, allChAT_NB] = selectChAT('CPu');

cells = 'HDB';
choosecb(cells);
[~, ~, allChAT_HDB] = selectChAT('Acb');

allChAT = [allChAT_NB allChAT_HDB];
global PATH;
PATH = '\CPu';

acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false); % RAW ACG

% SET SIGMA TO 0.01 (NORMALLY 0.002 IS USED)
psthExtract(allChAT, issave, actBurst, cells); % Runs psth for all events


% SPIKESHAPE ANALYSIS FOR RAWTRACE
[ChAT, pChAT, allChAT] = prepareScript;
NumChAT = length(allChAT);

spikeShape_wrapper(allChAT);


