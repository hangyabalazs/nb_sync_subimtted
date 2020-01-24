function BF_SYNC_MASTER
%BF_SYNC_MASTER Main wrapper for the figures presented in the Laszlovszky
%   et al. manuscript. Functions are custom written or modified from Hangya
%   matlab code, or from cellbase.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

%% Prepare
% Global variables
cb = 'NB'; % select cellbase for nucleus basalis
initGlobals(cb); % initialize global variables
pathSelect(cb); % select path for saving

% Select cholinergic cells
[ChAT, pChAT, allChAT] = selectChAT(cb);
NumChAT = length(allChAT);

%% Figure 1
%--------------------------------------------------------------------------
% PANEL A
% Confocal image from Cell paper

% PANEL B
reconstruction_TPcolors(allChAT);  % Reconstruction of the recording sites(Atlas figures)

% PANEL C
% Continous raw trace examples for bursting cell
BFexample = {'n029_120220a_3.1'}; % Bursting cell example
modeSelect = 1; % 1 - Bursting, 2 - Tonic
meanDiff = 0;
rawTracePlot(BFexample, modeSelect, meanDiff);
meanDiff = 1;
rawTracePlot(BFexample, modeSelect, meanDiff);

% PANEL D
% Continous raw trace example for tonic cell
LFexample = {'n046_130101a_6.1'}; % Tonic cell example
modeSelect = 2;
meanDiff = 0;
rawTracePlot(LFexample, modeSelect, meanDiff);
meanDiff = 1;
rawTracePlot(LFexample, modeSelect, meanDiff);

% PANEL E and SUPPLEMENTARY FIGURE S2
refractoryTestFit;  % Refractory GM fit model

% PANEL F-I
acg_NBsync(allChAT,'issave',true); % RAW ACG
nbtonicphasic_NBsync;  % ACG Tonic-Phasic properties

% Downsample based on in vitro data
poolDat = load(which('ACG_matrices_POOLED.mat'));
acg_NBsync(poolDat.cellids,'issave',true, 'inVitroCut',true,...
    'pooled', true); % RAW ACG

% PANEL J-K
parameterPlot;


%% Figure 2
%--------------------------------------------------------------------------

% IN VITRO DATA FROM DANI
% AVG Spikeshapes from parameterPlot.m

%% Figure 3
%--------------------------------------------------------------------------

% PANEL A
burstRatio; % Pie chart of burst/single spikes between groups; baseline FR

% PANEL B, C, D, E
psthExtract; % Runs psth for all events
eventDataPlot; % Plots psth of the selected event

%% Figure 4
%--------------------------------------------------------------------------
[allChAT, ~] = selectPCells(cb, allChAT);   % adding alternative sessions
sessiontypes = {'behavior'};   % also possible: 'ITI', 'sleep', 'freely moving', 'quiet wakefulness'
% SET 'RUNTIME' TO 'FIRST' IN CASE OF 1ST EXECUTION, AFTER THAT YOU CAN USE THE
% VALUE 'REPEATED' WHICH LOADS THE ALREADY DETECTED PAIRS
% PANEL A, LEFT: EXAMPLE PLOTS

poolCCG = load('D:\_MATLAB_DATA\POOLED\ccgPOOLED\CCG_matrices_POOLED.mat');
allChAT = [poolCCG.PairOfCells(:,1); poolCCG.PairOfCells(:,2)];
allChAT = unique(allChAT);

for i = 1:length(sessiontypes)
%     ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'none', 'runtime', 'repeated', 'sessiontype', sessiontypes{i});% RAW CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burst1', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURST1 CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'single', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % SINGLE CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burstall', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURSTALL CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burst1Single', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURSTALL CCG
end

ccgPATH = 'D:\_MATLAB_DATA\POOLED\ccg\CCG_matrices_POOLED.mat';


% PANEL B, C
% Makes imagesc, and CCG's 
groupCCG('ccgPath', ccgPATH);
groupCCG('burstfilter', 'burst1', 'ccgPath', ccgPATH);
groupCCG('burstfilter', 'burstall', 'ccgPath', ccgPATH);
groupCCG('burstfilter', 'single', 'ccgPath', ccgPATH);
groupCCG('burstfilter', 'burst1Single', 'ccgPath', ccgPATH);

% PANEL C, RIGHT: BAR GRAPH
% Bar plots for CCG ratios
ccgCorrBars;

%% Figure 5
%--------------------------------------------------------------------------
% Performs individual STA, ERS, Phase calculations
sta_ers
sta_ers('burstfilter', 'burst1');
sta_ers('burstfilter', 'single');
sta_ers('burstfilter', 'burst1Single', 'sync', 'sync')
sta_ers('burstfilter', 'burst1Single', 'sync', 'nonsync')
sta_ers('sessType', 'burstOn');

% Figure 9
sta_ers('sessType', 'fb_exclude_new');


sta_ers('sessType', 'stimHIT');
sta_ers('sessType', 'stimFA');
sta_ers('sessType', 'CR');
sta_ers('sessType', 'Miss');

sta_ers('sessType', 'CueHit')
sta_ers('sessType', 'CueFA')
sta_ers('sessType', 'CRhit');
sta_ers('sessType', 'Missfa');

sta_ers('sessType', 'fb_exclude');



% Group plots for STA, ERS, Phase, barplots, etc.
groupStaErs
fiveBar;
zeroLagSTA;



%% Figure 6
%--------------------------------------------------------------------------

choosecb('HDB'); % Change to HDB cellbase
cb = whichcb;
[~, ~, allChAT] = selectChAT(cb);

% PANEL A
acg_NBsync(allChAT,'issave',true); % RAW ACG

% PANEL B, C, D
nbtonicphasic_NBsync;  % ACG Tonic-Phasic properties
parameterPlot;

% see also parameterPlot_HDB.m

% PANEL E
[allChAT, ~] = selectPCells(cb, allChAT); % selectting putative cholinergic cells
sessiontypes = {'behavior'}; % also possible: 'ITI', 'sleep', 'freely moving', 'quiet wakefulness'
i=1;
% SET RUNTIME TO FIRST IN CASE OF 1ST EXECUTION, AFTER THAT YOU CAN USE THE
% VALUE REPEATED WHICH LOADS THE ALREADY DETECTED PAIRS
% PANEL A EXAMPLE PLOTS
for i=1:length(sessiontypes)
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'none', 'runtime', 'repeated', 'sessiontype', sessiontypes{i});% RAW CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burst1', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURST1 CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'single', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % SINGLE CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burstall', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURSTALL CCG
    ccg_NBsync(allChAT,'issave',true,'whichcells','allpairs', 'burstfilter', 'burst1Single', 'runtime', 'repeated', 'sessiontype', sessiontypes{i}); % BURSTALL CCG
end

% Makes imagesc, and CCG's 
groupCCG;
groupCCG('burstfilter', 'burst1');
groupCCG('burstfilter', 'burstall');
groupCCG('burstfilter', 'single');
groupCCG('burstfilter', 'burst1Single');

% Bar plots for CCG ratios
ccgCorrBars;


%% Figure 7
%--------------------------------------------------------------------------

choosecb('NB'); % Change to NB cellbase
[~, ~, allChAT] = selectChAT('other');
pathSelect('other');
% MODIFY ACG, NBTONICPHASIC, PARAMETERPLOT RESDIR AND DATAPATH FOR OTHER
% CELLS. USUALLY SAVED TO A FOLDER INSIDE TAGGED CELL LOCATIONS

% PANEL A
acg_NBsync(allChAT,'issave',true); % RAW ACG

% PANEL B, C, D
nbtonicphasic_NBsync;   % ACG Tonic-Phasic properties
parameterPlot;

%% Figure 9
%--------------------------------------------------------------------------
% Figure 9
sta_ers_fig9
sta_ers_fig9('burstfilter', 'burst1');
sta_ers_fig9('burstfilter', 'single');
sta_ers_fig9('burstfilter', 'burst1Single', 'sync', 'sync')
sta_ers_fig9('burstfilter', 'burst1Single', 'sync', 'nonsync')

% Figure 9
sta_ers_fig9('sessType', 'stimHIT');
sta_ers_fig9('sessType', 'stimFA');
sta_ers_fig9('sessType', 'CR');
sta_ers_fig9('sessType', 'Miss');

% Group plots for STA, ERS, Phase, barplots, etc.
groupStaErs
zeroLagSTA;

% SUPPLEMENTARY FIGURE1
%--------------------------------------------------------------------------
[~, ~, allChAT] = selectChAT('NB');
% Run separately acg, ccg, groupCCG, psth, for chat+ and pchat using the codes above
pathSelect('NB');
nbtonicphasic_pChAT;
psthExtract; % Runs psth for all events

cb = 'HDB'; % select cellbase for nucleus basalis
choosecb(cb);
initGlobals(cb); % initialize global variables
pathSelect(cb); % select path for saving

psthExtract; % Runs psth for all events
% eventDataPlot; % Plots psth of the selected event

cb = 'PannaHDB'; % select cellbase for nucleus basalis
choosecb(cb);
initGlobals(cb); % initialize global variables
pathSelect(cb); % select path for saving

psthExtract; % Runs psth for all events
% eventDataPlot; % Plots psth of the selected event


eventDataPlot_All;

% SUPPLEMENTARY FIGURE2
%--------------------------------------------------------------------------

refractoryTestFit;

% SUPPLEMENTARY FIGURE3
%--------------------------------------------------------------------------
% Run it for Hit trials
align = 'HIT'; % values: HIT, FA
eventDataPlot('align', align); % Plots psth of the selected event

% SUPPLEMENTARY FIGURE4
%--------------------------------------------------------------------------
choosecb('HDB'); % Change to HDB cellbase
[~, ~, allChAT] = selectChAT('otherHDB');
pathSelect('otherHDB');


% MODIFY ACG, NBTONICPHASIC, PARAMETERPLOT RESDIR AND DATAPATH FOR OTHER
% CELLS. USUALLY SAVED TO A FOLDER INSIDE TAGGED CELL LOCATIONS

% PANEL A
% ACGs
acg_NBsync(allChAT,'issave',true); % RAW ACG

% PANEL B, C, D
% ACG Tonic-Phasic properties
nbtonicphasic_NBsync;
parameterPlot;

% SUPPLEMENTARY FIGURE5
%--------------------------------------------------------------------------
choosecb('NB'); % Change to HDB cellbase
[~, ~, allChAT] = selectChAT('Acb');
pathSelect('CPu');

% MODIFY ACG, NBTONICPHASIC, PARAMETERPLOT RESDIR AND DATAPATH FOR OTHER
% CELLS. USUALLY SAVED TO A FOLDER INSIDE TAGGED CELL LOCATIONS

% PANEL A
% ACGs
acg_NBsync(allChAT,'issave',true); % RAW ACG

% PANEL B, C, D
% ACG Tonic-Phasic properties
nbtonicphasic_NBsync;
parameterPlot;

psthExtract(allChAT); % Runs psth for all events
eventDataPlot; % Plots psth of the selected event
align = 'HIT'; % values: HIT, FA
eventDataPlot('align', align); % Plots psth of the selected event


% SUPPLEMENTARY FIGURE6
%--------------------------------------------------------------------------
cb = 'NB'; % select cellbase
choosecb(cb); % Change to HDB cellbase
[~, ~, allChAT] = selectChAT('ACx');
pathSelect('ACx');
allChAT = [];
allChAT{1} = 'n023_111221a_5.1';
allChAT{2} = 'n023_111215b_5.1';
allChAT{3} = 'n040_121010a_5.1';
allChAT{4} = 'n043_121204a_5.1';
allChAT{5} = 'n043_121204a_5.2';
sta_ers;

% SUPPLEMENTARY FIGURE7
%--------------------------------------------------------------------------
cb = 'NB'; % select cellbase
choosecb(cb); % Change to HDB cellbase
[~, ~, allChAT] = selectChAT('CPu');
pathSelect('CPu');
fig_raster(allChAT)
