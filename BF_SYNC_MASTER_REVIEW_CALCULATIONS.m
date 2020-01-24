function BF_SYNC_MASTER_REVIEW
%BF_SYNC_MASTER_REVIEW Main wrapper for the figures presented in the Laszlovszky
%   et al. manuscript. Functions are custom written or modified from Hangya
%   matlab code, or from cellbase. This is an updated version of the
%   original wrapper for the Nature Neuroscience Review.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu
%
% See also BF_SYNC_MASTER


% Prepare variables for execution
[ChAT, pChAT, allChAT] = prepareScript;
NumChAT = length(allChAT);

%% Figure 1
%--------------------------------------------------------------------------
% PANEL A
% Confocal image from Cell paper

% PANEL B
BFexample = {'n029_120220a_3.1'}; % Bursting cell example
modeSelect = 1; % 1 - Bursting, 2 - Tonic, 3 - Poisson
meanDiff = 0;
rawTracePlot(BFexample, modeSelect, meanDiff);
meanDiff = 1;
rawTracePlot(BFexample, modeSelect, meanDiff);

% reconstruction_TPcolors(allChAT);  % Reconstruction of the recording sites(Atlas figures)


% PANEL C
% Continous raw trace example for poisson cell
PLexample = {'HDB36_190513a_3.3'};
modeSelect = 3; % 1 - Bursting, 2 - Tonic, 3 - Poisson
meanDiff = 0;
rawTracePlot(PLexample, modeSelect, meanDiff);


% PANEL D
% Continous raw trace example for tonic cell
LFexample = {'n046_130101a_6.1'}; % Tonic cell example
modeSelect = 2;
meanDiff = 0;
rawTracePlot(LFexample, modeSelect, meanDiff);
meanDiff = 1;
rawTracePlot(LFexample, modeSelect, meanDiff);

% PANEL E and SUPPLEMENTARY FIGURE S1 H,I
refractoryTestFit('POOLED');  % Refractory GM fit model

% PANEL F-I
% Also Fig7 C, Reviewer Figure2, Figure3 F
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', true); % RAW ACG
nbtonicphasic_NBsync('POOLED');  % ACG Tonic-Phasic properties

%% Figure 2
%--------------------------------------------------------------------------
% IN VITRO DATA FROM DANI

%% Figure 3
%--------------------------------------------------------------------------
% PANEL A
% Figure from Cell paper

% PANEL B
burstRatio('POOLED'); % Pie chart of burst/single spikes between groups

% PANEL C, D, E
% Without removing non-activating bursting cells
actBurst = false;
issave = true;
% NB
cells = 'NB';
choosecb(cells);
[~, ~, allChAT] = selectChAT(cells);
psthExtract(allChAT, issave, cells); % Runs psth for all events
eventDataPlot_All(cells, actBurst); % Plots psth of the selected event

% HDB
cells = 'HDB';
choosecb(cells);
[~, ~, allChAT] = selectChAT(cells);
psthExtract(allChAT, issave, cells); % Runs psth for all events
eventDataPlot_All(cells, actBurst); % Plots psth of the selected event

% PannaHDB
cells = 'PannaHDB';
choosecb(cells);
[~, ~, allChAT] = selectChAT(cells);
psthExtract(allChAT, issave, cells); % Runs psth for all events
eventDataPlot_All(cells, actBurst); % Plots psth of the selected event

% bothHDB
cells = 'bothHDB';
[~, ~, allChAT] = selectChAT(cells);
psthExtract(allChAT, issave, cells); % Runs psth for all events
eventDataPlot_All(cells, actBurst); % Plots psth of the selected event

% POOLED
psthExtract(allChAT, issave, 'POOLED'); % Runs psth for all events
eventDataPlot_All('POOLED', actBurst); % Plots psth of the selected event

% Removing non-activating bursting cells
actBurst = true;
issave = true;
psthExtract(allChAT, issave, 'POOLED'); % Runs psth for all events
eventDataPlot_All('POOLED', actBurst); % Plots psth of the selected event


%% Figure 4
%--------------------------------------------------------------------------
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

%% Figure 5
%--------------------------------------------------------------------------
% Group plots for STA, ERS, Phase, barplots, etc.
groupStaErs
fiveBar;
zeroLagSTA;

%% Figure 6
%--------------------------------------------------------------------------

%% Figure 7
%--------------------------------------------------------------------------

%% Figure 8
%--------------------------------------------------------------------------


% SUPPLEMENTARY FIGURE4
%--------------------------------------------------------------------------
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

psthExtract(allChAT, issave, actBurst, cells); % Runs psth for all events
eventDataPlot_All(cells, actBurst); % Plots psth of the selected event

cpuPlotter;

PATH = '';


% SUPPLEMENTARY FIGURE9
%--------------------------------------------------------------------------
peakLatencySTA;



% SUPPLEMENTARY FIGURE9
%--------------------------------------------------------------------------

choosecb('NB');

% Continous raw trace examples for poisson-like cell
BFexample = {'n077_141218a_5.1'}; % Bursting cell example
modeSelect = 1; % 1 - Bursting
meanDiff = 0;
rawTracePlot(BFexample, modeSelect, meanDiff);
meanDiff = 1;
rawTracePlot(BFexample, modeSelect, meanDiff);

% REVIEWER FIGURE1
%--------------------------------------------------------------------------
% Downsample based on in vitro data
poolDat = load(which('ACG_matrices_POOLED.mat'));
acg_NBsync(poolDat.cellids,'issave',true, 'inVitroCut',true,...
    'pooled', true); % RAW ACG

inVitro_BI_results;







% REVIEWER FIGURE2
%--------------------------------------------------------------------------




% REVIEWER FIGURE5
%--------------------------------------------------------------------------
sta_ers('sessType', 'fb_exclude_new')
sta_ers('sessType', 'fb_exclude_new', 'burstfilter', 'burst1');
sta_ers('sessType', 'fb_exclude_new', 'burstfilter', 'single');
sta_ers('sessType', 'fb_exclude_new', 'burstfilter', 'burst1Single', 'sync', 'sync')
sta_ers('sessType', 'fb_exclude_new', 'burstfilter', 'burst1Single', 'sync', 'nonsync')

% Group plots for STA, ERS, Phase, barplots, etc.
groupStaErs
fiveBar;
zeroLagSTA;

% REVIEWER FIGURE6
%--------------------------------------------------------------------------

intervalSpectMean_subplot;


actBurst = false;
issave = true;
% NB
cells = 'NB';
% Prepare variables for execution
[ChAT, pChAT, allChAT] = prepareScript;
NumChAT = length(allChAT);

psthExtract(allChAT, issave, actBurst, cells); % Runs psth for all events
eventDataPlot_All(cells, actBurst); % Plots psth of the selected event

keyboard;




%% EXTRA1
% Downsample based on in vitro data
acg_NBsync(allChAT,'issave',true, 'inVitroCut',true,...
    'pooled', true); % RAW ACG

% BothHDB
% bothHDB
% set acg resdir from Pooled to bothHDB
cells = 'bothHDB';
[~, ~, allChAT] = selectChAT(cells);
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false, 'cellbase', 'bothHDB'); % RAW ACG
nbtonicphasic_NBsync('bothHDB');  % ACG Tonic-Phasic properties


cells = 'otherBothHDB';
[~, ~, allChAT] = selectChAT(cells);
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', false, 'cellbase', 'bothHDB'); % RAW ACG
nbtonicphasic_NBsync('bothHDB');  % ACG Tonic-Phasic properties

% Reliability, Latency, Jitter
relJittLat_wrapper(allChAT);
