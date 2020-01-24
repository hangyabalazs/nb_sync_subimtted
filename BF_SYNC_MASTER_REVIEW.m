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

% Calculate it with +/- 500 ms window
acg_NBsync(allChAT,'issave',true, 'inVitroCut',false,...
    'pooled', true, 'longSeg', true); % RAW ACG
nbtonicphasic_NBsync('POOLED', 'longSeg', true);  % ACG Tonic-Phasic properties


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
% Without removing non activating bursting cells
actBurst = false;
issave = true;
doraster = true;
trialType = 'FB';
chatType = 'allC';
feedbackType = 'FA';
% POOLED
eventDataPlot_All('POOLED', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

%% Figure 4
%--------------------------------------------------------------------------

ccgPATH = 'D:\_MATLAB_DATA\POOLED\ccgPOOLED\CCG_matrices_POOLED.mat';
% PANEL A, B
% Makes imagesc, and CCG's 
groupCCG('ccgPath', ccgPATH, 'cellbase', 'POOLED');

% PANEL C
% Bar plots for CCG ratios
ccgCorrBars('POOLED');

%% Figure 5
%--------------------------------------------------------------------------
% Group plots for STA, ERS, Phase, barplots, etc.
groupStaErs('sessiontype', 'burstOn');
groupStaErs('sessiontype', 'behavior', 'burstfilter', 'all');

% UPDATED: PANEL F
fiveBar;

%% Figure 6
%--------------------------------------------------------------------------
% UNCHANGED: PANEL A
% Examples

% UNCHANGED: PANEL B
eventSTA;

% UPDATED: PANEL C
zeroLagSTA;

%% Figure 7
%--------------------------------------------------------------------------

% PANEL A
% Histology from cell paper

% PANEL B
choosecb('NB')
nbtonicphasic_NBsync('NB');  % ACG Tonic-Phasic properties
choosecb('HDB')
nbtonicphasic_NBsync('bothHDB');  % ACG Tonic-Phasic properties

% PANEL C
% Already done during Figure1 creation with POOLED cb

% PANEL D-E
% Dani's figures

%% Figure 8
%--------------------------------------------------------------------------
% Done by Balazs


% SUPPLEMENTARY FIGURE1
%--------------------------------------------------------------------------

% PANEL A, D, E, F
nbtonicphasic_pChAT('POOLED');  % ACG Tonic-Phasic properties

% PANEL B1
actBurst = false;
issave = true;
doraster = true;
trialType = 'FB';
chatType = 'ChAT';
feedbackType = 'FA';
% POOLED
eventDataPlot_All('NB', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

% PANEL B2
% POOLED
eventDataPlot_All('bothHDB', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

% PANEL C
chatType = 'pChAT';
% POOLED
eventDataPlot_All('NB', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

% PANEL G1-G2
actBurst = true;
chatType = 'ChAT';
% POOLED
eventDataPlot_All('POOLED', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

% PANEL G3-G4
chatType = 'pChAT';
% POOLED
eventDataPlot_All('POOLED', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event


% PANEL H
% Already done during Figure1 plotting


% SUPPLEMENTARY FIGURE2
%--------------------------------------------------------------------------
% PANEL A-G
pathSelect('other')
nbtonicphasic_NBsync('NB');  % ACG Tonic-Phasic properties
pathSelect('NB')

% PANEL H
chatProbabilityDist;

% SUPPLEMENTARY FIGURE3
%--------------------------------------------------------------------------
% Dani's figure

% SUPPLEMENTARY FIGURE4
%--------------------------------------------------------------------------

% PANEL A-B
% Done by Balazs

% PANEL C
issave = true;
actBurst = false;
doraster = true;
trialType = 'FB';
chatType = 'allC';
feedbackType = 'HIT';
global PATH;
PATH = '\CPu';

% POOLED
eventDataPlot_All('HDB', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

feedbackType = 'FA';

eventDataPlot_All('HDB', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

cpuMerge;

% PANEL D
cpuPlotter;
PATH = '';

% Run psthExtract with higher sigma to have a smoother plot
% for the -20-80ms psth. Etc: set sigma: from 0.002 to 0.01

% PANEL E-F-G
choosecb('NB')
feedbackType = 'FA';
eventDataPlot_All('POOLED', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

feedbackType = 'HIT';
eventDataPlot_All('POOLED', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event


% SUPPLEMENTARY FIGURE5
%--------------------------------------------------------------------------
% Plots in POOLED CCG folder. Assembled in illustrator.

% SUPPLEMENTARY FIGURE6
%--------------------------------------------------------------------------
% PANEL A
peakLatencySTA;

% PANEL B-C
% Dani's figure

% SUPPLEMENTARY FIGURE7
%--------------------------------------------------------------------------
% Example plotsin STA folder, assembled in illustrator

% SUPPLEMENTARY FIGURE8
%--------------------------------------------------------------------------

pathSelect('other')
nbtonicphasic_NBsync('bothHDB');  % ACG Tonic-Phasic properties


% REVIEWER FIGURE1
%--------------------------------------------------------------------------
inVitro_BI_results;

% REVIEWER FIGURE2
%--------------------------------------------------------------------------
% Already calculated during Fig1

% REVIEWER FIGURE3
%--------------------------------------------------------------------------
% Already calculated during Fig4


% REVIEWER FIGURE4
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


% REVIEWER FIGURE5
%--------------------------------------------------------------------------
choosecb('NB')
% PANEL A
intervalSpectMean_subplot;

% PANEL B
actBurst = false;
issave = true;
% NB
cells = 'NB';
% Prepare variables for execution
[ChAT, pChAT, allChAT] = prepareScript;
NumChAT = length(allChAT);

actBurst = false;
issave = true;
doraster = true;
trialType = 'FB';
chatType = 'allC';
feedbackType = 'FA';
eventDataPlot_All('NB', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event


% REVIEWER FIGURE6
%--------------------------------------------------------------------------

% PANEL A
actBurst = false;
issave = true;
doraster = true;
trialType = 'TONE';
chatType = 'allC';
feedbackType = 'HIT';
% POOLED
eventDataPlot_All('POOLED', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

% PANEL B
actBurst = false;
issave = true;
doraster = true;
trialType = 'TRIALSTART';
chatType = 'allC';
feedbackType = 'HIT';
% POOLED
eventDataPlot_All('POOLED', actBurst, trialType, chatType, feedbackType,...
    doraster); % Plots psth of the selected event

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
    'pooled', true); % RAW ACG
nbtonicphasic_NBsync('bothHDB');  % ACG Tonic-Phasic properties


% Reliability, Latency, Jitter
relJittLat_wrapper(allChAT);
