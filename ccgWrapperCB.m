function ccgWrapperCB
% CCGWRAPPERCB
% Loads CCG data from all the cellbases and saves POOLED data

global RESDIR;

% Load CCG matrices
NB = load('D:\_MATLAB_DATA\NB\ccgNB\behavior\noneNB\CCG_matrices_behavior_none_NB.mat');
HDB = load('D:\_MATLAB_DATA\HDB\ccgHDB\behavior\noneHDB\CCG_matrices_behavior_none_HDB.mat');
P_NTP = load('D:\_MATLAB_DATA\PannaHDB\ccgPannaHDB\nontetrodepairsCCG_matrices.mat');
ntpCells = P_NTP.PairOfCells;
[P_NTP.group1, P_NTP.group2] = defineGroups(ntpCells);


PairOfCells = [NB.PairOfCells; HDB.PairOfCells; P_NTP.PairOfCells];
CCR = [NB.CCR; HDB.CCR; P_NTP.CCR];
LCCR = [NB.LCCR; HDB.LCCR; P_NTP.LCCR];
UCCR = [NB.UCCR; HDB.UCCR; P_NTP.UCCR];
MeanH0 = [NB.MeanH0; HDB.MeanH0; P_NTP.MeanH0];
SDH0 = [NB.SDH0; HDB.SDH0; P_NTP.SDH0];
SegmentLength = [NB.SegmentLength; HDB.SegmentLength; P_NTP.SegmentLength];
group1 = [NB.group1'; HDB.group1'; P_NTP.group1];
group2 = [NB.group2'; HDB.group2'; P_NTP.group2];

resdir = [RESDIR 'POOLED\ccgPOOLED\'];
fnmm = 'CCG_matrices_POOLED.mat';
save(fullfile(resdir,fnmm),'PairOfCells','CCR','LCCR','UCCR','MeanH0',...
    'SDH0','SegmentLength', 'group1', 'group2');


%--------------------------------------------------------------------------
% LOCAL FUNCTIONS

function [group1, group2] = defineGroups(PairOfCells)

load('D:\_MATLAB_DATA\POOLED\acgPOOLED\ACG_matrices_POOLED.mat');

for iP = 1:size(PairOfCells,1)
    inx1 = find(strcmp(cellids, PairOfCells(iP,1)));
    inx2 = find(strcmp(cellids, PairOfCells(iP,2)));
    
    % Group1
    if Refractory(inx1)>= 40
        group1{iP,1} = {'tonic'};
    elseif Refractory(inx1) < 40 & BurstIndex(inx1) <= 0.2
        group1{iP,1} = {'poissonL'};
    elseif Refractory(inx1) < 40 & BurstIndex(inx1) > 0.2
        group1{iP,1} = {'phasicB'};
    else
        disp('Bug!!')
        keyboard;
    end
    
    % Group2
    if Refractory(inx2)>= 40
        group2{iP,1} = {'tonic'};
    elseif Refractory(inx2) < 40 & BurstIndex(inx2) <= 0.2
        group2{iP,1} = {'poissonL'};
    elseif Refractory(inx2) < 40 & BurstIndex(inx2) > 0.2
        group2{iP,1} = {'phasicB'};
    else
        disp('Bug!!')
        keyboard;
    end
    
end
