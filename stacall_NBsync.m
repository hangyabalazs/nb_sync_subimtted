function [st, sta, stats, H] = stacall_NBsync(vdisc,eeg,cellid, varargin)
%STACALL   Calls ASTANORM2 and STAFIG.
%   [ST STA STATS H] = STACALL(VDISC,EEG,SR,wn) calculates and plots Spike Triggered Average.
%   Input arguments:
%       VDISC: point process
%       EEG: continuous signal
%       CELLID: current cellid
%       SR (optional): sampling rate (default 1000 Hz)
%       WN (optional): window size in data points (default: 1 sec)
%       BURSTFILTER (optional): controls filterint properties (default none)
%       SYNC (optional): in case of cellpairs, filter for sync events (default: true)
%       CELLID2 (optional): in case of a pair the 2nd cellid (default: [])
%   Output argument:
%       ST: aligned data matrix
%       STA: spike triggered average
%       STATS: structure matrix of STA statistics with the followong fields:
%           CV: STA coefficient of variation (mean / SD)
%           STA_INDEX1: maximum STA
%           STA_INDEX2: maximum STA minus mean STA
%           STA_AMPLITUDE: STA amplitude
%           SPIKE_NUMBER: number of spikes used for the calculation
%           LOWER_LIMIT: lower 95% confidence limit
%           UPPER_LIMIT: upper 95% confidence interval
%       H: figure handle of the STA plot
%
%   See also ASTANORM2 and STAFIG.
% Default arguments
prs = inputParser;
addRequired(prs,'vdisc',@(s)isnumeric(s))
addRequired(prs,'eeg',@(s)isnumeric(s))
addRequired(prs,'cellid',@(s)iscellid(s)|issessionid(s))
addParameter(prs,'sr', 1000, @isnumeric)
addParameter(prs,'wn', 1,@isnumeric)
addParameter(prs,'burstfilter','none',@(s)ischar(s)&ismember(s,{'none',...
    'burst1', 'burstall', 'single', 'burst1Single'}))   % filter on burst first spikes or single spikes
addParameter(prs,'sync',true,@islogical)
addParameter(prs,'cellid2', [], @iscellid)
addParameter(prs,'sessType','behavior',@(s)ischar(s)&ismember(s,{'behavior',...
    'burstOn', 'stimFA', 'stimHIT', 'CR', 'Miss', 'CueHit', 'CueFA',...
    'CRhit', 'Missfa', 'fb_exclude', 'fb_exclude_new'}))   % session type filter
parse(prs,vdisc, eeg, cellid, varargin{:});
g = prs.Results;

cellbase = whichcb;
global RESDIR;
global PATH;
inverted = false;
isnorm = 0;
[sta, sta_cv, sta_index1, sta_index2, sta_amp, lenv, st] = astanorm2(vdisc,eeg,g.sr,g.wn,isnorm);
absWin = 1:length(sta);
if strcmp(g.sessType, 'behavior') && strcmp(g.burstfilter, 'none')
    if abs(max(sta(1,absWin))) > abs(min(sta(1,absWin))) % Invert the figure if max is bigger than min
        sta = sta*-1;
        st = st*-1;
        inverted = true;
    end
elseif strcmp(g.sessType, 'CueHit') && strcmp(g.burstfilter, 'none')
    if abs(max(sta(1,absWin))) > abs(min(sta(1,absWin))) % Invert the figure if max is bigger than min
        sta = sta*-1;
        st = st*-1;
        inverted = true;
    end
elseif strcmp(g.sessType, 'CueFA') && strcmp(g.burstfilter, 'none')
    if abs(max(sta(1,absWin))) > abs(min(sta(1,absWin))) % Invert the figure if max is bigger than min
        sta = sta*-1;
        st = st*-1;
        inverted = true;
    end
else
    if ~strcmp(g.cellid2, 'n028_120211a_8.2')
        load('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat');
        inxInv = find(cellfun(@(s)strcmp(s,cellid),cellids));
        invert = invertList(inxInv);
        if invert
            sta = sta*-1;
            st = st*-1;
            inverted = true;
        end
    else
        disp('re-invert')
    end
end
% Save into stats struct
stats.sta_cv = sta_cv;
stats.sta_index1 = sta_index1;
stats.sta_index2 = sta_index2;
stats.sta_amplitude = sta_amp;
stats.spike_number = lenv;

% Randomized STA
nrs = 200;
sta_rand = zeros(nrs,g.wn+1);
wb = waitbar(0,'Generating surrogate data set. Please wait...','Name','Running STACALL...');  % progress indicator
global WB
WB(end+1) = wb;
for rpp = 1:nrs
    vdisc_rand = randpoisson(length(vdisc),length(eeg));
    [sta_rand(rpp,1:g.wn+1), temp_cv, sta_index1_rand(rpp), sta_index2_rand(rpp),...
        sta_amp_rand(rpp)] = astanorm2(round(vdisc_rand),eeg,g.sr,g.wn,isnorm);
    waitbar(rpp/nrs)
end
close(wb)
cellidt = regexprep(cellid,'_',' ');

% ACG file for group identity
acgDir = [RESDIR cellbase '\acg' cellbase PATH '\'];
acgPath = [acgDir 'ACG_matrices_' cellbase '.mat'];
[group, groupID] = cellid2TPgroup(cellid, 'acgPath', acgPath);
stats.group1 = group;
stats.groupID1 = groupID;
if ~isempty(g.cellid2) % In case of pair do the grouping again
    cellidt2 = regexprep(g.cellid2,'_',' ');
    [group2, groupID2] = cellid2TPgroup(g.cellid2, 'acgPath', acgPath);
    stats.group2 = group2;
    stats.groupID2 = groupID2;
    [groupID, ~] = groupCellPair(group, group2);
    if iscell(group2)
        group2 = group2{1};
    end
else
    cellidt2 = '';
    group2 = '';
end
if iscell(group)
    group = group{1};
end

% Store if data was inverted
if inverted
    titlestr = {[cellidt ' - ' cellidt2 ];  [group '-' group2 ' ' g.burstfilter ' ' cellbase ' (Inverted)']};
else
    titlestr = {[cellidt ' - ' cellidt2 ];  [group '-' group2 ' ' g.burstfilter ' ' cellbase]};
end

% Store variables
stats.inverted = inverted;
sta_sd = std(st, [], 1);
sta_sd = sta_sd / sqrt(size(st, 1));
stats.se = sta_sd;
stats.sta_rand = sta_rand;
stats.sta_amp_rand = sta_amp_rand;
stats.sta_index2_rand = sta_index2_rand;
stats.wn = g.wn;
stats.sr = g.sr;

% Plot
[H, lwr, upr] = stafig(sta,sta_rand,sta_amp,sta_amp_rand,sta_index2,...
    sta_index2_rand,lenv,g.wn,g.sr,titlestr, g.sync, sta_sd, groupID);
line([0 0], get(gca, 'ylim'), 'color', 'r');
stats.lower_limit = lwr;
stats.upper_limit = upr;