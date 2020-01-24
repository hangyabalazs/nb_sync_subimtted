function [psth, spsth, spsth_se, tags, spt, stats] = ultimate_psth_NBsync(cellid,event_type,event,window,varargin)
%ULTIMATE_PSTH   Peri-stimulus time histogram.
%   [PSTH SPSTH SPSTH_SE] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   calculates peri-stimulus time histogram (PSTH) for the cell passed in
%   CELLID. Smoothed PSTH (SPSTH) and SE of smoothed PSTH (SPSTH_SE) are
%   also returned.
%
%   [PSTH SPSTH SPSTH_SE TAGS] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   returns partition tags (TAGS) corrsponding to PSTHs when trials are
%   partitioned; see PARTITION_TRIALS.
%
%   [PSTH SPSTH SPSTH_SE TAGS SPT] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   returns the bin raster (SPT); see STIMES2BINRASTER.
%
%   [PSTH SPSTH SPSTH_SE SPT TAGS STATS] = ULTIMATE_PSTH(CELLID,EVENT_TYPE,EVENT,WINDOW,VARARGIN)
%   calculates and returns test results for significant firing rate changes
%   after the event (see PSTH_STATS for details).
%
%   ULTIMATE_PSTH is also capable of using two different events for the
%   periods before and after 0, usefull for statistical testing with a
%   baseline period aligned to a different event than the test period (see
%   below and PSTH_STATS).
%
%   Mandatory input arguments:
%       CELLID: defines the cell (see CellBase documentation) or session
%           (for lick-PSTH)
%       EVENT: the event to which the PSTH is aligned; if EVENT is a cell
%           array of two strings, the first event is used for the PSTH
%           and binraster before 0 and the second event is used for the
%           PSTH and binraster after 0; if EVENT is a function handle, the
%           function is called for CELLID to define the aligning event
%           (dynamic event definition)
%       EVENT_TYPE: the type of event, 'stim', 'trial' or 'lick' (for
%           lick-PSTH)
%       WINDOW: window for calculation relative to the event in seconds
%
%   Default behavior of ULTIMATE_PSTH can be modified by using a set of
%   paramter-value pairs as optional input parameters. The following
%   parameters are implemented (with default values):
%   	'dt', 0.001 - time resolution in seconds
%       'sigma', 0.02 - smoothing kernel for the smoothed PSTH, in seconds
%       'margin',[-0.01 0.01] margins for PSTH calculation to get rid of
%           edge effect due to smoothing
%       'event_filter', 'none' - filter light-stimulation trials; see
%           FILTERTRIALS for implemented filter types
%       'filterinput',[] - some filters require additional input; see
%           FILTERTRIALS for details
%       'burst_filter', 'none' - filter for burst first spikes or single
%           spikes - see also BURSTDETECTION.
%       'maxtrialno', 5000 - maximal number of trials included; if ther are
%           more valid trials, they are randomly down-sampled
%       'first_event', [] - event name used to exclude spikes before
%           previous event
%       'last_event', [] - event name used to exclude spikes after
%           following event
%       'parts', 'all' - partitioning the set of trials; input to
%           PARTITION_TRIALS, see details therein (default, no
%           partitioning)
%       'isadaptive', 1 - 0, classic PSTH algorithm is applied; 1, adaptive
%           PSTH is calculated (see APSTH); 2, 'doubly adaptive' PSTH
%           algorithm is used (see DAPSTH)
%   	'baselinewin', [-0.25 0] - limits of baseline window for
%           statistical testing (see PSTH_STATS), time relative to 0 in
%           seconds
%   	'testwin', [0 0.1] - limits of test window for statistical testing
%           (see PSTH_STATS), time relative to 0 in seconds
%       'relative_threshold', 0.5 - threshold used to assess start and end
%           points of activation and inhibition intervals in PSTH_STATS; in
%           proportion of the peak-baseline difference (see PSTH_STATS)
%       'display', false - controls plotting
%
%   See also PSTH_STATS, STIMES2BINRASTER, BINRASTER2PSTH, BINRASTER2APSTH,
%   APSTH, VIEWCELL2B, PARTITION_TRIALS and FILTERTRIALS.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   07-May-2012

%   Edit log: BH 7/5/12, 8/12/12, 8/27/12, 12/12/13, 5/20/14

% Default arguments
prs = inputParser;
addRequired(prs,'cellid',@(s)iscellid(s)|issessionid(s))
addRequired(prs,'event_type',@ischar)   % event type ('stim' or 'trial')
addRequired(prs,'event',@(s)ischar(s)|...
    (iscellstr(s)&isequal(length(s),2))|...
    isa(s,'function_handle'))   % reference event
addRequired(prs,'window',@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event, in seconds
addParameter(prs,'event_filter','none',@(s)ischar(s)|iscellstr(s))   % filter events based on properties
addParameter(prs,'filterinput',[])   % some filters need additional input
addParameter(prs,'burst_filter','none',@(s)ischar(s)&ismember(s,{'none','burst1','single', 'burstall', 'burst1Single'}))   % filter on burst first spikes or single spikes
addParameter(prs,'maxtrialno',5000)   % downsample events if more than 'maxtrialno'
addParameter(prs,'first_event','',@(s)isempty(s)|ischar(s))   % exclude spikes before previous event
addParameter(prs,'last_event','',@(s)isempty(s)|ischar(s))   % exclude spikes after following events
addParameter(prs,'dt',0.001,@isnumeric)   % time resolution of the binraster, in seconds
addParameter(prs,'sigma',0.02,@isnumeric)     % smoothing kernel for the smoothed PSTH
addParameter(prs,'margin',[-0.1 0.1])  % margins for PSTH calculation to get rid of edge effect due to smoothing
addParameter(prs,'parts','all')   % partition trials
addParameter(prs,'isadaptive',true,@(s)islogical(s)|ismember(s,[0 1 2]))   % use adaptive PSTH algorithm
addParameter(prs,'baselinewin',[-0.25 0],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addParameter(prs,'testwin',[0 0.1],@(s)isnumeric(s)&isequal(length(s),2))  % time window relative to the event for stat. testing, in seconds
addParameter(prs,'relative_threshold',0.5,@(s)isnumeric(s)&s>=-1&s<=1)   % threshold used to assess interval limits in PSTH_STATS; negative thresholds selects the full window
addParameter(prs,'display',false,@(s)islogical(s)|ismember(s,[0 1]))   % control displaying rasters and PSTHs
addParameter(prs,'cellid2', [], @(s)iscell(s)|isempty(s)|iscellstr(s)|ischar(s))
addParameter(prs,'sync',false,@islogical)   % control filter for sync events
addParameter(prs,'downsample',false,@islogical)   % control filter for sync events
parse(prs,cellid,event_type,event,window,varargin{:})
g = prs.Results;
if nargout > 5   % statistics will only be calculted if asked for
    g.dostats = true;
else
    g.dostats = false;
end
global RESDIR;
cellbase = whichcb;
% Load event structure
if ~strcmp(event_type, 'trialstart')
    event_type = lower(event_type(1:4));
end
switch event_type
    case 'stim'
        % Load stimulus events
        try
            VE = loadcb(cellid,'StimEvents');   % load events
            VS = loadcb(cellid,'STIMSPIKES');   % load prealigned spikes
        catch ME
            disp('There was no stim protocol for this session.')
            error(ME.message)
        end
    case 'tria'
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');   % load events
            VS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
            if g.sync % PSTH for sync events
                VS2 = loadcb(g.cellid2, 'EVENTSPIKES');
            end
        catch ME
            disp('There was no behavioral protocol for this session.')
            error(ME.message)
        end
    case 'trialstart'
        % Load trial events
        try
            VE = loadcb(cellid,'TrialEvents');   % load events
            VS = loadcb(cellid,'EVENTSPIKES');   % load prealigned spikes
            if g.sync % PSTH for sync events
                VS2 = loadcb(g.cellid2, 'EVENTSPIKES');
            end
        catch ME
            disp('There was no behavioral protocol for this session.')
            error(ME.message)
        end
    case 'lick'
        % Load trial events (unsynchronized)
        animalID = cellid{1,1};   % session ID
        sessionID = cellid{1,2};
        cbdir = getpref('cellbase','datapath');
        datapath = fullfile(cbdir,animalID,sessionID,'TE.mat');
        VE = load(datapath);
    otherwise
        error('Input argument ''event_type'' should be either ''stim'' or ''trial''.')
end

% Events, time, valid trials
time = (window(1)+g.margin(1)):g.dt:(window(2)+g.margin(2));
if iscell(event)   % different ref. event for test and baseline window
    event1 = event{1};   % baseline event
    event2 = event{2};   % test event
    
    switch event_type
        case {'stim','tria'}
            event_pos1 = findcellstr(VS.events(:,1),event1);
            event_pos2 = findcellstr(VS.events(:,1),event2);
            if event_pos1 * event_pos2 == 0
                error('Event name not found');
            end
            stimes1 = VS.event_stimes{event_pos1};   % baseline spike times
            stimes2 = VS.event_stimes{event_pos2};   % test spike times
            if g.sync
                stimes3 = VS2.event_stimes{event_pos1};
                stimes4 = VS2.event_stimes{event_pos2};
                for i = 1:size(stimes1,2)
                    stimes1{i} = sort(cat(1, stimes1{i}, stimes3{i}), 'ascend');
                    stimes2{i} = sort(cat(1, stimes2{i}, stimes4{i}), 'ascend');
                end
            end
            triggerevent1 = VS.events{event_pos1,2};   % trigger event for baseline (event name may differ)
            triggerevent2 = VS.events{event_pos2,2};   % trigger event for test (event name may differ)
            valid_trials1 = filterTrials(cellid,'event_type',event_type,'event',event1,...
                'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
            valid_trials2 = filterTrials(cellid,'event_type',event_type,'event',event2,...
                'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
            if ~isequal(valid_trials1,valid_trials2)
                error('Valid trials should be the same for baseline and test period.')
            else
                valid_trials = valid_trials1;
            end
        case 'lick'
            NumTrials = length(VE.LickIn);
            stimes1 = arrayfun(@(s)VE.LickIn{s}-VE.(event1)(s),1:NumTrials,'UniformOutput',false);   % lick times
            stimes2 = arrayfun(@(s)VE.LickIn{s}-VE.(event2)(s),1:NumTrials,'UniformOutput',false);
            triggerevent1 = event1;
            triggerevent2 = event2;
            if ~isequal(g.event_filter,'none')
                error('Trial filters not implemented for lick PSTH.')
            end
            valid_trials = 1:NumTrials;
        case 'trialstart'
            NumTrials = length(VE.TrialStart);
            stimes1 = arrayfun(@(s)VE.TrialStart{s}-VE.(event1)(s),1:NumTrials,'UniformOutput',false);   % lick times
            stimes2 = arrayfun(@(s)VE.TrialStart{s}-VE.(event2)(s),1:NumTrials,'UniformOutput',false);
            triggerevent1 = event1;
            triggerevent2 = event2;
            if ~isequal(g.event_filter,'none')
                error('Trial filters not implemented for lick PSTH.')
            end
            valid_trials = 1:NumTrials;
    end
    [stimes1, starttimes1, endtimes1] = dynamicSpikeWindow(stimes1,VE,triggerevent1,g.first_event,g.last_event);   % restrict between previous and following event
    [stimes2, starttimes2, endtimes2] = dynamicSpikeWindow(stimes2,VE,triggerevent2,g.first_event,g.last_event);   % restrict between previous and following event
    
else
    if isa(event,'function_handle')
        event = feval(event,cellid);   % dynamic event definition
    end
    
    switch event_type
        case {'stim','tria'}
            event_pos = findcellstr(VS.events(:,1),event);
            if event_pos == 0
                error('Event name not found');
                psth = [];
                stats = [];
            end
            switch g.burst_filter
                case 'none'
                    stimes = VS.event_stimes{event_pos};
                    filtName = 'event_stimes';
                case 'burst1'
                    [stimes, ~, ~] = burstDetect(VS.event_stimes{event_pos});
                    filtName = 'event_bursttimes1';
                case 'single'
                    [~, stimes, ~] = burstDetect(VS.event_stimes{event_pos});
                    filtName = 'event_singlespikes';
                case 'burstall'
                    [~, ~, stimes] = burstDetect(VS.event_stimes{event_pos});
                    filtName = 'event_bursttimesall';
                case 'burst1Single' % Concatenate burst1 and Single spikes for the two cells
                    [~, sstimes, ~] = burstDetect(VS.event_stimes{event_pos});
                    [b1stimes, ~, ~] = burstDetect(VS.event_stimes{event_pos});
                    for i = 1:size(sstimes,2)
                        stimes{i} = sort(cat(1, sstimes{i}, b1stimes{i}), 'ascend');
                    end
                    filtName = 'event_burst1Single';
                    if g.sync % In case of sync events
                        [~, ss2times, ~] = burstDetect(VS2.event_stimes{event_pos});
                        [b12stimes, ~, ~] = burstDetect(VS2.event_stimes{event_pos});
                        for i = 1:size(sstimes,2)
                            stimes2{i} = sort(cat(1, ss2times{i}, b12stimes{i}), 'ascend');
                        end
                    end                    
            end
            % Loop for PSTH sync events
            if g.sync
                if ~strcmp(g.burst_filter, 'burst1Single')
                    stimes2 = VS2.(filtName){event_pos};
                end
                for i = 1:size(stimes,2)
                    if isempty(stimes{i}) && isempty(stimes2{i})
                        disp('Empty event for both cells')
                    else
                        Sstimes = sort(cat(1, stimes{i}, stimes2{i}), 'ascend'); %Concatenate the two cells
                        syncEvent = diff(Sstimes)<=0.01; % Sync events in a 10ms window
                        stimes{i} = Sstimes(syncEvent);
                        if g.downsample
                            downSamp = sum(double(syncEvent)); % Number of sync events
                            if downSamp == 0
                                stimes(i)={[]};
                            else
                                Sync1 = find(syncEvent==1)+1; % Second members of the sync events
                                syncEvent(Sync1) = 1; % Mark them also as part of sync events
                                nonSyncEvent = ~syncEvent;
                                stimes{i} = Sstimes(nonSyncEvent); % Take all the ones which are not part of a sync event
                            end
                        end
                    end
                end
            end
            triggerevent = VS.events{event_pos,2};   % trigger event (event name may differ)
            valid_trials = filterTrials(cellid,'event_type',event_type,'event',event,...
                'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
        case 'trialstart'
            NumTrials = length(VE.TrialStart);
            stimes = VE.TrialStart;
            triggerevent = event;
            %             if ~isequal(g.event_filter,'none')
            %                 error('Trial filters not implemented for lick PSTH.')
            %             end
            %             valid_trials = 1:NumTrials;
            valid_trials = filterTrials(cellid,'event_type',event_type,'event',event,...
                'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
        case 'lick'
            NumTrials = length(VE.LickIn);
            stimes = arrayfun(@(s)VE.LickIn{s}-VE.(event)(s),1:NumTrials,'UniformOutput',false);   % lick times
            triggerevent = event;
            %             if ~isequal(g.event_filter,'none')
            %                 error('Trial filters not implemented for lick PSTH.')
            %             end
            %             valid_trials = 1:NumTrials;
            valid_trials = filterTrials(cellid,'event_type',event_type,'event',event,...
                'event_filter',g.event_filter,'filterinput',g.filterinput);   % filter trials
     
    end
    [stimes, starttimes, endtimes] = dynamicSpikeWindow(stimes,VE,triggerevent,g.first_event,g.last_event);   % restrict between previous and following event
end

% Calculate bin rasters
if iscell(event)   % different ref. event for test and baseline window
    spt1 = stimes2binraster(stimes1,time,g.dt);
    spt1 = nanpadspt(time,spt1,starttimes1,endtimes1);  % replace zeros with NaNs outside the dynamic window
    spt2 = stimes2binraster(stimes2,time,g.dt);
    spt2 = nanpadspt(time,spt2,starttimes2,endtimes2);  % replace zeros with NaNs outside the dynamic window
    spt = [spt1(:,time<0) spt2(:,time>=0)];   % merge baseline and test raster
else
    spt = stimes2binraster(stimes,time,g.dt);
    spt = nanpadspt(time,spt,starttimes,endtimes);  % replace zeros with NaNs outside the dynamic window
end

% Partition trials
[COMPTRIALS, tags] = partition_trials(VE,g.parts);

% PSTH
switch g.isadaptive
    case {0,false}
        [psth, spsth, spsth_se] = binraster2psth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
    case {1, true}
        [psth, spsth, spsth_se] = binraster2apsth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
    case 2
        [psth, spsth, spsth_se] = binraster2dapsth(spt,g.dt,g.sigma,COMPTRIALS,valid_trials);
end
stm0 = abs(window(1)+g.margin(1)) * (1 / g.dt);   % zero point (diveding directly with 'g.dt' can result in numeric deviation from the desired integer result)
stm = round(stm0);   % still numeric issues
if abs(stm-stm0) > 1e-10
    error('Zero point is not an integer.')
end
inx = (stm+1+window(1)/g.dt):(stm+1+window(2)/g.dt);     % indices for cutting margin
psth = psth(:,inx);
spsth = spsth(:,inx);
spsth_se = spsth_se(:,inx);
NumPartitions = size(psth,1);
if NumPartitions > 1   % partitions
    pspt = spt;
    spt = cell(1,NumPartitions);
    for iP = 1:NumPartitions
        spt{iP} = pspt(intersect(valid_trials,COMPTRIALS{iP}),inx);
    end
else
    spt = spt(valid_trials,inx);
end

% Output statistics
global PATH;
if g.dostats
    stats = psth_stats(spt,psth,g.dt,window,...
        'baselinewin',g.baselinewin,'testwin',g.testwin,'display',g.display,...
        'relative_threshold',g.relative_threshold);
    statsS = psth_stats(spt,spsth,g.dt,window,...
        'baselinewin',g.baselinewin,'testwin',g.testwin,'display',g.display,...
        'relative_threshold',g.relative_threshold);
    if g.sync
        group1 = cellid2TPgroup(cellid, 'acgPath', [RESDIR cellbase '\acg' cellbase PATH '\ACG_matrices_' cellbase '.mat']);
        group2 = cellid2TPgroup(g.cellid2, 'acgPath', [RESDIR cellbase '\acg' cellbase PATH '\ACG_matrices_' cellbase '.mat']);
        [~, group] = groupCellPair(group1, group2);
        stats.group = group;
        statsS.group = group;
    else
        try
            [~, groupID] = cellid2TPgroup(cellid, 'acgPath', [RESDIR cellbase '\acg' cellbase PATH '\ACG_matrices_' cellbase '.mat']);
            stats.group = groupID;
            statsS.group = groupID;
        catch
            disp('GROUP ID IS MISSING')
        end
    end
    stats.SPSTH = statsS;
end

% -------------------------------------------------------------------------
function spt = nanpadspt(time,spt,starttimes,endtimes)

% For variable windows, change padding to NaN to ensure correct averaging
NUMtrials = size(spt,1);   % number of trials
for iT = 1:NUMtrials    % loop through trials
    inx = time < starttimes(iT);
    spt(iT,inx) = NaN;   % NaN trials before previous event
end
for iT = 1:NUMtrials    % loop through trials
    inx = time > endtimes(iT);
    spt(iT,inx) = NaN;   % NaN trials after following event
end