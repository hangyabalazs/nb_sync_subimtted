function psthExtract(I, issave, revCells, cb)
%PSTHEXTRACT   Runs ultimate psth for all the filters, with all smoothing.
%   PSTHEXTRACT(CELLIDS, ISSAVE) Psth wrapper for ultimate_psth
%
%   Output:
%       Eventdata.mat stores all data for every trialtype
%   See also ULTIMATE_PSTH_NBSYNC, EVENTDATAPLOT_ALL.

%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Pass the control to the user in case of error
dbstop if error;
global RESDIR;
global PATH;

% Input argument check
narginchk(0,4)
if nargin < 4
    issave = false;   % saving is controlled by the second input argument
    revCells = false;
    cb = whichcb;
end

if nargin < 1
    % Select cholinergic cells (n=45)
    [ChAT, pChAT, allChAT] = selectChAT(cb);
    I = allChAT;
end

% List of cellIDs
I = I(:)';   % convert to row vector

% Cells which hasn't got burst in the activation window
if revCells
    % Load removecell list
    rName = 'D:\_MATLAB_DATA\POOLED\psthPOOLED\removeCells\removeCells.mat';
    if exist(rName, 'file')
        rCells = load('D:\_MATLAB_DATA\POOLED\psthPOOLED\removeCells\removeCells.mat');
        revInx = cellfun(@(x) any(contains(rCells.nonActCells, x)), I);
        I = I(~revInx);
        if any(revInx)
            revFile = 'removed';
        else
            revFile = [];
        end
    else
        disp('NO REMOVECELLS FILE FOUND')
        revCells = false;
        keyboard;
    end
else
    revFile = [];
end

% Control 'poppsth' behavior
doraster = false;

align1 = 'tone';  % controls alignment: 'tone' or 'response' or 'feedback
align2 = 'feedback';
align3 = 'trialstart';
align5 = 'Miss';
align6 = 'CorrectRejection';

% Directories
fs = filesep;
resdir = [RESDIR cb fs 'PSTH_newdata' cb PATH fs];   % results directory

% PSTH
% Trialstart
[baseline_trialstart, maxvalue_ts, EVENTDATA_TRIALSTART_FA] = poppsth(I,...
    resdir,align3,'fa',doraster,issave); %#ok<*ASGLU>
disp('TRIALSTART FA DONE')
[baseline_trialstart, maxvalue_ts, EVENTDATA_TRIALSTART_HIT] = poppsth(I,...
    resdir,align3,'hit',doraster,issave); %#ok<*ASGLU>
disp('TRIALSTART HIT DONE')

% Cue
[baseline_tone, maxvalue_tone, EVENTDATA_TONE_FA] = poppsth(I,resdir,...
    align1,'fa',doraster,issave); %#ok<*ASGLU>
disp('TONE FA DONE')
[baseline_tone, maxvalue_tone, EVENTDATA_TONE_HIT] = poppsth(I,resdir,...
    align1,'hit',doraster,issave); %#ok<*ASGLU>
disp('TONE HIT DONE')
% Feedback
[baseline_fb, maxvalue_fb, EVENTDATA_FB_FA] = poppsth(I,resdir,align2,...
    'fa',doraster,issave); %#ok<*ASGLU>
disp('FEEDBACK FA DONE')
[baseline_fb, maxvalue_fb, EVENTDATA_FB_HIT] = poppsth(I,resdir,align2,...
    'hit',doraster,issave); %#ok<*ASGLU>
disp('FEEDBACK HIT DONE')
% Miss
[baseline_miss, maxvalue_miss, EVENTDATA_MISS] = poppsth(I,resdir,align5,...
    'fa',doraster,issave); %#ok<*ASGLU>
disp('MISS DONE')
% CorrectRejection
[baseline_cr, maxvalue_cr, EVENTDATA_CR] = poppsth(I,resdir,align6,'hit',...
    doraster,issave); %#ok<*ASGLU>
disp('CR DONE')

% Save results into struct
EventData.FB.FA = EVENTDATA_FB_FA.feedback.fa;
EventData.FB.HIT = EVENTDATA_FB_HIT.feedback.hit;
EventData.TONE.FA = EVENTDATA_TONE_FA.tone.fa;
EventData.TONE.HIT = EVENTDATA_TONE_HIT.tone.hit;
EventData.TRIALSTART.FA = EVENTDATA_TRIALSTART_FA.trialstart.fa;
EventData.TRIALSTART.HIT = EVENTDATA_TRIALSTART_HIT.trialstart.hit;
EventData.MISS = EVENTDATA_MISS;
EventData.CR = EVENTDATA_CR;

% Save struct
save([resdir 'EventData_review_' revFile '_' date '.mat'], 'EventData',...
    '-v7.3');
if issave
    fnm = fullfile(resdir,['stats' revFile '.mat']);
    save(fnm,'baseline_fb','maxvalue_fb');
end
keyboard;

% -------------------------------------------------------------------------
function [baseline, maxvalue, EVENTDATA] = poppsth(I,resdir,align,...
    trialtype,doraster,issave)

dbstop if error;
% Load CellBase if indices to CELLIDLIST are passed
if isnumeric(I)
    loadcb
    I = CELLIDLIST(I);
end

% Time window
wn = [-1 1];   % in seconds
dt = 0.001;

% Call 'ultimate_psth'
% Preallocate
problem_ids = {};   % cellIDs for which PSTH failed
problem_msg = {};   % corresponding error messages
[baseline, maxvalue] = deal([]);  % PSTH stats
NumCell = length(I);  % number of cells
% adapFilter = {'NotAdaptive', 'Adaptive', 'DoublyAdaptive'};
adapFilter = {'NotAdaptive'};
burstFilter = {'none', 'burst1', 'single', 'burstall', 'burst1Single'};

for i = 1:length(adapFilter)
    for l = 1:length(burstFilter)
        for k = 1:NumCell
            % Special case for CPu and Acb cholinergic cells
            cellid = I{k};
            disp(cellid)
            cbSwitcher(cellid);
            cb = whichcb;
            
            % Control PSTH and raster alignment
            switch align
                case 'tone'
                    hitevent = 'StimulusOn'; %#ok<NASGU>
                    faevent = 'StimulusOn'; %#ok<NASGU>
                    sevent = 'StimulusOff';
                case 'response'
                    hitevent = 'LeftWaterValveOn'; %#ok<NASGU>
                    faevent = 'LeftPortIn'; %#ok<NASGU>
                    sevent = 'StimulusOn';
                case 'feedback'
                    if ~strcmp(cb, 'PannaHDB')
                        hitevent = findAlignEvent_posfeedback_gonogo(cellid); %#ok<NASGU>
                        faevent = findAlignEvent_negfeedback_gonogo(cellid); %#ok<NASGU>
                        sevent = 'StimulusOn';
                    else
                        hitevent = 'DeliverAllFeedback';
                        faevent = 'DeliverAllFeedback';
                        sevent = 'DeliverAllFeedback';
                    end
                case 'trialstart'
                    hitevent = 'TrialStart';
                    faevent = 'TrialStart';
                    sevent = 'TrialStart';
                case 'Miss'
                    hitevent = 'StimulusOn'; %#ok<NASGU>
                    faevent = 'StimulusOn'; %#ok<NASGU>
                    sevent = 'StimulusOff';
                case 'CorrectRejection'
                    hitevent = 'StimulusOn'; %#ok<NASGU>
                    faevent = 'StimulusOn'; %#ok<NASGU>
                    sevent = 'StimulusOff';
            end
            cevent = eval([trialtype 'event']);  % hitevent or faevent
            if strcmp(align, 'Miss') || strcmp(align, 'CorrectRejection')
                hitfilter = 'CorrectRejection==1'; %#ok<NASGU>
                fafilter = 'Miss==1'; %#ok<NASGU>
            else
                fafilter = 'FalseAlarm==1'; %#ok<NASGU>
                hitfilter = 'Hit==1'; %#ok<NASGU>
            end
            cfilter = eval([trialtype 'filter']);   % hitfilter or fafilter
            if strcmp(cb, 'NB')
                switch trialtype
                    case 'hit'
                        tonelastev = findAlignEvent_posfeedback_gonogo(cellid); %#ok<NASGU>
                    case 'fa'
                        tonelastev = findAlignEvent_negfeedback_gonogo(cellid); %#ok<NASGU>
                end
            else
                tonelastev = [];
            end
            feedbacklastev = []; %#ok<NASGU>
            if strcmp(align, 'trialstart') || strcmp(align, 'ITI')...
                    || strcmp(align, 'Miss') || strcmp(align, 'CorrectRejection')...
                    || strcmp(align, 'punishment')
                clastev = [];
            else
                clastev = eval([align 'lastev']);   % last event is feedback for tone alignment
            end
            
            switch adapFilter{i}
                case 'NotAdaptive'
                    adapFilt = 0;
                case 'Adaptive'
                    adapFilt = 1;
                case 'DoublyAdaptive'
                    adapFilt = 2;
            end
            trType = 'trial';
            try
                % Calcualte PSTH for current parameters
                [psth, spsth, spsth_se, tags, spt, stats] = ultimate_psth_NBsync(cellid,trType,cevent,wn,...
                    'dt',dt,'display',false,'sigma',0.002,'parts','all',...
                    'isadaptive',adapFilt,'event_filter','custom','filterinput',cfilter,...
                    'first_event',[],'last_event',clastev,...
                    'maxtrialno',Inf,'baselinewin',[-0.5 0],'testwin',[0 0.25],...
                    'relative_threshold',0.1, 'burst_filter', burstFilter{l});
                
                % Store data in struct
                EVENTDATA.(align).(trialtype).(adapFilter{i}).(burstFilter{l}).psth(k,:) = psth;
                EVENTDATA.(align).(trialtype).(adapFilter{i}).(burstFilter{l}).spsth(k,:) = spsth;
                EVENTDATA.(align).(trialtype).(adapFilter{i}).(burstFilter{l}).spsth_se(k,:) = spsth_se;
                EVENTDATA.(align).(trialtype).(adapFilter{i}).(burstFilter{l}).raster{k} = spt;
                EVENTDATA.(align).(trialtype).(adapFilter{i}).(burstFilter{l}).stats{k} = stats;
                EVENTDATA.(align).(trialtype).(adapFilter{i}).(burstFilter{l}).cellid{k} = cellid;
                baseline(k) = stats.baseline;
                maxvalue(k) = stats.maxvalue;
                
                % Raster plot
                if doraster
                    H = figure;
                    viewcell2b(cellid,'TriggerName',cevent,'SortEvent',sevent,...
                        'eventtype','behav','ShowEvents',{{sevent}},...
                        'LastEvents',clastev,'Partitions','#ResponseType','window',wn)
                    maximize_figure(H)
                    
                    % Save raster plot figure
                    if issave
                        fnm = [resdir cellidt '_raster_CR_MISS_' align '.fig'];
                        saveas(H,fnm)
                        fnm = [resdir cellidt '_raster_CR_MISS_' align '.jpg'];
                        saveas(H,fnm)
                    end
                    close(H)
                end
            catch ME
                % Error handling
                problem_ids{end+1} = cellid;  %#ok<AGROW> % collect cellIDs resulting in error
                problem_msg{end+1} = ME.message;  %#ok<AGROW> % store corresponding error messages
                disp(['Something went wrong for cell ' cellid '.'])
                disp(ME.message)   % display error message
            end
        end
    end
end