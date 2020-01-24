function relJittLat_wrapper(cellids)
% update


for iC=1:numel(cellids)
    currCell = cellids{iC};
    cbSwitcher(currCell);
    
    [currEvent, cFilter] = eventSelect(currCell);
    
    % Efficiency, latency and jitter for 'PulseOn'
    [R, L, J, ~, ~, ~, ~, ~, ~] = reliability_latency_jitter(currCell,...
        'event',currEvent, 'event_type', 'trial', 'event_filter', 'custom',...
        'filterinput', cFilter, 'isadaptive', 2,...
        'baselinewin', [-0.02 0], 'testwin', [0 0.1], 'relative_threshold',...
        0.05, 'jitterdefinition', 'burst', 'display', false, 'window',...
        [-1 1]);
    Reliability(iC) = R;
    Latency(iC) = L;
    Jitter(iC) = J;
    
    
end




medLat = nanmedian(Latency);
medRel = nanmedian(Reliability);
medJit = nanmedian(Jitter);

SE_medLat = nanse(Latency);
SE_medRel = nanse(Reliability);
SE_medJit = nanse(Jitter);

keyboard;

%--------------------------------------------------------------------------
function [faevent, cfilter] = eventSelect(cellid)

cb = whichcb;
align = 'feedback';
% Control PSTH and raster alignment
if ~strcmp(cb, 'PannaHDB')
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
            hitevent = findAlignEvent_posfeedback_gonogo(cellid); %#ok<NASGU>
            faevent = findAlignEvent_negfeedback_gonogo(cellid); %#ok<NASGU>
            sevent = 'StimulusOn';
        case 'trialstart'
            hitevent = 'TrialStart';
            faevent = 'TrialStart';
            sevent = 'StimulusOn';
        case 'ITI'
            hitevent = 'ITI';
            faevent = 'ITI';
            sevent = 'StimulusOn';
        case 'Miss'
            hitevent = 'StimulusOn'; %#ok<NASGU>
            faevent = 'StimulusOn'; %#ok<NASGU>
            sevent = 'StimulusOff';
        case 'CorrectRejection'
            hitevent = 'StimulusOn'; %#ok<NASGU>
            faevent = 'StimulusOn'; %#ok<NASGU>
            sevent = 'StimulusOff';
    end
else
    hitevent = 'DeliverAllFeedback';
    faevent = 'DeliverAllFeedback';
    sevent = 'DeliverAllFeedback';
end

trialtype = 'fa';

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