function sta_ers_fig9(varargin)
%STA_ERS   Spike triggered average, event related spectrogram
%   STA_ERS(SESSTYPE, SYNC, BURSTFILTER) STA and spectrum for LFP data.
%   It generates individual sta plots and spectrums. Data saved to struct
%   for group analysis.
%
%   See also STACALL, STAERS, FIVEBAR, ERS, EEGWAVELET

% Input arguments
prs = inputParser;
addParameter(prs,'sessType','behavior',@(s)ischar(s)&ismember(s,{'behavior',...
    'burstOn', 'stimFA', 'stimHIT', 'CR', 'Miss'}))   % session type filter
addParameter(prs,'sync','original',...
    @(s)ismember(s,{'original','sync','nonsync'}))   % sync filter
addParameter(prs,'burstfilter','none',@(s)ischar(s)&ismember(s,{'none',...
    'burst1', 'burstall', 'single', 'burst1Single'}))   % burstfilter
parse(prs,varargin{:})
g = prs.Results;

% Pass the control to the user in case of error
dbstop if error;
cb = whichcb;
fs = filesep;
global RESDIR;



% Initialize current cellids
[ChAT, pChAT, allChAT] = selectChAT(cb);
NumChAT = length(allChAT);
Cellids = allChAT;
%Check input datatypes, conversion to cell
[sessType, burstfilter, sync] = cellChecker(g.sessType, g.burstfilter, g.sync);

% % Numver options
% if strcmp(g.sessType, 'behavior')
%     switch g.burstfilter
%         case 'none'
%             numVer = '1'; % RAW
%         case 'burst1' 
%             numVer = '2'; % BURST
%         case 'single'
%             numVer = '3'; % SINGLE
%         case 'burst1Single'
%             if strcmp(g.sync, 'sync')
%                 numVer = '4';
%             else
%                 numVer = '5';
%             end
%     end
% else
%     switch g.sessType
%         case 'stimHIT'
%             numVer = '6';
%         case 'stimFA'
%             numVer = '7';
%         case 'CR'
%             numVer = '8';
%         case 'Miss'
%             numVer = '9';
%         case 'burstOn'
%             numVer = '10';
%     end
% end
numVer = '13';

for si = 1:length(sessType) % Sessiontype loop
    for bi = 1:length(burstfilter) % Burstfilter loop
        burstfilter = burstfilter{bi};
        for iS = 1:length(sync) % Sync loop
            sMode = sync{iS};
            % STA loop
            if ~strcmp('original', sMode) % Sync, NonSync loop
                % Load ccg pairs for sync calculations
                load([RESDIR cb fs 'ccg' cb fs sessType{si} fs burstfilter cb fs...
                    'CCG_matrices_' sessType{si} '_' burstfilter '_' cb '.mat']);
                cellidS = PairOfCells(:,1);
                cellidS2 = PairOfCells(:,2);
                for iC = 1:length(cellidS)
                    cellid = cellidS{iC};
                    if strcmp(cellid(3), num2str(7)) % no LFP for these cells
                        disp('No LFP for the current cell')
                        if strcmp(numVer, '1')
                            load('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat');
                            invertList(iC) = false;
                            cellids{iC} = cellid;
                            save('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat', 'invertList', 'cellids');
                        end
                    else
                        stamain(cellid, burstfilter, sessType{si},...
                            cellidS2{iC}, iC, sMode, numVer);   % calculate STA in a loop
                    end
                    disp(['Cell: ' num2str(iC)])
                end
            else % Original loop
                for iC = 1:NumChAT
                    cellid = Cellids{iC};
                    if strcmp(cellid(3), num2str(7))
                        disp('No LFP for the current cell')
                        if strcmp(numVer, '1')
                            load('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat');
                            invertList(iC) = false;
                            cellids{iC} = cellid;
                            save('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat', 'invertList', 'cellids');
                        end
                    else
                        stamain(cellid, burstfilter, sessType{si},...
                            '', iC, sMode, numVer);   % calculate STA in a loop
                    end
                    disp(['Cell: ' num2str(iC)])
                end
            end
            disp([sMode ' done'])
        end
        disp([burstfilter ' done'])
    end
end

% -------------------------------------------------------------------------
function stamain(cellid, burstfilter, sessType,cellid2, index, sMode, numVer)

cellbase = whichcb;
global RESDIR;
global PATH;
resdir = [RESDIR cellbase '\sta' cellbase '\' burstfilter cellbase ];  % results directory
global cCode;

% LFP full path
[r, s, tetrodename] = cellid2tags(cellid);   %#ok<*ASGLU> % get tetrode number
matname = cellid2fnames(cellid,'cont',tetrodename);   % .mat filename for LFP
[pathname, filename, extension] = fileparts(matname);   % parse filename

% for n037 only CSC8 containd LFP data
if strcmp('n037', r)
    filename = 'CSC8';
else
    filename = 'CSC5';
end

cscname = fullfile(pathname,[filename '.ncs']);   % filename for the Neuralynx CSC file

% Load LFP
try
[TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples, NlxHeader] = ...
    Nlx2MatCSC(cscname,[1 1 1 1 1],1,1,1);  % A1 LFP
catch
   disp('No CSC data for current cell!')
   return
end
lfp_orig = -1 * Samples(:);   % invert Neuralynx data
sr_orig = 512 / mean(diff(TimeStamp)) * 10^6;    % original sampling rate (changed later)
dt_orig = 1 / sr_orig;   % original timestep
starttime = TimeStamp(1) / 10^6;   % start time of recording
time_orig = repmat(TimeStamp/10^6,512,1) + repmat((0:511)'*dt_orig,1,length(TimeStamp));
time_orig = time_orig(:)';

% Load Events
fn = [pathname '\EVENTS.mat'];
load(fn);         % load converted Neuralynx events

% Use the longest recording from file
start_recording = find(strcmp(Events_EventStrings,'Starting Recording'));   % find recording start time(s)
if length(start_recording) > 1   % more than one starts
    disp([cellid ': Multiple recording sessions found. Longest recording is used.'])
    endinx = [start_recording(2:end)-1 length(Events_EventIDs)];   % last event of each recording
    lenrec = Events_TimeStamps(endinx) - Events_TimeStamps(start_recording);  % length of the recording
    chinx = lenrec==max(lenrec);   % index for the longest recording
    pst = Events_TimeStamps(start_recording(chinx));   % use longest recording
    starttime = time_orig(find(time_orig>=pst,1,'first'));   % start time of last recording
    linx = time_orig > pst;  % indices of last recording
    lfp_orig = lfp_orig(linx);   % restrict LFP
    time_orig = time_orig(linx);   % restrict time vector
end

% Resample
sr = 1000;  % new sampling rate
dt = 1 / sr;   % new time step
[p q] = rat(sr/sr_orig);    % resample LFP at 1000 Hz
lfp = resample(lfp_orig,p,q);
time = (0:length(lfp)-1) * dt + starttime;   % generate new time vector

switch sessType
    case 'behavior'
        % Determine time window
        wn = 1000;    % +/- 0.5s window
    case 'burstOn'
        % Determine time window
        wn = 10000;    % +/- 5s window
    otherwise
        % Determine time window
        wn = 1000;    % +/- 0.5s window
end
% Load spike train
limit_spikes = [0 50000];   % include max 100000 spikes

try
    switch sessType
        case 'behavior' % behavior
            tseg = findSegs3(cellid,'segfilter','stim_excl_nb');  % find time segments
            unit = extractSegSpikes_NBsync(cellid,tseg, 'burstfilter', burstfilter);   % find spikes in the time segments
        case 'stimHIT' % behavior
            tseg = findSegs3(cellid,'segfilter','hit_incl_nb',...
                'feedback_duration', [-0.5 0.5],'margins',[0 0]);  % find time segments
            unit = extractSegSpikes_NBsync(cellid,tseg, 'burstfilter', burstfilter);   % find spikes in the time segments
        case 'stimFA' % behavior
            tseg = findSegs3(cellid,'segfilter','fa_incl_nb',...
                'feedback_duration', [-0.5 0.5],'margins',[0 0]);  % find time segments
            unit = extractSegSpikes_NBsync(cellid,tseg, 'burstfilter', burstfilter);   % find spikes in the time segments
        case 'CR' % behavior
            tseg = findSegs3(cellid,'segfilter','cr_incl_nb',...
                'feedback_duration', [-0.5 0.5],'margins',[0 0]);  % find time segments
            unit = extractSegSpikes_NBsync(cellid,tseg, 'burstfilter', burstfilter);   % find spikes in the time segments
        case 'Miss' % behavior
            tseg = findSegs3(cellid,'segfilter','miss_incl_nb',...
                'feedback_duration', [-0.5 0.5],'margins',[0 0]);  % find time segments
            unit = extractSegSpikes_NBsync(cellid,tseg, 'burstfilter', burstfilter);   % find spikes in the time segments
        case 'ITI' % ITI
            tseg = findSegs3(cellid,'segfilter','ITI');  % find time segments
            unit = extractSegSpikes_NBsync(cellid,tseg, 'burstfilter', burstfilter);   % find spikes in the time segments
        case 'burstOn' % burstOn
            try
                % Load stimulus events
                SE = loadcb(cellid,'StimEvents');
                % Pre-stimulus segments: from the offset of a pulse to the onset of the next
                if ~isempty(SE.BurstOn)
                    unit = [SE.BurstOn];
                    unit = unit(~isnan(unit));
                else
                    unit = [];
                end
            catch
                disp([cellid ': There was no stim protocol for this session.'])
                unit = [];
            end
        case 'burstOff' % burstOff
            try
                % Load stimulus events
                SE = loadcb(cellid,'StimEvents');
                % Pre-stimulus segments: from the offset of a pulse to the onset of the next
                if ~isempty(SE.BurstOff)
                    unit = [SE.BurstOff];
                else
                    unit = [];
                end
            catch
                disp([cellid ': There was no stim protocol for this session.'])
                unit = [];
            end
    end
catch ME
    disp('findSegs3 exited with the following error.')
    disp(ME.message)
    return;
end

if any(unit>time(end)) | any(unit<time(1))   %#ok<OR2> % drop spikes out of LFP time range
    noso = sum(unit>time(end)|unit<time(1));   % number of out-of-range spikes
    unit(unit>time(end)|unit<time(1)) = [];   % restrict to LFP time range
    disp([cellid ': ' num2str(noso) ' spikes out of LFP time range.'])
end

if ~isempty(cellid2) % calculates sync events for cellpairs from the same session
    sync = true;
    unit2 = extractSegSpikes_NBsync(cellid2,tseg, 'burstfilter', burstfilter);   % find spikes in the time segments
    try
        unitAll = [unit; unit2];
    catch
        unitAll = [unit'; unit2];
    end
    unitAll = sort(unitAll, 'ascend');
    switch sMode
        case 'sync'
            syncEvent = diff(unitAll)<=0.01; % Sync events in a 10ms window
            unit = unitAll(syncEvent);
        case 'nonsync'
            syncEvent = diff(unitAll)<=0.01;
            Sync1 = find(syncEvent==1)+1; % Second members of the sync events
            syncEvent(Sync1) = 1; % Mark them also as part of sync events
            nonSyncEvent = ~syncEvent;
            unit = unitAll(nonSyncEvent); % Take all the ones which are not part of a sync event
            inx = randperm(length(unit));
            mixUnit = unit(inx); % Randomly shuffle the timestamps
            unit = mixUnit;
    end
else
    sync = false;
end

if length(unit) > limit_spikes(2)      % crop if too long to avoid out of memory
    unit = unit(1:limit_spikes(2));
end

NumSpikes = length(unit);
vdisc = round((unit-starttime)*sr) + 1;
if NumSpikes < 5 % Minimum spike number for running STA
    disp(['Not enought spikes found in the current segment for ' cellid])
    return;
end

% STA CALCULATION, STAFIG
[st, sta, stats, H] = stacall_NBsync(vdisc,lfp,cellid, 'sr', sr, 'wn', wn,...
    'burstfilter', burstfilter, 'sync', sync, 'cellid2', cellid2, 'sessType',...
    sessType);   % call STA code

% Saving STA
cellidt = regexprep(cellid,'\.','_');
cellidt2 = regexprep(cellid2,'\.','_');
% ACG file for group identity
group = stats.group1{1};
groupID = stats.groupID1;
if ~isempty(cellid2)
    group2 = stats.group2{1};
    [groupID, Label] = groupCellPair(group, group2);
end

H2 = gcf;
ff1 = [resdir '\STAplots\' regexprep(sessType,' ','_') '\' sMode PATH '\'...
    group '_' cellidt '_' cellidt2 '_' num2str(index) '_' 'STA_new' numVer '.fig'];   % save STA figure
ff2 = [resdir '\STAjpegs\'  regexprep(sessType,' ','_') '\'  sMode PATH '\'...
    group '_' cellidt '_' cellidt2 '_' num2str(index) '_' 'STAimage_new' numVer '.jpg'];   % save STA jpeg
ff3 = ['D:\_MATLAB_DATA\NB\staNB\_allSTA\' cellidt '_' cellidt2 '_' numVer '.jpg'];   % save STA jpeg
saveas(H2,ff2)
saveas(H2,ff1)
saveas(H2,ff3)
close(H2);

% % Spectrogram
% try
%     [pow, phase, f] = eegwavelet(lfp,sr, 0.15); % EEG wavelet
%     [ersp, erspa] = ers(vdisc,pow,sr, wn); % Power spectrum
%     [ersph, erspha] = ers(vdisc,phase,sr, wn); % Phase spectrum
% catch % In case of running out of RAM (increase scale parameter)
%     try
%         keyboard;
%         [pow, phase, f] = eegwavelet(lfp,sr, 0.22);        % EEG wavelet
%         [ersp, erspa] = ers(vdisc,pow,sr, wn);
%     catch
%         keyboard;
%         [pow, phase, f] = eegwavelet(lfp,sr, 0.25);        % EEG wavelet
%         [ersp, erspa] = ers(vdisc,pow,sr, wn);
%     end
% end
% 
% % Convert to decibel
% I = log(ersp) - log(repmat(mean(ersp,2),1,size(ersp,2)));
% 
% % Store variables in struct
% stALL.(sessType).(burstfilter).allErsp = ersp;
% stALL.(sessType).(burstfilter).allPhase = ersph;
% stALL.(sessType).(burstfilter).erspLog = I;
% stALL.(sessType).(burstfilter).Grouping = groupID;
% stALL.(sessType).(burstfilter).scaleVector = f;
% stALL.(sessType).(burstfilter).allChAT = cellid;
% stALL.(sessType).(burstfilter).phases = phase;
% stALL.(sessType).(burstfilter).stats = stats;
stALL.(sessType).(burstfilter).st = st;
stALL.(sessType).(burstfilter).sta = sta;
SE = stats.se;
inverted = stats.inverted;

% [~,uppLim] = min(abs(f-100));
% fRange = size(I,1);
% fValues = f(uppLim:fRange); % cuts frequency range at 100Hz
% fValues = round(fValues, 1); % Frequency scale
% [~,uppGamma] = min(abs(fValues-40));
% [~,uppTheta] = min(abs(fValues-12));
% [~,uppDelta] = min(abs(fValues-4));
% fPos = [1 uppGamma uppTheta uppDelta length(fValues)];
% 
% I2 = I(uppLim:end,:);
% limit2 = size(I2,2)-1;
% H3 = figure;
% imagesc(I2)
% yticks(fPos);
% b_rescaleaxis('Y',round(fValues))
% setappdata(gca,'scaley',round(fValues))
% b_zoomset_for_wavelet
% xlim([0 limit2]);
% xticks([0 limit2/2 limit2]);
% xticklabels({num2str(-limit2/2), '0', num2str(limit2/2)});
% title([regexprep(cellid,'_',' ') ' , group: ' group])
% colorbar;
% colormap jet;
% cValues = caxis;
% cLimit{index} = cValues;
% ff1 = [resdir '\STAspect\' regexprep(sessType,' ','_') '\' sMode PATH '\' group '_' cellidt...
%     '_' cellidt2 '_' num2str(index) '_' 'STAspect_new' numVer '.fig'];   % save spectrum
% ff2 = [resdir '\STAjpegs\' regexprep(sessType,' ','_') '\' sMode PATH '\' group '_' cellidt...
%     '_' cellidt2 '_' num2str(index) '_' 'STAspect_new' numVer '.jpg'];   % save spectrum jpg
% saveas(H3,ff1)
% saveas(H3,ff2)
% close(H3);
% 
% % Plot phase
% ersph2 = ersph(uppLim:end,:);
% H4 = figure;
% imagesc(ersph2)
% yticks(fPos);
% b_rescaleaxis('Y',round(fValues))
% setappdata(gca,'scaley',round(fValues))
% b_zoomset_for_wavelet
% xlim([0 limit2]);
% xticks([0 limit2/2 limit2]);
% xticklabels({num2str(-limit2/2), '0', num2str(limit2/2)});
% title([regexprep(cellid,'_',' ') ' , group: ' group])
% colorbar;
% colormap jet;
% cValues = caxis;
% cLimit{index} = cValues;
% ff3 = [resdir '\STAspect\' regexprep(sessType,' ','_') '\' sMode PATH '\' group '_' cellidt...
%     '_' cellidt2 '_' num2str(index) '_' 'STAphase_new' numVer '.fig'];   % save spectrum
% ff4 = [resdir '\STAjpegs\' regexprep(sessType,' ','_') '\' sMode PATH '\' group '_' cellidt...
%     '_' cellidt2 '_' num2str(index) '_' 'STAphase_new' numVer '.jpg'];   % save spectrum jpg
% saveas(H4,ff3)
% saveas(H4,ff4)
% close(H4);
if strcmp(numVer, '1')
    if index==1
        invertList(1) = inverted;
        cellids{1} = cellid;
        save('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat', 'invertList', 'cellids');
    else
        load('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat');
        invertList(index) = inverted;
        cellids{index} = cellid;
        save('D:\_MATLAB_DATA\NB\staNB\_ALLDATA\invertList.mat', 'invertList', 'cellids');
    end
end
staFile2 = ['D:\_MATLAB_DATA\NB\staNB\_allSTAmatrix\' cellidt '_' numVer '.mat'];
staFile = [resdir '\SpectData\' regexprep(sessType,' ','_') '\' sMode PATH...
    '\' cellidt 'staMatrix_new' '12' '.mat'];
% phaseFile = [resdir '\SpectData\' regexprep(sessType,' ','_') '\' sMode PATH...
%     '\' cellidt 'phaseMatrix_new' numVer '.mat'];
% spectFile = [resdir '\SpectData\' regexprep(sessType,' ','_') '\' sMode PATH...
%     '\' cellidt 'spectMatrix_new' numVer '.mat'];
save(staFile, 'sta','SE','cellid','groupID','inverted', '-v7.3');
save(staFile2, 'sta','SE','cellid','groupID','inverted', '-v7.3');
% save(phaseFile, 'ersph','cellid','groupID','f', '-v7.3');
% save(spectFile, 'I','cellid','groupID', 'f', '-v7.3');

% -------------------------------------------------------------------------
function [sessType, burstfilter, sync] = cellChecker(sessType, burstfilter, sync)

if ~iscell(sessType)
    sessType={sessType};
end
if ~iscell(burstfilter)
    burstfilter={burstfilter};
end
if ~iscell(sync)
    sync={sync};
end