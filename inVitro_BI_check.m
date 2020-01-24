function newSpikeTimes = inVitro_BI_check(cellid, spTimes)
% INVITRO_BI_CHECK
% Downsamples in vivo spikes and recalculates BurstIndex
% for the shorter time windows. Time windows are selected
% to match the in vitro firing rate, and recording length.
%
% Check also the corresponding plot function for the results:
% See also INVITRO_BI_RESULTS

% Check if data already exists
dataName = which('daniData.mat');

%% Data loading/structuring
if ~exist(dataName)
    % Load Dani's BI data
    fName = which('burst_index_Spikenum_rec_length_export.xlsx');
    [resdir,~,~] = fileparts(fName);
    [numData, txtData, ~] = xlsread(fName);
    
    % Define Bursting and Late firing groups
    BFfilt = contains(txtData(2:end,1), 'BF');
    LFfilt = ~BFfilt;
    
    % Code groupID: bursting - 1, late firing - 2
    gID = BFfilt+(LFfilt*2);
    
    % Store data in a .mat file
    daniData.bIndex = numData(:,1);
    daniData.spikeNum = numData(:,2);
    daniData.recLength = numData(:,3);
    daniData.spikeFreq = daniData.spikeNum./daniData.recLength;
    daniData.groupID = gID;
    
    % save the .mat file
    save(fullfile(resdir, 'daniData.mat'), 'daniData')
else
    load(dataName);
end

% Load inVivo acg data
ivName = which('ACG_matrices_POOLED.mat');
ivData = load(ivName);


%% Match inVivo and inVitro data

cellMatch = contains(ivData.cellids, cellid);
ivType = ivData.groupID(cellMatch);

invivoLength = spTimes(end)-spTimes(1);
iVoFreq = numel(spTimes)/invivoLength;

% Select in vitro filter based on in vivo groupID
if ivType<3 % in vitro type: BF
    filtR = daniData.groupID==1;
else % in vitro type LF
    filtR = daniData.groupID==2;
end

Freq = daniData.spikeFreq(filtR);
Spikes = daniData.spikeNum(filtR);
Length = daniData.recLength(filtR);

% Find in vitro cell with the closest freq to the in vivo cell
[~, cellPos] = min(abs(Freq-iVoFreq));

% select in vitro cells recording properties
vitroSpikes = Spikes(cellPos);
vitroLength = Length(cellPos);


if vitroSpikes<length(spTimes)
    % Define in vivo segments based on matching spike number
    segStarts = 1:vitroSpikes:length(spTimes);
    segEnds = segStarts(2:end)-1;
    segs = [segStarts(1:end-1); segEnds];
    
    % Find the closest in vivo segment matching in vitro length
    segLengths = spTimes(segs(2,:))-spTimes(segs(1,:));
    [minDiff, segPos] = min(abs(segLengths-vitroLength));
    
    % save segment length difference
    [dirName,~,~] = fileparts(dataName);
    save([dirName filesep regexprep(cellid,'\.','_') '_lengthDiff.mat'],...
        'minDiff');
    
    % Downsample
    newSpikeTimes = spTimes(segs(1,segPos):segs(2,segPos));
else
    newSpikeTimes = spTimes;
end
