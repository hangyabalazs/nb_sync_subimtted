function rawTrace_Pplotter(cellid)

dbstop if error;
global RESDIR;
fs = filesep;
cellbase = whichcb;
cbSwitcher(cellid);
[r, s, tetrodename] = cellid2tags(cellid);   %#ok<*ASGLU> % get tetrode number
matname = cellid2fnames(cellid,'cont',tetrodename);   % .mat filename for LFP
[pathname, filename, extension] = fileparts(matname);   % parse filename

filename = ['CSC' num2str(tetrodename)];

chanNum = str2num(cellid(end));

rawData = load([RESDIR 'POOLED' fs 'rawTracePOOLED' fs 'TT' num2str(tetrodename) '_'...
    s '_rawdata.mat']);

currData = rawData.data(:,chanNum)*-1;
rawts = rawData.ts{chanNum};
currTS = load(['G:\pavlovian_cholinergic_cellbase\HDB36\' s '\TT' num2str(tetrodename)...
    '_' cellid(end) '.mat']);

ts10 = currTS.TS/10000;

% Detect bursts
[bursts, singleSpikes, burstTimes] = burstDetect(ts10);
burstTimes = burstTimes{1};
singleSpikes = singleSpikes{1};

figure;
hold on;
plot(rawts, currData);
plot(burstTimes, ones(1, numel(burstTimes)), '*', 'Color', [0.8 0 0]);
plot(singleSpikes, ones(1, numel(singleSpikes))+10, '*', 'Color', [0 0.8 0]);

