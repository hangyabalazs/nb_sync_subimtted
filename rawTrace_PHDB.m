function rawTrace_PHDB(cellid)



dbstop if error;
global RESDIR;
cellbase = whichcb;
cbSwitcher(cellid);
[r, s, tetrodename] = cellid2tags(cellid);   %#ok<*ASGLU> % get tetrode number
matname = cellid2fnames(cellid,'cont',tetrodename);   % .mat filename for LFP
[pathname, filename, extension] = fileparts(matname);   % parse filename

filename = ['CSC' num2str(tetrodename)];
cscname = fullfile(pathname,[filename '.ncs']);   % filename for the Neuralynx CSC file

% Load spikes
cellSpikes = loadcb(cellid);

% Detect bursts
[bursts, singleSpikes, burstTimes] = burstDetect(cellSpikes);
burstTimes = burstTimes{1};
burstTimes = burstTimes *1e5;
cellSpikes = cellSpikes *1e5;

[tdata, ts, info] = load_open_ephys_data([ pathname '\' '100_CH' num2str(19)  '.continuous']);

figure;
hold on;
plot(tdata)
oneStars = ones(numel(burstTimes),1);
plot(burstTimes, oneStars, '*', 'Color', [1 0 0])
plot(cellSpikes, ones(numel(cellSpikes*1e6),1), '*', 'Color', [0 1 0])