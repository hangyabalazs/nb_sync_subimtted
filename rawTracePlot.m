function rawTracePlot(cellid, mode, meanDiff)
%RAWTRACEPLOT   Selects a typical 1s long  episode from a continous LFP
%   recording.
%   RAWTRACEPLOT(CELLID, MODE, MEANDIFF) Raw trace plot from continous
%   Neuralynx or OE file.
%
%   Input arguments:
%       CELLID - cellID of example cell
%       MODE - 1: Bursting, 2: Tonic 3:Poisson
%       MEANDIFF - if true - subtract noise from another tetrode
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Initialize variables
dbstop if error;
global RESDIR;
exploreM = false; % Controlling exploratory plotting
fs = filesep;

% Select appropriate cellbase based on cellid
cbSwitcher(cellid);
cb = whichcb;

% Prepare filenames
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
try
    if ~strcmp('PannaHDB', cb)
        [TimeStamp, ChanNum, SampleFrequency, NumValSamples, Samples,...
            NlxHeader] = Nlx2MatCSC(cscname,[1 1 1 1 1],1,1,1);  %  LFP
        data = Samples(1:end)';   % convert Neuralynx data to column vector
    else
        dataDir = [RESDIR 'POOLED' fs 'rawTracePOOLED'];
        read_openephys('datadir', pathname, 'resdir',dataDir, 'TTspec',...
            tetrodename);
    end
catch
    disp(['No CSC file for the current cell: ' cellid])
    return;
end

if meanDiff
    noise = 'subtracted_';
else
    noise = '';
end

% Only for NLX data
if mode < 3
    %% Loading continous data from another tetrode to subtract common noise
    if meanDiff
        switch mode
            case 1
                tetrodename2 = 8;
            case 2
                tetrodename2 = 4;
            case 3
                tetrodename2 = 8;
        end
        filename2 = ['CSC' num2str(tetrodename2)];
        cscname2 = fullfile(pathname,[filename2 '.ncs']);   % filename for the Neuralynx CSC file
        
        
        [TimeStamp2, ChanNum2, SampleFrequency2, NumValSamples2, Samples2, NlxHeader2] = ...
            Nlx2MatCSC(cscname2,[1 1 1 1 1],1,1,1);  %  LFP
        
        
        data2 = Samples2(1:end)';   % convert Neuralynx data to column vector
        data = data-data2;
    end
    
    TimeStamps_CSC = TimeStamp / 1e6;   % convert time to seconds
    NHs = NlxHeader{cellfun(@(s)~isempty(s),strfind(NlxHeader,'ADBitVolts'))};   % header scale variable
    NHsv = regexp(NHs,'\d\.\d*','match');
    ADBitVolts = str2double(NHsv{1});
    data = data * ADBitVolts * 1e06; % Convert data to Volts
    sr_orig = 512 / (mean(diff(TimeStamps_CSC)));   % sampling rate
    dt_orig = 1 /sr_orig;  % sampling time
    time_orig0 = TimeStamps_CSC(1):dt_orig:TimeStamps_CSC(1)+dt_orig*(length(data)-1);  % time vector (even)
    starttime = TimeStamps_CSC(1);   % first time stamp
    
    numData = length(time_orig0);
    pointer = 1:numData;
    Divisions = pointer(rem(numData,pointer)==0);
    
    % Exploratory mode, plot entire trace to select appropriate segment
    if exploreM
        % Plot the entire trace
        H0 = figure;
        hold on;
        plot(time_orig0, data); % raw data
        plot(burstTimes, ones(length(burstTimes), 1)*50, '*', 'Color', [0.8 0 0]) % burst positions
        plot(cellSpikes, ones(length(cellSpikes), 1)*1, '*', 'Color', [0 0.8 0]) % allspike positions
        close(H0);
        return;
    end
    
end

% Final plots
switch mode
    case 1 % Bursting cell
        windows = numData/Divisions(end-9);
        wins = [1 linspace(Divisions(end-9), numData, windows)];
        nW = 12; % Segment number to plot
        
        % Saving of Fig1. Panel B1, B2
        H1 = figure;
        hold on;
        plot(time_orig0(wins(nW):wins(nW+1)),...
            data(wins(nW):wins(nW+1)), 'Color', [0 0 0]);
        x_lim = [10435.2 10436.2]; % 1s window with bursts
        xlim(x_lim);
        ylim([-300 300])
        scalebar;
        setmyplot_tamas;
        axis off;
        %1s plot
        savefig(H1,[RESDIR cb '\rawtrace' cb '\' noise...
            'bursting_1s_example.fig'],'compact')
        %150 ms plot
        xlim([10435.45 10435.6]);
        scalebar;
        savefig(H1,[RESDIR cb '\rawtrace' cb '\' noise...
            'bursting_150ms_example.fig'],'compact')
        close(H1);
        
    case 2 % Tonic cell
        windows = numData/Divisions(end-12);
        wins = [1 linspace(Divisions(end-12), numData, windows)];
        nW=15; % Segment number to plot
        
        % Saving of Fig1. Panel C1
        H2 = figure;
        hold on;
        plot(time_orig0(wins(nW):wins(nW+1)),...
            data(wins(nW):wins(nW+1)), 'Color', [0 0 0]);
        x_lim = [2622.5 2623.5]; % 1s window with bursts
        xlim(x_lim);
        ylim([-300 300])
        scalebar;
        setmyplot_tamas;
        axis off;
        %1s plot
        savefig(H2,[RESDIR cb '\rawtrace' cb '\' noise...
            'tonic_1s_example.fig'],'compact')
        close(H2);
        
    case 3 % Poisson cell
        % Load data saved by read_openephys
        rawData = load([RESDIR 'POOLED' fs 'rawTracePOOLED' fs 'TT'...
            num2str(tetrodename) '_' s '_rawdata.mat']);
        chanNum = str2double(cellid{1}(end)); % Channel number
        currData = rawData.data(:,chanNum)*-1; % Invert data
        rawts = rawData.ts{chanNum}; % timestamps
        % Load Timestamps
        currTS = load(['G:' fs 'pavlovian_cholinergic_cellbase' fs...
            'HDB36' fs s '\TT' num2str(tetrodename) '_' cellid{1}(end)...
            '.mat']);
        sr = 1e4;
        ts10 = currTS.TS/sr;
        
        % Detect bursts
        [bursts, singleSpikes, burstTimes] = burstDetect(ts10);
        burstTimes = burstTimes{1};
        singleSpikes = singleSpikes{1};
        
        % Define plotWindow
        [minN, minI] = min(diff(singleSpikes));
        minPos = singleSpikes(minI);
        x_lim = [1763.5 1764.5];
        plotWin(1) = find(rawts==x_lim(1));
        plotWin(2) = find(rawts==x_lim(2));
        plotTS = plotWin(1):plotWin(2);
        
        % Plot figure
        H3 = figure;
        hold on;
        plot(rawts(plotTS), currData(plotTS), 'Color', [0 0 0]);
        x_lim = xlim;
        ylim([-300 300]);
        scalebar;
        setmyplot_tamas;
        axis off;
        %1s plot
        savefig(H3,[RESDIR 'NB\rawtraceNB\' noise...
            'poisson_1s_example.fig'],'compact')
        close(H3);   
end


