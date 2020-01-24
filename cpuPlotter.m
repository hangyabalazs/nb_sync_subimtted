function cpuPlotter

global RESDIR;
global PATH;
fs = filesep;
cb = 'HDB';
dataDir = [RESDIR cb fs 'PSTH_newdata' cb fs 'CPu' fs];
listDir = dir([dataDir 'EventData*']);
resdir = [RESDIR cb fs 'psth' cb fs 'figures' fs 'CPu' fs];
% Load all psth the data for every cellbase
% Find latest EventData file
edFiles = {listDir.name};

if ~isempty(edFiles)
    remInx = find(cellfun('isempty',strfind(edFiles,'removed')));
    currFiles = {edFiles{remInx}};
    dates = [listDir(remInx).datenum];
    latestD = find(dates==max(dates));
    allData = load([dataDir currFiles{latestD}]);  
else
    allData = loadSeparate(outLiers);
end

EventData = allData.EventData.FB;

groups = EventData.FA.NotAdaptive.none.stats;

% Define TP groups
for i = 1:length(groups)
    switch groups{i}.group
        case 1
            groupID(i) = 1;
        case 2
            groupID(i) = 2;
        case 3
            groupID(i) = 3;
    end
end

phasicB = groupID==1;
poissonL = groupID==2;
tonic = groupID==3;

datSize = EventData.FA.NotAdaptive.none.spsth;
lag = size(datSize,2);
middle = round(size(datSize,2)/2);
actWin = [-250 750];
plotWin = [middle+actWin(1):middle+actWin(2)];
tType = {'FA', 'HIT'};


for iT = 1:numel(tType)
    currData = EventData.(tType{iT}).NotAdaptive.none.spsth;
    tonicDat2 = currData(tonic,:);
    tonicDat = [tonicDat2(3,:); tonicDat2(2,:); tonicDat2(1,:);...
        tonicDat2(4,:); tonicDat2(5,:)];
    H1 = figure;
    imagesc(tonicDat)
    colormap jet
    axis square;
    setmyplot_tamas;
    xticks([751 1001 1251 1751])
%     xticklabels({num2str(actWin(1)), '0', num2str(actWin(2)/2),...
%         num2str(actWin(2))})
    xticklabels({num2str(actWin(1)), '0', '250',...
        num2str(actWin(2))})

    if strcmp('FA',tType{iT})
        xlabel('Time from punishment (ms)')
    else
        xlabel('Time from reward (ms)')
    end
    ylabel('Cells')
    title([(tType{iT})])
    xlim([middle+actWin(1) 1751])
    fnmS4 = [resdir 'striatal_TAN_' tType{iT} 'win250_750.fig'];
    fnmSE4 = [resdir 'striatal_TAN_' tType{iT} 'win250_750.jpeg'];
    saveas(H1, fnmS4);
    saveas(H1, fnmSE4);
    close(H1);
end











