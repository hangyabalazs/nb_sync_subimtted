function burstRatio(cb)
%BURSTRATIO  Pie charts comparing the number of burst spikes vs single
% spikes for all cholinergic cells.
%
%   See also BURSTDETECT.

%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Initialize
dbstop if error;
global RESDIR;
fs = filesep;
% Load ACG for current cellbase
acgData = load([RESDIR cb fs 'acg' cb fs 'ACG_matrices_' cb '.mat']);

% Load data from acg data
allChAT = acgData.cellids;
NumChAT = length(allChAT);
BurstIndex = acgData.BurstIndex;
Refractory = acgData.Refractory;

% TP group indexes
phasicB = Refractory < 40 & BurstIndex > 0.2;   % phasic, bursting
poissonL = Refractory < 40 & BurstIndex <= 0.2;   % poisson-like
tonic = Refractory >= 40;   % refractory above 40 ms

% Preallocate
[spikeNum, burst1Num, burstAllNum, allNum] = deal(zeros(NumChAT, 1));
for i = 1:NumChAT
    cellid = allChAT(i);
    cbSwitcher(cellid);
    spikeTimes = loadcb(cellid,'SPIKES');
    
    % Burst Detection
    [burst1Spikes, singleSpikes, burstAllSpikes] = burstDetect(spikeTimes);
    spikeNum(i) = size(cell2mat(singleSpikes),1); % n of single spikes
    burst1Num(i) = size(cell2mat(burst1Spikes),1); % n of burst1 spikes
    burstAllNum(i) = size(cell2mat(burstAllSpikes),1); % n of burstall spikes
    allNum(i) = size(spikeTimes,1); % n of all spikes
end


%% Pie chart of Burstall vs Allspikes
% Subplots
%Bursting
H1 = figure;
ax1=subplot(1,3,1);
plotDat1 = [sum(burstAllNum(phasicB)), sum(spikeNum(phasicB))];
pie(ax1, plotDat1)
title(ax1,'Bursting');

% Poisson-like
ax2 = subplot(1,3,2);
plotDat2 = [sum(burstAllNum(poissonL)), sum(spikeNum(poissonL))];
pie(ax2, plotDat2)
title(ax2,'Poisson');

% Tonic
if sum(tonic)>0
ax3 = subplot(1,3,3);
plotDat3 = [sum(burstAllNum(tonic)), sum(spikeNum(tonic))];
pie(ax3, plotDat3)
title(ax3,'Tonic');
saveas(H1,[RESDIR cb '\burstRatio' cb '\pieChart_3plot.fig'])   % save plot
saveas(H1,[RESDIR cb '\burstRatio' cb '\pieChart_3plot.jpeg'])   % save plot
end
close(H1);
%% Plot them individually
% Separate plots for burst ratio

% Bursting
H2 = figure;
plotDat1 = [sum(burstAllNum(phasicB)), sum(spikeNum(phasicB))];
pie(plotDat1)
title('Bursting');
axis square;
setmyplot_tamas;
saveas(H2,[RESDIR cb '\burstRatio' cb '\pieChart_burst.fig'])   % save plot
saveas(H2,[RESDIR cb '\burstRatio' cb '\pieChart_burst.jpeg'])   % save plot
close(H2);

% Poisson-like
H3 = figure;
plotDat2 = [sum(burstAllNum(poissonL)), sum(spikeNum(poissonL))];
pie(plotDat2)
title('Poisson');
axis square;
setmyplot_tamas;
saveas(H3,[RESDIR cb '\burstRatio' cb '\pieChart_poissonL.fig'])   % save plot
saveas(H3,[RESDIR cb '\burstRatio' cb '\pieChart_poissonL.jpeg'])   % save plot
close(H3);

% Tonic
if sum(tonic)>0
H4 = figure;
plotDat3 = [sum(burstAllNum(tonic)), sum(spikeNum(tonic))];
pie(plotDat3)
title('Tonic');
axis square;
setmyplot_tamas;
saveas(H4,[RESDIR cb '\burstRatio' cb '\pieChart_tonic.fig'])   % save plot
saveas(H4,[RESDIR cb '\burstRatio' cb '\pieChart_tonic.jpeg'])   % save plot
close(H4);
end
