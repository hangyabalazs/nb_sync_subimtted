function peakLatencySTA
%PEAKLATENCYSTA   Plots STA min distance from zero lag.
%   See also STA_ERS, GROUPSTAERS

dbstop if error;
global RESDIR;
global cCode;
choosecb('NB');
cb = whichcb;
fs = filesep;
resdir = [RESDIR cb '\sta' cb '\_ALLDATA'];
fnm = [resdir fs 'masterData_new3.mat'];
load(fnm);

% verType = {'stimHIT', 'stimFA', 'CR', 'Miss'};
groupID = master.behavior.none.original.groups;
phasicB = groupID == 1;
poissonL = groupID == 2;
tonic = groupID == 3;

% for iV = 1:length(verType)
% Store all sta
allsta = master.behavior.none.original.sta;
centerP = round(size(allsta,2)/2); % Zero lag
numCell = size(allsta,1);
minpos = nan(1,numCell);
for iC = 1:numCell
    currCell = allsta(iC,:);
    minpos(iC) = find(currCell==min(currCell));
end

% Find STA min distance from Zero lag
latency = minpos - centerP;
lat_mean = round(mean(latency),2);
lat_med = round(median(latency),2);
se_lat = round(std(latency) / sqrt(numel(latency)));
se_of_median_lat = se_of_median(latency);


% Significance test
% Boxplots
[H1, p1] = boxstat(latency(phasicB),latency(poissonL),'Burst-SB','Burst-PL',0.05,'nonpaired');
close(H1);
[H2, p2] = boxstat(latency(phasicB),latency(tonic),'Burst-SB','Regular',0.05,'nonpaired');
close(H2);
[H3, p3] = boxstat(latency(poissonL),latency(tonic),'Burst-PL','Regular',0.05,'nonpaired');
close(H3);

p(1) = round(p1,3);
p(2) = round(p2,3);
p(3) = round(p3,3);

groupInc = [p1<=0.05, p2<=0.05, p3<0.05];

groups={[1,2],[1,3],[2,3]};
groups = groups(groupInc);
pAll = p(groupInc);

% Barplot for Peak latencies
H4 = figure;
hold on;
bar(1,nanmedian(latency(phasicB)), 'FaceColor','none','EdgeColor',cCode(1,:), 'LineWidth',3)
bar(2,nanmedian(latency(poissonL)), 'FaceColor','none','EdgeColor',cCode(2,:), 'LineWidth',3)
bar(3,nanmedian(latency(tonic)), 'FaceColor','none','EdgeColor',cCode(3,:), 'LineWidth',3)


plot(1, latency(phasicB), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7], 'LineWidth',1)
plot(2, latency(poissonL), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7], 'LineWidth',1)
plot(3, latency(tonic), 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', [0.7 0.7 0.7], 'LineWidth',1)

sigstar(groups,pAll);


ylabel('STA peak latency (ms)');
plotText = {['median: ' num2str(lat_med) ' ms + SEM: ' num2str(se_of_median_lat)],...
    ['mean: ' num2str(lat_mean) ' ms + SEM: ' num2str(se_lat)]};
axis square;
%     xlim([0 2]);
x_lim = xlim;
y_lim = ylim;
yticks([y_lim(1), 0, y_lim(2)/2 y_lim(2)])
%     text((x_lim(2)*0.35), y_lim(2)*0.85, plotText, 'Color', 'k');
setmyplot_tamas;
saveas(H4, 'D:\_MATLAB_DATA\NB\STApeakLatency\latencyPlot.fig');
saveas(H4, 'D:\_MATLAB_DATA\NB\STApeakLatency\latencyPlot.jpeg');
close(H4);

% % Hist for latencies
% figure;
% hist(latency)
% xlabel('STA Peak latency (ms)');
% setmyplot_tamas;

% end