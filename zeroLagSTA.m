function zeroLagSTA
%EVENTSTA   STA aligned to StimulusOn for Hit, FA, CR, and Miss
%   It loads datamatrix from groupStaErs
%
%   See also STA_ERS, GROUPSTAERS
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

choosecb('NB');
cb = whichcb;
fs = filesep;
global RESDIR;
global cCode;
%
% NUMVER8 WAS USED FOR THE MANUSCRIPT
%


numVer = '8';
resdir = [RESDIR cb fs 'sta' cb fs '_ALLDATA'];
maxMode = false;
load([resdir fs 'masterData_new' numVer '.mat']);
verType = {'stimHIT', 'stimFA', 'CR', 'Miss'};
groupID = master.(verType{1}).none.original.groups;
phasicB = groupID == 1;
poissonL = groupID == 2;
tonic = groupID == 3;

% STA Inverted to be positive on bar plots
hitSTA = master.(verType{1}).none.original.sta*-1; % HIT
faSTA = master.(verType{2}).none.original.sta*-1; % FA
crSTA = master.(verType{3}).none.original.sta*-1; % CR
missSTA = master.(verType{4}).none.original.sta*-1; % Miss

% Plot parameters
wn = size(hitSTA,2)-1; % Plotwindow
time = linspace(-wn/2,wn/2,length(hitSTA)); % Time vector
curLag = round(size(time,2)/2); % Zero lag for current sta

% Use min positions insted of 0 lag
if maxMode
    % Keep only min values
    hitSTA = max(abs(hitSTA), [], 2);
    faSTA = max(abs(faSTA), [], 2);
    crSTA = max(abs(crSTA), [], 2);
    missSTA = max(abs(missSTA), [], 2);
    curLag = 1;
end

% STA value at curlag
hitB = hitSTA(phasicB,curLag);
hitP = hitSTA(poissonL,curLag);
hitT = hitSTA(tonic,curLag);
faB = faSTA(phasicB,curLag);
faP = faSTA(poissonL,curLag);
faT = faSTA(tonic,curLag);
crB = crSTA(phasicB,curLag);
crP = crSTA(poissonL,curLag);
crT = crSTA(tonic,curLag);
missB = missSTA(phasicB,curLag);
missP = missSTA(poissonL,curLag);
missT = missSTA(tonic,curLag);

% SE
SE_hitB = std(hitB) / sqrt(numel(hitB));
SE_hitP = std(hitP) / sqrt(numel(hitP));
SE_hitT = std(hitT) / sqrt(numel(hitT));
SE_faB = std(faB) / sqrt(numel(faB));
SE_faP = std(faP) / sqrt(numel(faP));
SE_faT = std(faT) / sqrt(numel(faT));
SE_crB = std(crB) / sqrt(numel(crB));
SE_crP = std(crP) / sqrt(numel(crP));
SE_crT = std(crT) / sqrt(numel(crT));
SE_missB = std(missB) / sqrt(numel(missB));
SE_missP = std(missP) / sqrt(numel(missP));
SE_missT = std(missT) / sqrt(numel(missT));





%% PhasicB

xmin=-0.25;
xmax=0.25;
% x coordinate generation to scatter the plot
x1 = xmin+rand(numel(hitB),1)*(xmax-xmin)+1;
x2 = xmin+rand(numel(faB),1)*(xmax-xmin)+2;
x3 = xmin+rand(numel(crB),1)*(xmax-xmin)+3;
x4 = xmin+rand(numel(missB),1)*(xmax-xmin)+4;

% Boxplots
[H1, pB(1)] = boxstat(hitB,faB,'Hit','FA',0.05,'paired');
close(H1);
[H1, pB(2)] = boxstat(crB,missB,'CR','Miss',0.05,'paired');
close(H1);
[H1, pB(3)] = boxstat(hitB,crB,'Hit','CR',0.05,'paired');
close(H1);
[H1, pB(4)] = boxstat(hitB,missB,'Hit','Miss',0.05,'paired');
close(H1);
[H1, pB(5)] = boxstat(faB,crB,'FA','CR',0.05,'paired');
close(H1);
[H1, pB(6)] = boxstat(faB,missB,'FA','Miss',0.05,'paired');
close(H1);

% Significance
p(1) = round(pB(1),3);
p(2) = round(pB(2),3);
p(3) = round(pB(3),3);
p(4) = round(pB(4),3);
p(5) = round(pB(5),3);
p(6) = round(pB(6),3);

groupInc = [p(1)<=0.05, p(2)<=0.05, p(3)<0.05, p(4)<0.05, p(5)<0.05, p(6)<0.05];

groups={[1,2],[3,4],[1,3],[1,4],[2,3],[2,4]};
groups = groups(groupInc);
pAll = p(groupInc);


% Hit-FA-CR-Miss plot phasicB
tickPos = [1 2 3 4]; % Bar x positions
H3 = figure;
hold on;
plot(x1, hitB, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(3,:), 'LineWidth',3)
plot(x2, faB, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(1,:), 'LineWidth',3)
plot(x3, crB, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0.75], 'LineWidth',3)
plot(x4, missB, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth',3)
TG = bar(tickPos, [mean(hitB), mean(faB), mean(crB), mean(missB)], 'FaceColor','none', 'LineWidth',2);
TG.EdgeColor = 'flat';
TG.CData(1,:) = cCode(3,:);
TG.CData(2,:) = cCode(1,:);
TG.CData(3,:) = [0 0 0.6];
TG.CData(4,:) = [0 0 0];
ylabel('mV');
xlabel('Groups');
xticks(tickPos)
xticklabels({'Hit' 'FA', 'CR', ' Miss'});
xlim([0 5])
ylim([-2000 7400])
yticks([-2000 0 2000 4000 6000])
yticklabels({'-2', '0', '2', '4', '6'});
sigstar(groups,pAll);
setmyplot_tamas;
axis square;
saveas(H3, [resdir '\ExtraPlots\fourBars\PhasicB_Bar_new_' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\fourBars\PhasicB_Bar_new_' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\fourBars\PhasicB_Bar_new_' numVer '.jpeg']);
close(H3);

%% PoissonL
% x coordinate generation to scatter the plot
x1 = xmin+rand(numel(hitP),1)*(xmax-xmin)+1;
x2 = xmin+rand(numel(faP),1)*(xmax-xmin)+2;
x3 = xmin+rand(numel(crP),1)*(xmax-xmin)+3;
x4 = xmin+rand(numel(missP),1)*(xmax-xmin)+4;


% Boxplots
[H1, pP(1)] = boxstat(hitP,faP,'Hit','FA',0.05,'paired');
close(H1);
[H1, pP(2)] = boxstat(crP,missP,'CR','Miss',0.05,'paired');
close(H1);
[H1, pP(3)] = boxstat(hitP,crP,'Hit','CR',0.05,'paired');
close(H1);
[H1, pP(4)] = boxstat(hitP,missP,'Hit','Miss',0.05,'paired');
close(H1);
[H1, pP(5)] = boxstat(faP,crP,'FA','CR',0.05,'paired');
close(H1);
[H1, pP(6)] = boxstat(faP,missP,'FA','Miss',0.05,'paired');
close(H1);
% p1 = round(p1,3);

% Significance
p(1) = round(pP(1),3);
p(2) = round(pP(2),3);
p(3) = round(pP(3),3);
p(4) = round(pP(4),3);
p(5) = round(pP(5),3);
p(6) = round(pP(6),3);

groupInc = [p(1)<=0.05, p(2)<=0.05, p(3)<0.05, p(4)<0.05, p(5)<0.05, p(6)<0.05];

groups={[1,2],[3,4],[1,3],[1,4],[2,3],[2,4]};
groups = groups(groupInc);
pAll = p(groupInc);



% Hit-FA-CR-Miss plot poissonL
tickPos = [1 2 3 4]; % Bar x positions
H3 = figure;
hold on;
plot(x1, hitP, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(3,:), 'LineWidth',3)
plot(x2, faP, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(1,:), 'LineWidth',3)
plot(x3, crP, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0.75], 'LineWidth',3)
plot(x4, missP, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth',3)
TG = bar(tickPos, [mean(hitP), mean(faP), mean(crP), mean(missP)], 'FaceColor','none', 'LineWidth',2);
TG.EdgeColor = 'flat';
TG.CData(1,:) = cCode(3,:);
TG.CData(2,:) = cCode(1,:);
TG.CData(3,:) = [0 0 0.6];
TG.CData(4,:) = [0 0 0];
ylabel('mV');
xlabel('Groups');
xticks(tickPos)
xticklabels({'Hit' 'FA', 'CR', ' Miss'});
xlim([0 5])
ylim([-2000 7400])
yticks([-2000 0 2000 4000 6000])
yticklabels({'-2', '0', '2', '4', '6'});
sigstar(groups,pAll);
setmyplot_tamas;
axis square;
saveas(H3, [resdir '\ExtraPlots\fourBars\PoissonL_Bar_new_' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\fourBars\PoissonL_Bar_new_' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\fourBars\PoissonL_Bar_new_' numVer '.jpeg']);
close(H3);

%% Tonic

% x coordinate generation to scatter the plot
x1 = xmin+rand(numel(hitT),1)*(xmax-xmin)+1;
x2 = xmin+rand(numel(faT),1)*(xmax-xmin)+2;
x3 = xmin+rand(numel(crT),1)*(xmax-xmin)+3;
x4 = xmin+rand(numel(missT),1)*(xmax-xmin)+4;

% Boxplots
[H1, pT(1)] = boxstat(hitT,faT,'Hit','FA',0.05,'paired');
close(H1);
[H1, pT(2)] = boxstat(crT,missT,'CR','Miss',0.05,'paired');
close(H1);
[H1,  pT(3)] = boxstat(hitT,crT,'Hit','CR',0.05,'paired');
close(H1);
[H1,  pT(4)] = boxstat(hitT,missT,'Hit','Miss',0.05,'paired');
close(H1);
[H1,  pT(5)] = boxstat(faT,crT,'FA','CR',0.05,'paired');
close(H1);
[H1,  pT(6)] = boxstat(faT,missT,'FA','Miss',0.05,'paired');
close(H1);
% p1 = round(p1,3);

% Significance
p(1) = round(pT(1),3);
p(2) = round(pT(2),3);
p(3) = round(pT(3),3);
p(4) = round(pT(4),3);
p(5) = round(pT(5),3);
p(6) = round(pT(6),3);

groupInc = [p(1)<=0.05, p(2)<=0.05, p(3)<0.05, p(4)<0.05, p(5)<0.05, p(6)<0.05];

groups={[1,2],[3,4],[1,3],[1,4],[2,3],[2,4]};
groups = groups(groupInc);
pAll = p(groupInc);



% Hit-FA-CR-Miss plot poissonL
tickPos = [1 2 3 4]; % Bar x positions
H3 = figure;
hold on;
plot(x1, hitT, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(3,:), 'LineWidth',3)
plot(x2, faT, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', cCode(1,:), 'LineWidth',3)
plot(x3, crT, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0.75], 'LineWidth',3)
plot(x4, missT, 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'LineWidth',3)
TG = bar(tickPos, [mean(hitT), mean(faT), mean(crT), mean(missT)], 'FaceColor','none', 'LineWidth',2);
TG.EdgeColor = 'flat';
TG.CData(1,:) = cCode(3,:);
TG.CData(2,:) = cCode(1,:);
TG.CData(3,:) = [0 0 0.6];
TG.CData(4,:) = [0 0 0];
ylabel('mV');
xlabel('Groups');
xticks(tickPos)
xticklabels({'Hit' 'FA', 'CR', ' Miss'});
xlim([0 5])
ylim([-2000 7400])
yticks([-2000 0 2000 4000 6000])
yticklabels({'-2', '0', '2', '4', '6'});
sigstar(groups,pAll);
setmyplot_tamas;
axis square;
saveas(H3, [resdir '\ExtraPlots\fourBars\Tonic_Bar_new_' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\fourBars\Tonic_Bar_new_' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\fourBars\Tonic_Bar_new_' numVer '.jpeg']);
close(H3);


