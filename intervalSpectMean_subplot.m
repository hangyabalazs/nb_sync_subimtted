function intervalSpectMean_subplot



cb = whichcb;
fs = filesep;
global RESDIR;
global cCode;
%
% NUMVER8 WAS USED FOR THE MANUSCRIPT
%
fileDir = 'D:\_MATLAB_DATA\NB\intervalSpectMean';

numVer = '15';
resdir = [RESDIR cb fs 'sta' cb fs '_ALLDATA'];
maxMode = true;
load([resdir fs 'masterData_new' numVer '.mat']);
% load([resdir fs 'invertList' numVer '.mat']);
verType = {'CueHit', 'CueFA', 'CRhit', 'Missfa'};
% verType = {'stimHIT', 'stimFA', 'CR', 'Miss'};
groupID = master.(verType{1}).none.original.groups;
phasicB = groupID == 1;
poissonL = groupID == 2;
tonic = groupID == 3;

% STA Inverted to be positive on bar plots
hitSpect = master.(verType{1}).none.original.spect; % HIT
faSpect = master.(verType{2}).none.original.spect; % FA
crSpect = master.(verType{3}).none.original.spect; % CR
missSpect = master.(verType{4}).none.original.spect; % Miss

% Plot parameters
wn = size(hitSpect{1},2)-1; % Plotwindow
time = linspace(-wn/2,wn/2,length(hitSpect{1})); % Time vector
curLag = round(size(time,2)/2); % Zero lag for current sta

spectGroups = {hitSpect, faSpect, crSpect, missSpect};

deltainx = 33:52;
thetainx = 21:32;
gammainx = 10:20;
for iV = 1:length(spectGroups)
    currVer = spectGroups{iV};
    for aP = 1:length(currVer)
        currHz = currVer{aP}; % spect data
        centerP = round(size(currHz, 2)/2); % center for plotting
        deltaS.(verType{iV})(aP,:) = mean(currHz(deltainx,:)); % delta range <4Hz
        thetaS.(verType{iV})(aP,:) = mean(currHz(thetainx,:)); % theta range ~4Hz-12Hz
        gammaS.(verType{iV})(aP,:) = mean(currHz(gammainx,:)); % gamma range ~12Hz-41Hz
    end
end

%% Figure By groups for all trialtypes

% PhasicB MEAN
dHit = mean(deltaS.CueHit(phasicB,:),1);
tHit = mean(thetaS.CueHit(phasicB,:),1);
gHit = mean(gammaS.CueHit(phasicB,:),1);

dFA = mean(deltaS.CueFA(phasicB,:),1);
tFA = mean(thetaS.CueFA(phasicB,:),1);
gFA = mean(gammaS.CueFA(phasicB,:),1);

dCR = mean(deltaS.CRhit(phasicB,:),1);
tCR = mean(thetaS.CRhit(phasicB,:),1);
gCR = mean(gammaS.CRhit(phasicB,:),1);

dMiss = mean(deltaS.Missfa(phasicB,:),1);
tMiss = mean(thetaS.Missfa(phasicB,:),1);
gMiss = mean(gammaS.Missfa(phasicB,:),1);


% PhasicB SE
dHitS = std(deltaS.CueHit(phasicB,:), [], 1)/sqrt(sum(phasicB));
tHitS = std(thetaS.CueHit(phasicB,:), [], 1)/sqrt(sum(phasicB));
gHitS = std(gammaS.CueHit(phasicB,:), [], 1)/sqrt(sum(phasicB));

dFAS = std(deltaS.CueFA(phasicB,:), [], 1)/sqrt(sum(phasicB));
tFAS = std(thetaS.CueFA(phasicB,:), [], 1)/sqrt(sum(phasicB));
gFAS = std(gammaS.CueFA(phasicB,:), [], 1)/sqrt(sum(phasicB));

dCRS = std(deltaS.CRhit(phasicB,:), [], 1)/sqrt(sum(phasicB));
tCRS = std(thetaS.CRhit(phasicB,:), [], 1)/sqrt(sum(phasicB));
gCRS = std(gammaS.CRhit(phasicB,:), [], 1)/sqrt(sum(phasicB));

dMissS = std(deltaS.Missfa(phasicB,:), [], 1)/sqrt(sum(phasicB));
tMissS = std(thetaS.Missfa(phasicB,:), [], 1)/sqrt(sum(phasicB));
gMissS = std(gammaS.Missfa(phasicB,:), [], 1)/sqrt(sum(phasicB));

LWidth = 2;

% H1 = figure;
% 
% 
% % subplot(3,3,1);
% hold on;
% errorshade(time,dHit, dHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
% errorshade(time,dFA,dFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
% errorshade(time,dCR,dCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
% errorshade(time,dMiss,dMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
% title(['BFCN-burstSB Delta (<4Hz)'])
% % legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% % legend('boxoff');
% % setmyplot_tamas;
% ylim([-1 1])
% % xlim(plotTime);
% % xticks(x_ticks);
% % xticklabels(xtick_labels);
% xlabel(['Time from event (ms)']);
% axis square;
% % fName = [fileDir 'phasicB_Delta.fig'];   % save PSTH figure
% % fNameJ = [fileDir 'phasicB_Delta.jpeg'];   % save PSTH figure
% % saveas(H1,fName);
% % saveas(H1,fNameJ);

H1 = figure;
% subplot(3,3,2);
hold on;
errorshade(time,tHit, tHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
errorshade(time,tFA,tFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
errorshade(time,tCR,tCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
errorshade(time,tMiss,tMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
title(['BFCN-burstSB Theta (4-12 Hz)'])
% legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% legend('boxoff');
setmyplot_tamas;
ylim([-1.5 1])
% xlim(plotTime);
% xticks(x_ticks);
% xticklabels(xtick_labels);
xlabel(['Time from event (ms)']);
axis square;

fName = [fileDir fs 'allPlot_phasicB_theta.fig'];   % save PSTH figure
fNameJ = [fileDir fs 'allPlot_phasicB_theta.jpeg'];   % save PSTH figure
saveas(H1,fName);
saveas(H1,fNameJ);
% fName = [fileDir 'phasicB_Theta.fig'];   % save PSTH figure
% fNameJ = [fileDir 'phasicB_Theta.jpeg'];   % save PSTH figure
% saveas(H12,fName);
% saveas(H12,fNameJ);

% H1 = figure;
% subplot(3,3,3);
% hold on;
% errorshade(time,gHit, gHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
% errorshade(time,gFA,gFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
% errorshade(time,gCR,gCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
% errorshade(time,gMiss,gMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
% title(['BFCN-burstSB Gamma (12-41 Hz)'])
% % legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% % legend('boxoff');
% % setmyplot_tamas;
% ylim([-1 1])
% % xlim(plotTime);
% % xticks(x_ticks);
% % xticklabels(xtick_labels);
% xlabel(['Time from event (ms)']);
% axis square;
% % fName = [fileDir 'phasicB_Gamma.fig'];   % save PSTH figure
% % fNameJ = [fileDir 'phasicB_Gamma.jpeg'];   % save PSTH figure
% % saveas(H13,fName);
% % saveas(H13,fNameJ);


%% PoissonL
dHit = mean(deltaS.CueHit(poissonL,:),1);
tHit = mean(thetaS.CueHit(poissonL,:),1);
gHit = mean(gammaS.CueHit(poissonL,:),1);

dFA = mean(deltaS.CueFA(poissonL,:),1);
tFA = mean(thetaS.CueFA(poissonL,:),1);
gFA = mean(gammaS.CueFA(poissonL,:),1);

dCR = mean(deltaS.CRhit(poissonL,:),1);
tCR = mean(thetaS.CRhit(poissonL,:),1);
gCR = mean(gammaS.CRhit(poissonL,:),1);

dMiss = mean(deltaS.Missfa(poissonL,:),1);
tMiss = mean(thetaS.Missfa(poissonL,:),1);
gMiss = mean(gammaS.Missfa(poissonL,:),1);
LWidth = 2;


% PoissonL SE
dHitS = std(deltaS.CueHit(poissonL,:), [], 1)/sqrt(sum(poissonL));
tHitS = std(thetaS.CueHit(poissonL,:), [], 1)/sqrt(sum(poissonL));
gHitS = std(gammaS.CueHit(poissonL,:), [], 1)/sqrt(sum(poissonL));

dFAS = std(deltaS.CueFA(poissonL,:), [], 1)/sqrt(sum(poissonL));
tFAS = std(thetaS.CueFA(poissonL,:), [], 1)/sqrt(sum(poissonL));
gFAS = std(gammaS.CueFA(poissonL,:), [], 1)/sqrt(sum(poissonL));

dCRS = std(deltaS.CRhit(poissonL,:), [], 1)/sqrt(sum(poissonL));
tCRS = std(thetaS.CRhit(poissonL,:), [], 1)/sqrt(sum(poissonL));
gCRS = std(gammaS.CRhit(poissonL,:), [], 1)/sqrt(sum(poissonL));

dMissS = std(deltaS.Missfa(poissonL,:), [], 1)/sqrt(sum(poissonL));
tMissS = std(thetaS.Missfa(poissonL,:), [], 1)/sqrt(sum(poissonL));
gMissS = std(gammaS.Missfa(poissonL,:), [], 1)/sqrt(sum(poissonL));



% subplot(3,3,4);
% hold on;
% errorshade(time,dHit, dHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
% errorshade(time,dFA,dFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
% errorshade(time,dCR,dCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
% errorshade(time,dMiss,dMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
% title(['BFCN-burstPL Delta (<4Hz)'])
% % legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% % legend('boxoff');
% % setmyplot_tamas;
% ylim([-1 1])
% % xlim(plotTime);
% % xticks(x_ticks);
% % xticklabels(xtick_labels);
% xlabel(['Time from event (ms)']);
% axis square;
% % fName = [fileDir 'poissonL_Delta.fig'];   % save PSTH figure
% % fNameJ = [fileDir 'poissonL_Delta.jpeg'];   % save PSTH figure
% % saveas(H2,fName);
% % saveas(H2,fNameJ);

H2 = figure;
% subplot(3,3,5);
hold on;
errorshade(time,tHit, tHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
errorshade(time,tFA,tFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
errorshade(time,tCR,tCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
errorshade(time,tMiss,tMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
title(['BFCN-burstPL Theta (4-12 Hz)'])
% legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% legend('boxoff');
setmyplot_tamas;
ylim([-1.5 1])
% xlim(plotTime);
% xticks(x_ticks);
% xticklabels(xtick_labels);
xlabel(['Time from event (ms)']);
axis square;
fName = [fileDir fs 'allPlot_poissonL_theta.fig'];   % save PSTH figure
fNameJ = [fileDir fs 'allPlot_poissonL_theta.jpeg'];   % save PSTH figure
saveas(H2,fName);
saveas(H2,fNameJ);

% fName = [fileDir 'poissonL_Theta.fig'];   % save PSTH figure
% fNameJ = [fileDir 'poissonL_Theta.jpeg'];   % save PSTH figure
% saveas(H22,fName);
% saveas(H22,fNameJ);

% subplot(3,3,6);
% hold on;
% errorshade(time,gHit, gHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
% errorshade(time,gFA,gFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
% errorshade(time,gCR,gCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
% errorshade(time,gMiss,gMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
% 
% plot(time,gHit,'Color', cCode(3,:), 'LineWidth', LWidth)
% plot(time,gFA,'Color', cCode(1,:), 'LineWidth', LWidth)
% plot(time,gCR,'Color', [0 0 0.75], 'LineWidth', LWidth)
% plot(time,gMiss,'Color', [0 0 0], 'LineWidth', LWidth)
% title(['BFCN-burstPL Gamma (12-41 Hz)'])
% % legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% % legend('boxoff');
% % setmyplot_tamas;
% ylim([-1 1])
% % xlim(plotTime);
% % xticks(x_ticks);
% % xticklabels(xtick_labels);
% xlabel(['Time from event (ms)']);
% axis square;
% % fName = [fileDir 'poissonL_Gamma.fig'];   % save PSTH figure
% % fNameJ = [fileDir 'poissonL_Gamma.jpeg'];   % save PSTH figure
% % saveas(H23,fName);
% % saveas(H23,fNameJ);


%% Tonic
dHit = mean(deltaS.CueHit(tonic,:),1);
tHit = mean(thetaS.CueHit(tonic,:),1);
gHit = mean(gammaS.CueHit(tonic,:),1);

dFA = mean(deltaS.CueFA(tonic,:),1);
tFA = mean(thetaS.CueFA(tonic,:),1);
gFA = mean(gammaS.CueFA(tonic,:),1);

dCR = mean(deltaS.CRhit(tonic,:),1);
tCR = mean(thetaS.CRhit(tonic,:),1);
gCR = mean(gammaS.CRhit(tonic,:),1);

dMiss = mean(deltaS.Missfa(tonic,:),1);
tMiss = mean(thetaS.Missfa(tonic,:),1);
gMiss = mean(gammaS.Missfa(tonic,:),1);
LWidth = 2;

% Tonic SE
dHitS = std(deltaS.CueHit(tonic,:), [], 1)/sqrt(sum(tonic));
tHitS = std(thetaS.CueHit(tonic,:), [], 1)/sqrt(sum(tonic));
gHitS = std(gammaS.CueHit(tonic,:), [], 1)/sqrt(sum(tonic));

dFAS = std(deltaS.CueFA(tonic,:), [], 1)/sqrt(sum(tonic));
tFAS = std(thetaS.CueFA(tonic,:), [], 1)/sqrt(sum(tonic));
gFAS = std(gammaS.CueFA(tonic,:), [], 1)/sqrt(sum(tonic));

dCRS = std(deltaS.CRhit(tonic,:), [], 1)/sqrt(sum(tonic));
tCRS = std(thetaS.CRhit(tonic,:), [], 1)/sqrt(sum(tonic));
gCRS = std(gammaS.CRhit(tonic,:), [], 1)/sqrt(sum(tonic));

dMissS = std(deltaS.Missfa(tonic,:), [], 1)/sqrt(sum(tonic));
tMissS = std(thetaS.Missfa(tonic,:), [], 1)/sqrt(sum(tonic));
gMissS = std(gammaS.Missfa(tonic,:), [], 1)/sqrt(sum(tonic));




% subplot(3,3,7);
% hold on;
% errorshade(time,dHit, dHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
% errorshade(time,dFA,dFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
% errorshade(time,dCR,dCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
% errorshade(time,dMiss,dMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
% title(['BFCN-REG Delta (<4Hz)'])
% % legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% % legend('boxoff');
% % setmyplot_tamas;
% ylim([-1 1])
% % xlim(plotTime);
% % xticks(x_ticks);
% % xticklabels(xtick_labels);
% xlabel(['Time from event (ms)']);
% axis square;
% % fName = [fileDir 'tonic_Delta.fig'];   % save PSTH figure
% % fNameJ = [fileDir 'tonic_Delta.jpeg'];   % save PSTH figure
% % saveas(H3,fName);
% % saveas(H3,fNameJ);

H3 = figure;
% subplot(3,3,8);
hold on;
errorshade(time,tHit, tHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
errorshade(time,tFA,tFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
errorshade(time,tCR,tCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
errorshade(time,tMiss,tMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
title(['BFCN-REG Theta (4-12 Hz)'])
% legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% legend('boxoff');
setmyplot_tamas;
ylim([-1.5 1])
% xlim(plotTime);
% xticks(x_ticks);
% xticklabels(xtick_labels);
xlabel(['Time from event (ms)']);
axis square;
fName = [fileDir fs 'allPlot_tonic_theta.fig'];   % save PSTH figure
fNameJ = [fileDir fs 'allPlot_tonic_theta.jpeg'];   % save PSTH figure
saveas(H3,fName);
saveas(H3,fNameJ);

% fName = [fileDir 'tonic_Theta.fig'];   % save PSTH figure
% fNameJ = [fileDir 'tonic_Theta.jpeg'];   % save PSTH figure
% saveas(H32,fName);
% saveas(H32,fNameJ);

% subplot(3,3,9);
% hold on;
% errorshade(time,gHit, gHitS,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',LWidth)
% errorshade(time,gFA,gFAS,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',LWidth)
% errorshade(time,gCR,gCRS,'LineColor',[0 0 0.75],'ShadeColor',[0 0 0.75],'LineWidth',LWidth)
% errorshade(time,gMiss,gMissS,'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',LWidth)
% title(['BFCN-REG Gamma (12-41 Hz)'])
% % legend({'HIT', 'FA', 'CR','MISS'}, 'FontSize', 12, 'Location','southeast')
% % legend('boxoff');
% % setmyplot_tamas;
% ylim([-1 1])
% % xlim(plotTime);
% % xticks(x_ticks);
% % xticklabels(xtick_labels);
% xlabel(['Time from event (ms)']);
% axis square;
% % fName = [fileDir 'tonic_Gamma.fig'];   % save PSTH figure
% % fNameJ = [fileDir 'tonic_Gamma.jpeg'];   % save PSTH figure
% % saveas(H33,fName);
% % saveas(H33,fNameJ);

% maximize_figure;
% fName = [fileDir fs 'allPlot_allEvent.fig'];   % save PSTH figure
% fNameJ = [fileDir fs 'allPlot_allEvent.jpeg'];   % save PSTH figure
% saveas(H1,fName);
% saveas(H1,fNameJ);
