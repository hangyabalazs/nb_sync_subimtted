function eventSTA
%EVENTSTA   STA aligned to StimulusOn for Hit, FA, CR, and Miss
%   It loads datamatrix from groupStaErs
%
%   See also STA_ERS, GROUPSTAERS
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu
cb = whichcb;
fs = filesep;
global RESDIR;
global cCode;
resdir = [RESDIR cb fs 'sta' cb fs '_ALLDATA'];
numVer = '8';
load([resdir fs 'masterData_new' numVer '.mat']);

groupID = master.stimHIT.none.original.groups;
phasicB = groupID == 1;
poissonL = groupID == 2;
tonic = groupID == 3;


% STA
hitSTA = master.stimHIT.none.original.sta; % HIT
faSTA = master.stimFA.none.original.sta; % FA
crSTA = master.CR.none.original.sta; % CR
missSTA = master.Miss.none.original.sta; % Miss

wn = size(hitSTA,2)-1;
time = linspace(-wn/2,wn/2,length(hitSTA));

% STA event plot
% Bursting
H1 = figure;
hold on;
errorshade(time,mean(hitSTA(phasicB,:), 1), std(hitSTA(phasicB,:), [], 1)/sqrt(sum(phasicB)),'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',1.5)
errorshade(time,mean(faSTA(phasicB,:), 1), std(faSTA(phasicB,:), [], 1)/sqrt(sum(phasicB)),'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',1.5)
errorshade(time,mean(crSTA(phasicB,:), 1), std(crSTA(phasicB,:), [], 1)/sqrt(sum(phasicB)),'LineColor',[0 0 1],'ShadeColor',[0 0 1],'LineWidth',1.5)
errorshade(time,mean(missSTA(phasicB,:), 1), std(missSTA(phasicB,:), [], 1)/sqrt(sum(phasicB)),'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',1.5)
xlim([-wn/2 wn/2]);
xticks([-wn/2 0 wn/2]);
xlabel('Time from tone (ms)');
y_lim = ylim;
y_lim = round(y_lim,1,'significant');
ylim(y_lim);
yticks([y_lim(1) 0 y_lim(2)]);
ylabel('uV');
title('   ');
axis square;
legend({'Hit', 'FA', 'CR', 'Miss'}, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
saveas(H1, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_phasicB_new' numVer '.fig']);
saveas(H1, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_phasicB_new' numVer '.pdf']);
saveas(H1, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_phasicB_new' numVer '.jpeg']);
close(H1);

% PoissonL
H2 = figure;
hold on;
errorshade(time,mean(hitSTA(poissonL,:), 1), std(hitSTA(poissonL,:), [], 1)/sqrt(sum(poissonL)),'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',1.5)
errorshade(time,mean(faSTA(poissonL,:), 1), std(faSTA(poissonL,:), [], 1)/sqrt(sum(poissonL)),'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',1.5)
errorshade(time,mean(crSTA(poissonL,:), 1), std(crSTA(poissonL,:), [], 1)/sqrt(sum(poissonL)),'LineColor',[0 0 1],'ShadeColor',[0 0 1],'LineWidth',1.5)
errorshade(time,mean(missSTA(poissonL,:), 1), std(missSTA(poissonL,:), [], 1)/sqrt(sum(poissonL)),'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',1.5)
xlim([-wn/2 wn/2]);
xticks([-wn/2 0 wn/2]);
xlabel('Time from tone (ms)');
ylim(y_lim);
yticks([y_lim(1) 0 y_lim(2)]);
ylabel('uV');
title('   ');
axis square;
legend({'Hit', 'FA', 'CR', 'Miss'}, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
saveas(H2, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_poissonL_new' numVer '.fig']);
saveas(H2, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_poissonL_new' numVer '.pdf']);
saveas(H2, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_poissonL_new' numVer '.jpeg']);
close(H2);

% Tonic
H3 = figure;
hold on;
errorshade(time,mean(hitSTA(tonic,:), 1), std(hitSTA(tonic,:), [], 1)/sqrt(sum(tonic)),'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',1.5)
errorshade(time,mean(faSTA(tonic,:), 1), std(faSTA(tonic,:), [], 1)/sqrt(sum(tonic)),'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',1.5)
errorshade(time,mean(crSTA(tonic,:), 1), std(crSTA(tonic,:), [], 1)/sqrt(sum(tonic)),'LineColor',[0 0 1],'ShadeColor',[0 0 1],'LineWidth',1.5)
errorshade(time,mean(missSTA(tonic,:), 1), std(missSTA(tonic,:), [], 1)/sqrt(sum(tonic)),'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',1.5)
xlim([-wn/2 wn/2]);
xticks([-wn/2 0 wn/2]);
xlabel('Time from tone (ms)');
ylim(y_lim);
yticks([y_lim(1) 0 y_lim(2)]);
ylabel('uV');
title('   ');
axis square;
legend({'Hit', 'FA', 'CR', 'Miss'}, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
saveas(H3, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_tonic_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_tonic_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\_new\Hit_vs_FA_cr_miss_tonic_new' numVer '.jpeg']);
close(H3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([resdir fs 'barData' numVer '.mat'])

% PhasicB THETA-GAMMA
hitB_T = barData.stimHIT.none.original.Bursting(2,:); % Raw
faB_T = barData.stimFA.none.original.Bursting(2,:); % Burst1
hitG_T = barData.stimHIT.none.original.Bursting(3,:); % Single
faG_T = barData.stimFA.none.original.Bursting(3,:); % Sync

% FourBar Theta-Gamma
eventT = [mean(hitB_T), mean(faB_T), mean(hitG_T),...
    mean(faG_T)];

% SE: Theta-Gamma for phasicB
SE_hitBT = std(hitB_T) / sqrt(numel(hitB_T));
SE_faBT = std(faB_T) / sqrt(numel(faB_T));
SE_hitGT = std(hitG_T) / sqrt(numel(hitG_T));
SE_faGT = std(faG_T) / sqrt(numel(faG_T));
SE_T = [SE_hitBT SE_faBT SE_hitGT SE_faGT];

% Boxplots
[H1, p1] = boxstat(hitB_T,faB_T,'Hit','FA',0.05,'paired');
close(H1);
[H2, p2] = boxstat(hitG_T,faG_T,'Hit','FA',0.05,'paired');
close(H2);
p1 = round(p1,3);
p2 = round(p2,3);

% Theta fourBar plot
tickPos = [1 2 3.5 4.5]; % Bar x positions
H3 = figure;
hold on;
TG = bar(tickPos, eventT, 'FaceColor','none', 'LineWidth',2);
E = errorbar(tickPos,eventT,SE_T, '.');
E.LineWidth = 2;
E.Color = [0.6 0.6 0.6];
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 1 0];
TG.CData(2,:) = [1 0 0];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
axis square;
ylabel('Decibel');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
% xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 5])
title('PhasicB Theta-Gamma')
% ylim([0 0.4])
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fourBar_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fourBar_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fourBar_new' numVer '.jpeg']);
close(H3);

keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PoissonL THETA-GAMMA
hitB_T = barData.stimHIT.none.original.Poisson(2,:); % Raw
faB_T = barData.stimFA.none.original.Poisson(2,:); % Burst1
hitG_T = barData.stimHIT.none.original.Poisson(3,:); % Single
faG_T = barData.stimFA.none.original.Poisson(3,:); % Sync

% FourBar Theta-Gamma
eventT = [mean(hitB_T), mean(faB_T), mean(hitG_T),...
    mean(faG_T)];

% SE: Theta-Gamma for phasicB
SE_hitBT = std(hitB_T) / sqrt(numel(hitB_T));
SE_faBT = std(faB_T) / sqrt(numel(faB_T));
SE_hitGT = std(hitG_T) / sqrt(numel(hitG_T));
SE_faGT = std(faG_T) / sqrt(numel(faG_T));
SE_T = [SE_hitBT SE_faBT SE_hitGT SE_faGT];

% Boxplots
[H1, p1] = boxstat(hitB_T,faB_T,'Hit','FA',0.05,'paired');
close(H1);
[H2, p2] = boxstat(hitG_T,faG_T,'Hit','FA',0.05,'paired');
close(H2);
p1 = round(p1,3);
p2 = round(p2,3);

% Theta fiveBar plot
tickPos = [1 2 3.5 4.5]; % Bar x positions
H3 = figure;
hold on;
TG = bar(tickPos, eventT, 'FaceColor','none', 'LineWidth',2);
E = errorbar(tickPos,eventT,SE_T, '.');
E.LineWidth = 2;
E.Color = [0.6 0.6 0.6];
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 1 0];
TG.CData(2,:) = [1 0 0];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
axis square;
ylabel('Decibel');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
% xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 5])
title('PoissonL Theta-Gamma')
% ylim([0 0.4])
saveas(H3, [resdir '\ExtraPlots\PoissonL_Theta_fourBar_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\PoissonL_Theta_fourBar_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\PoissonL_Theta_fourBar_new' numVer '.jpeg']);
close(H3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tonic THETA-GAMMA
hitB_T = barData.stimHIT.none.original.Tonic(2,:); % Raw
faB_T = barData.stimFA.none.original.Tonic(2,:); % Burst1
hitG_T = barData.stimHIT.none.original.Tonic(3,:); % Single
faG_T = barData.stimFA.none.original.Tonic(3,:); % Sync

% FourBar Theta-Gamma
eventT = [mean(hitB_T), mean(faB_T), mean(hitG_T),...
    mean(faG_T)];

% SE: Theta-Gamma for phasicB
SE_hitBT = std(hitB_T) / sqrt(numel(hitB_T));
SE_faBT = std(faB_T) / sqrt(numel(faB_T));
SE_hitGT = std(hitG_T) / sqrt(numel(hitG_T));
SE_faGT = std(faG_T) / sqrt(numel(faG_T));
SE_T = [SE_hitBT SE_faBT SE_hitGT SE_faGT];

% Boxplots
[H1, p1] = boxstat(hitB_T,faB_T,'Hit','FA',0.05,'paired');
close(H1);
[H2, p2] = boxstat(hitG_T,faG_T,'Hit','FA',0.05,'paired');
close(H2);
p1 = round(p1,3);
p2 = round(p2,3);

% Theta fiveBar plot
tickPos = [1 2 3.5 4.5]; % Bar x positions
H3 = figure;
hold on;
TG = bar(tickPos, eventT, 'FaceColor','none', 'LineWidth',2);
E = errorbar(tickPos,eventT,SE_T, '.');
E.LineWidth = 2;
E.Color = [0.6 0.6 0.6];
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 1 0];
TG.CData(2,:) = [1 0 0];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
axis square;
ylabel('Decibel');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
% xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 5])
title('Tonic Theta-Gamma')
% ylim([0 0.4])
saveas(H3, [resdir '\ExtraPlots\Tonic_Theta_fourBar_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\Tonic_Theta_fourBar_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\Tonic_Theta_fourBar_new' numVer '.jpeg']);
close(H3);
keyboard;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CR - MISS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STA event plot
% Bursting
H1 = figure;
hold on;
errorshade(time,mean(crSTA(phasicB,:), 1), std(crSTA(phasicB,:), [], 1)/sqrt(sum(phasicB)),'LineColor',[0 0 1],'ShadeColor',[0 0 1],'LineWidth',1.5)
errorshade(time,mean(missSTA(phasicB,:), 1), std(missSTA(phasicB,:), [], 1)/sqrt(sum(phasicB)),'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',1.5)
xlim([-wn/2 wn/2]);
xticks([-wn/2 0 wn/2]);
xlabel('Lag (ms)');
y_lim = ylim;
y_lim = round(y_lim,1,'significant');
ylim(y_lim);
yticks([y_lim(1) 0 y_lim(2)]);
ylabel('uV');
title('   ');
axis square;
legend({'CR', 'Miss'}, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
saveas(H1, [resdir '\ExtraPlots\cr_vs_miss_phasicB_new' numVer '.fig']);
saveas(H1, [resdir '\ExtraPlots\cr_vs_miss_phasicB_new' numVer '.pdf']);
saveas(H1, [resdir '\ExtraPlots\cr_vs_miss_phasicB_new' numVer '.jpeg']);
close(H1);

% PoissonL
H2 = figure;
hold on;
errorshade(time,mean(crSTA(poissonL,:), 1), std(crSTA(poissonL,:), [], 1)/sqrt(sum(poissonL)),'LineColor',[0 0 1],'ShadeColor',[0 0 1],'LineWidth',1.5)
errorshade(time,mean(missSTA(poissonL,:), 1), std(missSTA(poissonL,:), [], 1)/sqrt(sum(poissonL)),'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',1.5)
xlim([-wn/2 wn/2]);
xticks([-wn/2 0 wn/2]);
xlabel('Lag (ms)');
ylim(y_lim);
yticks([y_lim(1) 0 y_lim(2)]);
ylabel('uV');
title('   ');
axis square;
legend({'CR', 'Miss'}, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
saveas(H2, [resdir '\ExtraPlots\cr_vs_miss_poissonL_new' numVer '.fig']);
saveas(H2, [resdir '\ExtraPlots\cr_vs_miss_poissonL_new' numVer '.pdf']);
saveas(H2, [resdir '\ExtraPlots\cr_vs_miss_poissonL_new' numVer '.jpeg']);
close(H2);

% Tonic
H3 = figure;
hold on;
errorshade(time,mean(crSTA(tonic,:), 1), std(crSTA(tonic,:), [], 1)/sqrt(sum(tonic)),'LineColor',[0 0 1],'ShadeColor',[0 0 1],'LineWidth',1.5)
errorshade(time,mean(missSTA(tonic,:), 1), std(missSTA(tonic,:), [], 1)/sqrt(sum(tonic)),'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',1.5)
xlim([-wn/2 wn/2]);
xticks([-wn/2 0 wn/2]);
xlabel('Lag (ms)');
ylim(y_lim);
yticks([y_lim(1) 0 y_lim(2)]);
ylabel('uV');
title('   ');
axis square;
legend({'CR', 'Miss'}, 'Location','northeast');
legend('boxoff');
setmyplot_tamas;
saveas(H3, [resdir '\ExtraPlots\cr_vs_miss_tonic_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\cr_vs_miss_tonic_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\cr_vs_miss_tonic_new' numVer '.jpeg']);
close(H3);


% PhasicB THETA-GAMMA
hitB_T = barData.CR.none.original.Bursting(2,:); % Raw
faB_T = barData.Miss.none.original.Bursting(2,:); % Burst1
hitG_T = barData.CR.none.original.Bursting(3,:); % Single
faG_T = barData.Miss.none.original.Bursting(3,:); % Sync

% FourBar Theta-Gamma
eventT = [mean(hitB_T), mean(faB_T), mean(hitG_T),...
    mean(faG_T)];

% SE: Theta-Gamma for phasicB
SE_hitBT = std(hitB_T) / sqrt(numel(hitB_T));
SE_faBT = std(faB_T) / sqrt(numel(faB_T));
SE_hitGT = std(hitG_T) / sqrt(numel(hitG_T));
SE_faGT = std(faG_T) / sqrt(numel(faG_T));
SE_T = [SE_hitBT SE_faBT SE_hitGT SE_faGT];

% Boxplots
[H1, p1] = boxstat(hitB_T,faB_T,'CR','Miss',0.05,'paired');
close(H1);
[H2, p2] = boxstat(hitG_T,faG_T,'CR','Miss',0.05,'paired');
close(H2);
p1 = round(p1,3);
p2 = round(p2,3);

% Theta fourBar plot
tickPos = [1 2 3.5 4.5]; % Bar x positions
H3 = figure;
hold on;
TG = bar(tickPos, eventT, 'FaceColor','none', 'LineWidth',2);
E = errorbar(tickPos,eventT,SE_T, '.');
E.LineWidth = 2;
E.Color = [0.6 0.6 0.6];
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 1 0];
TG.CData(2,:) = [1 0 0];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
axis square;
ylabel('Decibel');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
% xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 5])
title('PhasicB Theta-Gamma')
% ylim([0 0.4])
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fourBar_CRMiss_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fourBar_CRMiss_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\PhasicB_Theta_fourBar_CRMiss_new' numVer '.jpeg']);
close(H3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PoissonL THETA-GAMMA
hitB_T = barData.CR.none.original.Poisson(2,:); % Raw
faB_T = barData.Miss.none.original.Poisson(2,:); % Burst1
hitG_T = barData.CR.none.original.Poisson(3,:); % Single
faG_T = barData.Miss.none.original.Poisson(3,:); % Sync

% FourBar Theta-Gamma
eventT = [mean(hitB_T), mean(faB_T), mean(hitG_T),...
    mean(faG_T)];

% SE: Theta-Gamma for phasicB
SE_hitBT = std(hitB_T) / sqrt(numel(hitB_T));
SE_faBT = std(faB_T) / sqrt(numel(faB_T));
SE_hitGT = std(hitG_T) / sqrt(numel(hitG_T));
SE_faGT = std(faG_T) / sqrt(numel(faG_T));
SE_T = [SE_hitBT SE_faBT SE_hitGT SE_faGT];

% Boxplots
[H1, p1] = boxstat(hitB_T,faB_T,'CR','Miss',0.05,'paired');
close(H1);
[H2, p2] = boxstat(hitG_T,faG_T,'CR','Miss',0.05,'paired');
close(H2);
p1 = round(p1,3);
p2 = round(p2,3);

% Theta fiveBar plot
tickPos = [1 2 3.5 4.5]; % Bar x positions
H3 = figure;
hold on;
TG = bar(tickPos, eventT, 'FaceColor','none', 'LineWidth',2);
E = errorbar(tickPos,eventT,SE_T, '.');
E.LineWidth = 2;
E.Color = [0.6 0.6 0.6];
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 1 0];
TG.CData(2,:) = [1 0 0];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
axis square;
ylabel('Decibel');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
% xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 5])
title('PoissonL Theta-Gamma')
% ylim([0 0.4])
saveas(H3, [resdir '\ExtraPlots\PoissonL_Theta_fourBar_CRMiss_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\PoissonL_Theta_fourBar_CRMiss_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\PoissonL_Theta_fourBar_CRMiss_new' numVer '.jpeg']);
close(H3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tonic THETA-GAMMA
hitB_T = barData.CR.none.original.Tonic(2,:); % Raw
faB_T = barData.Miss.none.original.Tonic(2,:); % Burst1
hitG_T = barData.CR.none.original.Tonic(3,:); % Single
faG_T = barData.Miss.none.original.Tonic(3,:); % Sync

% FourBar Theta-Gamma
eventT = [mean(hitB_T), mean(faB_T), mean(hitG_T),...
    mean(faG_T)];

% SE: Theta-Gamma for phasicB
SE_hitBT = std(hitB_T) / sqrt(numel(hitB_T));
SE_faBT = std(faB_T) / sqrt(numel(faB_T));
SE_hitGT = std(hitG_T) / sqrt(numel(hitG_T));
SE_faGT = std(faG_T) / sqrt(numel(faG_T));
SE_T = [SE_hitBT SE_faBT SE_hitGT SE_faGT];

% Boxplots
[H1, p1] = boxstat(hitB_T,faB_T,'CR','Miss',0.05,'paired');
close(H1);
[H2, p2] = boxstat(hitG_T,faG_T,'CR','Miss',0.05,'paired');
close(H2);
p1 = round(p1,3);
p2 = round(p2,3);

% Theta fiveBar plot
tickPos = [1 2 3.5 4.5]; % Bar x positions
H3 = figure;
hold on;
TG = bar(tickPos, eventT, 'FaceColor','none', 'LineWidth',2);
E = errorbar(tickPos,eventT,SE_T, '.');
E.LineWidth = 2;
E.Color = [0.6 0.6 0.6];
TG.EdgeColor = 'flat';
TG.CData(1,:) = [0 1 0];
TG.CData(2,:) = [1 0 0];
TG.CData(3,:) = cCode(3,:);
TG.CData(4,:) = cCode(1,:);
axis square;
ylabel('Decibel');
xlabel('Groups');
setmyplot_tamas;
xticks(tickPos)
% xticklabels({'AllSpike' 'Burst1', 'Single', 'Sync', 'Async'});
xlim([0.5 5])
title('Tonic Theta-Gamma')
% ylim([0 0.4])
saveas(H3, [resdir '\ExtraPlots\Tonic_Theta_fourBar_CRMiss_new' numVer '.fig']);
saveas(H3, [resdir '\ExtraPlots\Tonic_Theta_fourBar_CRMiss_new' numVer '.pdf']);
saveas(H3, [resdir '\ExtraPlots\Tonic_Theta_fourBar_CRMiss_new' numVer '.jpeg']);
close(H3);

