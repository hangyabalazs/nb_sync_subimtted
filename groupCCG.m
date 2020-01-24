function groupCCG(varargin)
%GROUPCCG   Cross correlation averages by TP groups.
%   GROUPCCG(BURSTFILTER, CCGPATH, ACGPATH, SESSIONTYPE) Calculates and
%   plots average cross-correlations for TP groups.
%
%   See also CCG, ACG and CELLID2TPGROUP.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Input arguments
prs = inputParser;
addParameter(prs,'ccgPath',[],@(s)isempty(s)|ischar(s))   % CCG data path
addParameter(prs,'acgPath',[],@(s)isempty(s)|ischar(s))   % ACG data path
addParameter(prs,'sessiontype','behavior',@(s)ischar(s)&ismember(s,{'behavior',...
    'ITI', 'sleep', 'quiet wakefulness', 'freely moving', 'selfCCG'}))   % Plot CCG and corresponding ACG's together
addParameter(prs,'burstfilter','none',@(s)ischar(s)&ismember(s,{'none',...
    'burst1', 'burstall', 'single' 'burst1Single'}))   % Plot CCG and corresponding ACG's together
addParameter(prs,'pairTypes','normal',@(s)ischar(s)&ismember(s,{'normal',...
    'putative', 'untagged', 'tagged', 'untaggedMix'}))
addParameter(prs,'cellbase',[],@(s)ischar(s)&ismember(s,{'NB',...
    'HDB', 'PannaHDB', 'POOLED'}))   % select cellbase
parse(prs,varargin{:})
g = prs.Results;

dbstop if error;
if isempty(g.cellbase)
    cb = 'NB';
else
    cb = g.cellbase;
end
global RESDIR;
global cCode;
global PATH;

% Normalization
mode = {'intragroup'};  % use intragroup for Z-score

% Load CCG matrix
ccgPath = g.ccgPath;
if ~isempty(ccgPath)
    load(ccgPath);   % Load CCG matrix
else
    ccgPath = [RESDIR cb '\ccg' cb '\' regexprep(g.sessiontype,' ','_')...
        '\' g.burstfilter cb PATH '\'];
    fnm_ccg = [ccgPath 'CCG_matrices_' regexprep(g.sessiontype,' ','_') '_'...
        g.burstfilter '_' cb '.mat'];
    try
        load(fnm_ccg);   % Load CCG matrix
    catch
        [ccg_file, ccg_filedir] = uigetfile([RESDIR cb '\ccg' cb '\'...
            regexprep(g.sessiontype,' ','_') '\' g.burstfilter cb PATH '\'],...
            'Choose CCG file!');
        if ccg_file < 1
            error('No ccg file has been selected. Program aborts.')
        end
        fnm_ccg = fullfile(ccg_filedir, ccg_file);
        load(fnm_ccg);   % Load CCG matrix
    end
end

% Finds which CCG matrix is selected (look for burstfilter)
if isempty(g.burstfilter)
    uPos = strfind(ccg_file, '_');
    bfilt = [ccg_file(uPos(2)+1:uPos(3)-1) ' ' cb ' spikes'];
    burstfilter = ccg_file(uPos(2)+1:uPos(3)-1);
    resdir1 = [RESDIR cb '\groupCCG_' cb '\' regexprep(g.sessiontype,' ','_') '\' burstfilter cb '\'];
else
    bfilt = [g.burstfilter ' ' cb ' spikes'];
    burstfilter = g.burstfilter;
    resdir1 = [RESDIR cb '\groupCCG' cb '\' regexprep(g.sessiontype,' ','_') '\' burstfilter cb '\'];
end

% Categorize the pairs based on burstiness
[groupID, yLabel, pairPath, chInx] = categorizeTP(PairOfCells, group1, group2, g.pairTypes);
groupID = groupID';
for iT = 1:length(groupID)
    ttInx(iT) = strcmp(PairOfCells{iT,1}(end-2), PairOfCells{iT,2}(end-2));
end


%%
% Load untagged ccg

if ~strcmp(whichcb, 'NB')
    choosecb('NB');
end

UT = load([RESDIR 'NB\ccgNB\behavior\noneNB\other\CCG_matrices_behavior_none_NB.mat']);
chatType1 = cellfun(@(x) getvalue('ChAT+', x), UT.PairOfCells(:,1));
chatType2 = cellfun(@(x) getvalue('ChAT+', x), UT.PairOfCells(:,2));

pairType = chatType1 | chatType2;
untagged = find(pairType == 0);
numPair_UT = numel(untagged);
CCR_UT = UT.CCR(untagged,:);
PairOfCells_UT = UT.PairOfCells(untagged,:);
for iC = 1:size(PairOfCells_UT,1)
   tetrodeP(iC) = strcmp(PairOfCells_UT{iC,1}(end-2), PairOfCells_UT{iC,2}(end-2)); 
end

numAll = numel(tetrodeP);
numTP = sum(tetrodeP);
numNTP = sum(~tetrodeP);

% DELETING AND AVERAGING CENTRAL DATABINS
plotWin = [round(((size(CCR_UT,2)-1)/2)*-1) round(((size(CCR_UT,2)-1)/2))]; % +/- 250ms lag
lags = [plotWin(1):1:-1, 1:1:plotWin(2)]; % +/- 1 s window, 1 ms sampling rate
centerP = round(((size(CCR_UT,2)-1)/2))+1; % Center point of the current data
CCR_UT(:,centerP) = mean([CCR_UT(:,centerP-1), CCR_UT(:,centerP+2)], 2); % Average the first and third databin
CCR_UT(:,centerP+1) = []; % deletes 2nd databin

CCRnorm_UT = [];
for k = 1:numPair_UT
    sccr = smooth(CCR_UT(k,:),'linear',15);   % smooth
    b = mean(sccr);
    normValue(k) = b;
    stdValue(k) = std(sccr);
    CCRnorm_UT(k,:) = (sccr - normValue(k)) / stdValue(k);  % mean-normalize
end

% All untagged
mn_all = median(CCRnorm_UT, 1);   % average CCG, tetrodepairs
mn_all = normalizem(mn_all,centerP);
se_all = se_of_median(CCRnorm_UT,1);   % SE, tetrodepairs
se_all = normalizem(se_all,centerP);


% Tedrodepairs
mn_all_TP = median(CCRnorm_UT(tetrodeP,:), 1);   % average CCG, tetrodepairs
mn_all_TP = normalizem(mn_all_TP,centerP);
se_all_TP = se_of_median(CCRnorm_UT(tetrodeP,:),1);   % SE, tetrodepairs
se_all_TP = normalizem(se_all_TP,centerP);
% Non-tetrodepairs
mn_all_NTP = median(CCRnorm_UT(~tetrodeP,:), 1);   % average CCG, non-tetrodepairs
mn_all_NTP = normalizem(mn_all_NTP,centerP);
se_all_NTP = se_of_median(CCRnorm_UT(~tetrodeP,:),1);   % SE, non-tetrodepairs
se_all_NTP = normalizem(se_all_NTP,centerP);

%% Synchrony index

actRatioTP = mean(CCR_UT(tetrodeP, (centerP-30):(centerP+30)),2);
baseRatioTP = mean(CCR_UT(tetrodeP, (centerP+100):end),2);
actRatioNTP = mean(CCR_UT(~tetrodeP, (centerP-30):(centerP+30)),2);
baseRatioNTP = mean(CCR_UT(~tetrodeP, (centerP+100):end),2);

ratioTP = (actRatioTP-baseRatioTP)./(actRatioTP+baseRatioTP);
ratioNTP = (actRatioNTP-baseRatioNTP)./(actRatioNTP+baseRatioNTP);


actRatioUT = mean(CCR_UT(:, (centerP-30):(centerP+30)),2);
baseRatioUT = mean(CCR_UT(:, (centerP+100):end),2);
ratioUT = (actRatioUT-baseRatioUT)./(actRatioUT+baseRatioUT);



%%

tetTypes = {'all', 'nontetrodepairs', 'tetrodepairs'};

for iP = 1:length(tetTypes)
    
    switch tetTypes{iP}
        case 'all'
            currTT = true(1,length(groupID));
        case 'nontetrodepairs'
            currTT = ~ttInx;
        case 'tetrodepairs'
            currTT = ttInx;
    end
    
    % Tetrode vs Non-tetrode pairs
    CCR = CCR(currTT,:);
    groupID = groupID(currTT);
    
%     % Subselect CCR
%     CCR = CCR(chInx,:);
    

    numPair = length(groupID);
    pairB = groupID == 1; % Bursting
    pairPL = groupID == 2; % Poisson-like
    pairT = groupID == 3; % Tonic
    pairEF = groupID == 1 | groupID == 2 | groupID == 4;
    pairB_PL = groupID == 4; % Bursting + Poisson-like
    pairB_T = groupID == 5; % Bursting + Tonic
    pairPL_T = groupID == 6; % Poisson-like + Tonic
    pairPcont = groupID == 2 | groupID == 4; % Poisson-like + Poisson-like-Burtsing mixed
    pairTcont = groupID == 3 | groupID == 5 | groupID == 6;
    
    if iP==1
        % DELETING AND AVERAGING CENTRAL DATABINS
        plotWin = [round(((size(CCR,2)-1)/2)*-1) round(((size(CCR,2)-1)/2))]; % +/- 1s lag
        lags = [plotWin(1):1:-1, 1:1:plotWin(2)]; % +/- 1 s window, 1 ms sampling rate
        centerP = round(((size(CCR,2)-1)/2))+1; % Center point of the current data
        CCR(:,centerP) = mean([CCR(:,centerP-1), CCR(:,centerP+2)], 2); % Average the first and third databin
        CCR(:,centerP+1) = []; % deletes 2nd databin
    end
    CCRnorm = [];
    if strcmp(burstfilter, 'none') % Intragroup: Normalize only for the current burstfilter
        for k = 1:numPair
            sccr = smooth(CCR(k,:),'linear',15);   % smooth
            b = mean(sccr);
            normValue(k) = b;
            stdValue(k) = std(sccr);
            CCRnorm(k,:) = (sccr - normValue(k)) / stdValue(k);  % mean-normalize
        end
        save([resdir1 PATH pairPath '\' 'normValue_' tetTypes{iP} '.mat'], 'normValue');
    else
        %     load([RESDIR cellbase '\groupCCG' cellbase '\' regexprep(g.sessiontype,' ','_') '\none' cellbase PATH '\normValue.mat']);
        for k = 1:numPair
            sccr = smooth(CCR(k,:),'linear',15);   % smooth
            b = mean(sccr);
            normValue(k) = b;
            stdValue(k) = std(sccr);
            CCRnorm(k,:) = (sccr - normValue(k)) / stdValue(k);  % mean-normalize
        end
        
    end
    
    % Sort CCGs based on categories of pairs
    tonicCCR = CCRnorm(pairT,:);  % CCG for tonic cells
    phasicBCCR = CCRnorm(pairB,:);  % CCG for phasic bursting cells
    poissonLCCR = CCRnorm(pairPL,:);  % CCG for phasic non-bursting cells
    earlyFiringCCR = CCRnorm(pairEF,:); % CCG for early firing cells
    tContCCR = CCRnorm(pairTcont,:);% CCG for tonic-containing cells
    pContCCR = CCRnorm(pairPcont,:);% CCG for 'Poisson-containing' cells
    B_PL_CCR = CCRnorm(pairB_PL,:);  % CCG for B_PL cells
    B_T_CCR = CCRnorm(pairB_T,:);  % CCG for B_T cells
    PL_T_CCR = CCRnorm(pairPL_T,:);  % CCG for PL_T cells
    
    mn_tonic = median(tonicCCR, 1);   % average CCG, tonic cells
    mn_tonic = normalizem(mn_tonic,centerP);
    se_tonic = se_of_median(tonicCCR);   % SE, tonic cells
    se_tonic = normalizem(se_tonic,centerP);
    mn_phasicB = median(phasicBCCR, 1);   % average CCG, phasic, bursting cells
    mn_phasicB = normalizem(mn_phasicB,centerP);
    se_phasicB = se_of_median(phasicBCCR, [], 1);   % SE, phasic, bursting cells
    se_phasicB = normalizem(se_phasicB,centerP);
    mn_poissonL = median(poissonLCCR, 1);   % average CCG, phasic, non-bursting cells
    mn_poissonL = normalizem(mn_poissonL,centerP);
    se_poissonL = se_of_median(poissonLCCR, [], 1);   % SE, non-bursting cells
    se_poissonL = normalizem(se_poissonL,centerP);
    mn_B_PL = median(B_PL_CCR, 1);   % average CCG, phasicB-poissonL pairs
    mn_B_PL = normalizem(mn_B_PL,centerP);
    se_B_PL = se_of_median(B_PL_CCR, [], 1);   % SE, phasicB-poissonL pairs
    se_B_PL = normalizem(se_B_PL,centerP);
    mn_B_T = median(B_T_CCR, 1);   % average CCG, phasicB-tonic pairs
    mn_B_T = normalizem(mn_B_T,centerP);
    se_B_T = se_of_median(B_T_CCR, [], 1);   % SE, phasicB-tonic pairs
    se_B_T = normalizem(se_B_T,centerP);
    mn_PL_T = median(PL_T_CCR, 1);   % average CCG, poissonL-tonic pairs
    mn_PL_T = normalizem(mn_PL_T,centerP);
    se_PL_T = se_of_median(PL_T_CCR, [], 1);   % SE, poissonL-tonic pairs
    se_PL_T = normalizem(se_PL_T,centerP);
    mn_EF = median(earlyFiringCCR, 1);   % average CCG, EF pairs
    mn_EF = normalizem(mn_EF,centerP);
    se_EF = se_of_median(earlyFiringCCR, [], 1);   % SE, EF pairs
    se_EF = normalizem(se_EF,centerP);
    mn_tCont = median(tContCCR, 1);   % average CCG, Tonic containing pairs
    mn_tCont = normalizem(mn_tCont,centerP);
    se_tCont = se_of_median(tContCCR, [], 1);   % SE, phasic, Tonic containing pairs
    se_tCont = normalizem(se_tCont,centerP);
    mn_pCont = median(pContCCR, 1);   % average CCG, Poisson containing pairs
    mn_pCont = normalizem(mn_pCont,centerP);
    se_pCont = se_of_median(pContCCR, [], 1);   % SE, phasic, Poisson containing pairs
    se_pCont = normalizem(se_pCont,centerP);
    
    
    numB = sum(pairB);
    numPL = sum(pairPL);
    numT = sum(pairT);
    numB_PL = sum(pairB_PL);
    numB_T = sum(pairB_T);
    numPL_T = sum(pairPL_T);
    numPcont = sum(pairPcont);
    numTcont = sum(pairTcont);
    numEF = sum(pairEF);
    
    % REVIEWER FIGURE 3
    % Plot EF and Tonic containing
    H6 = figure;
    hold on
    errorshade(lags,mn_phasicB,se_phasicB,'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',2)
    errorshade(lags,mn_pCont,se_pCont,'LineColor',cCode(2,:),'ShadeColor',cCode(2,:),'LineWidth',2)
    errorshade(lags,mn_tCont,se_tCont,'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',2)
    errorshade(lags,mn_all,se_all,'LineColor',[0.6 0.6 0.6],'ShadeColor',[0.6 0.6 0.6],'LineWidth',2)
    xlim([plotWin(1) plotWin(2)]);
    xticks([plotWin(1) 0 plotWin(2)]);
    xlabel('Lag (ms)')
    ylabel('Normalized count')
    title(['B, P-cont, T-cont, untagged cells ' bfilt ' ' tetTypes{iP}])
    plotStr= {['BURST-SB n=',num2str(numB)];
        ['Pcont n=',num2str(numPcont)];
        ['Tcont n=',num2str(numTcont)];
        ['Untagged n=',num2str(numAll)]};
    legend(plotStr, 'FontSize', 12, 'Location','northeast')
    legend('boxoff');
    setmyplot_tamas;
    axis square;
    fnm = ['B_PC_TC_UT_' burstfilter '_' cb '_' tetTypes{iP} '.fig'];
    fnm2 = ['B_PC_TC_UT_' burstfilter '_' cb '_' tetTypes{iP} '.jpeg'];
    saveas(H6,fullfile([resdir1  PATH pairPath  '\'],fnm))   % save plot
    saveas(H6,fullfile([resdir1  PATH pairPath  '\'],fnm2))   % save jpeg
    close(H6);
    
    % Activation Index
    ccrB = CCR(pairB,:);   % from non-normalized CCG
    ccrP = CCR(pairPL,:);   % from non-normalized CCG
    ccrT = CCR(pairT,:);   % from non-normalized CCG
    ccrTC = CCR(pairTcont,:);
    ccrPC = CCR(pairPcont,:);
    actRatioB = mean(ccrB(:, (centerP-30):(centerP+30)),2);
    baseRatioB = mean(ccrB(:, (centerP+100):end),2);
    actRatioP = mean(ccrP(:, (centerP-30):(centerP+30)),2);
    baseRatioP = mean(ccrP(:, (centerP+100):end),2);
    actRatioT = mean(ccrT(:, (centerP-30):(centerP+30)),2);
    baseRatioT = mean(ccrT(:, (centerP+100):end),2);
    actRatioTC = mean(ccrTC(:, (centerP-30):(centerP+30)),2);
    baseRatioTC = mean(ccrTC(:, (centerP+100):end),2);
    actRatioPC = mean(ccrPC(:, (centerP-30):(centerP+30)),2);
    baseRatioPC = mean(ccrPC(:, (centerP+100):end),2);
    actRatio = mean(CCR(:, (centerP-30):(centerP+30)),2);
    baseRatio = mean(CCR(:, (centerP+100):end),2);
    ratioB = (actRatioB-baseRatioB)./(actRatioB+baseRatioB);
    ratioP = (actRatioP-baseRatioP)./(actRatioP+baseRatioP);
    ratioT = (actRatioT-baseRatioT)./(actRatioT+baseRatioT);
    ratioTC = (actRatioTC-baseRatioTC)./(actRatioTC+baseRatioTC);
    ratioPC = (actRatioPC-baseRatioPC)./(actRatioPC+baseRatioPC);
    ratioALL = (actRatio-baseRatio)./(actRatio+baseRatio);
    
    % Boxplot
    [H7, Wp_B_TC] = boxstat(ratioB,ratioTC,'Bursting','Tonic-containing',0.05,'nonpaired');
    setmyplot_tamas;
    storeP.B_TC = round(Wp_B_TC,4);
    close(H7);
    [H7, Wp_B_PC] = boxstat(ratioB,ratioPC,'Bursting','Poisson-containing',0.05,'nonpaired');
    setmyplot_tamas;
    storeP.B_PC = round(Wp_B_PC,4);
    close(H7);
    [H7, Wp_PC_TC] = boxstat(ratioPC,ratioTC,'Poisson-containing','Tonic-containing',0.05,'nonpaired');
    setmyplot_tamas;
    storeP.PC_TC = round(Wp_PC_TC,4);
    close(H7);
    
    [H7, Wp_B_UT] = boxstat(ratioB,ratioUT,'Bursting','Untagged',0.05,'nonpaired');
    setmyplot_tamas;
    storeP.B_UT = round(Wp_B_UT,4);
    close(H7);
    [H7, Wp_PC_UT] = boxstat(ratioPC,ratioUT,'Poisson-containing','Untagged',0.05,'nonpaired');
    setmyplot_tamas;
    storeP.PC_UT = round(Wp_PC_UT,4);
    close(H7);
    [H7, Wp_TC_UT] = boxstat(ratioTC,ratioUT,'Tonic-containing','Untagged',0.05,'nonpaired');
    setmyplot_tamas;
    storeP.TC_UT = round(Wp_TC_UT,4);
    close(H7);

    
    coords = [1, 2, 3.5, 5];
    
    % Significance test
    [groups, pAll, pVals] = sigTest(ratioALL, ratioUT, pairB, pairPcont, pairTcont,...
        true(numel(ratioUT),1), coords);
    
    % x coordinate generation to scatter the plot
    x1 = linspace(coords(1)-0.2, coords(1)+0.2, sum(pairB));
    x2 = linspace(coords(2)-0.2, coords(2)+0.2, sum(pairPcont));
    x3 = linspace(coords(3)-0.2, coords(3)+0.2, sum(pairTcont));
    
    % FIGURE 4 PANEL B
    % Bargraph of EF and Tonic containing group - MEDIAN
    H8 = figure;
    hold on;
    plot(x1, ratioB, 'o', 'MarkerFaceColor','none','MarkerEdgeColor',cCode(1,:), 'LineWidth',3)
    plot(x2, ratioPC, 'o', 'MarkerFaceColor','none','MarkerEdgeColor',cCode(2,:), 'LineWidth',3)
    plot(x3, ratioTC, 'o', 'MarkerFaceColor','none','MarkerEdgeColor',cCode(3,:), 'LineWidth',3)
    bar(coords, [median(ratioB) median(ratioPC) median(ratioTC) median(ratioUT)], 'FaceColor', 'none', 'LineWidth',2)
    x_lim = xlim;
    y_lim = ylim;
    sigstar(groups,pAll);
    plotStr= {['BURST-SB n=',num2str(numB)];
        ['Pcont n=',num2str(numPcont)];
        ['Tcont n=',num2str(numTcont)];
        ['Untagged n=',num2str(numAll)]
        };
    annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
    'EdgeColor', 'none');
    xticks([]);
    setmyplot_tamas;
    axis square;
    title(['B, P-cont and T-cont cells ' bfilt ' ' tetTypes{iP}])
    ylabel('Synchrony Index')
    fnm = ['B_PC_TC_TPairs_' burstfilter '_' cb '_' tetTypes{iP} 'barplot_median.fig'];
    fnm2 = ['B_PC_TC_TPairs_' burstfilter '_' cb '_' tetTypes{iP} 'barplot_median.jpeg'];
    saveas(H8,fullfile([resdir1  PATH pairPath  '\'],fnm))   % save plot
    saveas(H8,fullfile([resdir1  PATH pairPath  '\'],fnm2))   % save jpeg
    close(H8);
    
    
    ratios.ratioB.(tetTypes{iP}) = {ratioB};
    ratios.ratioPC.(tetTypes{iP}) = {ratioPC};
    ratios.ratioTC.(tetTypes{iP}) = {ratioTC};
    ratios.ratioAll.(tetTypes{iP}) = ratioALL;
    ratios.pairB.(tetTypes{iP}) = pairB;
    ratios.pairPC.(tetTypes{iP}) = pairPcont;
    ratios.pairTC.(tetTypes{iP}) = pairTcont;
    ratios.ttInx.(tetTypes{iP}) = currTT;
    
    % IMAGESC'S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIGURE 4 PANEL A PLOT
    % Plot CCG image sorted by activation index (bursting, then rest)
    meanAct = median(CCRnorm(:, (centerP-30):(centerP+30)),2);
    [~,index3] = sortrows(ratioALL,'descend');
    
    % Plot image
    H10 = figure;
    imagesc(lags(lags>=plotWin(1)&lags<=plotWin(2)), 1:size(CCRnorm,1), CCRnorm(index3, lags>=plotWin(1)&lags<=plotWin(2)))
    yticks = 1:length(CCRnorm);
    groupids = groupID(index3);
    yMatrix = nan(numPair,2);
    yMatrix(cellfun(@(s)strcmp(s,'phasicB'),group1(index3)),1) = 1;
    yMatrix(cellfun(@(s)strcmp(s,'poissonL'),group1(index3)),1) = 2;
    yMatrix(cellfun(@(s)strcmp(s,'tonic'),group1(index3)),1) = 3;
    yMatrix(cellfun(@(s)strcmp(s,'phasicB'),group2(index3)),2) = 1;
    yMatrix(cellfun(@(s)strcmp(s,'poissonL'),group2(index3)),2) = 2;
    yMatrix(cellfun(@(s)strcmp(s,'tonic'),group2(index3)),2) = 3;
    set(gca,'ytick',yticks);
    yticklabels({' '});
    title(['Normalized CCG of Cholinergic pairs ' tetTypes{iP}])
    xlabel('Lag (ms)')
    ylabel('Group ID')
    colormap jet
    % caxis([0 3]);
    xticks([plotWin(1) 0 plotWin(2)])
    fnm = ['sorted_normalized_avgCCG_allGroup_' burstfilter '_' cb '_' tetTypes{iP} '_allRatio.fig'];
    fnm2 = ['sorted_normalized_avgCCG_allGroup_' burstfilter '_' cb '_' tetTypes{iP} '_allRatio.jpeg'];
    axis square;
    colorbar;
    setmyplot_tamas;
    saveas(H10,fullfile([resdir1  PATH pairPath  '\'],fnm))   % save plot
    saveas(H10,fullfile([resdir1  PATH pairPath  '\'],fnm2))   % save jpeg
    close(H10);
    
    % Colormap plot for replacing yLabels
    map = [cCode(1,:); cCode(2,:); cCode(3,:)];
    H11 = figure;
    imagesc(yMatrix);
    colormap(map);
    line([1.5 1.5], ylim, 'Color', [0 0 0], 'LineWidth',1)
    for i = 1:length(groupids)
        line([0.5 2.5], [i+0.5 i+0.5], 'Color', [0 0 0], 'LineWidth',2)
    end
    axis off;
    axis square;
    fnm = ['colormap_boxes_' burstfilter '_' cb '_' tetTypes{iP} '.fig'];
    fnm2 = ['colormap_boxes_' burstfilter '_' cb '_' tetTypes{iP} '.jpeg'];
    saveas(H11,fullfile([resdir1  PATH pairPath  '\'],fnm))   % save plot
    saveas(H11,fullfile([resdir1  PATH pairPath  '\'],fnm2))   % save jpeg
    close(H11);
end
keyboard;

% Extra
%% Compare tetrode and non-tetrodepairs

ratioB_nt = ratios.ratioB.nontetrodepairs;
ratioPC_nt = ratios.ratioPC.nontetrodepairs;
ratioTC_nt = ratios.ratioTC.nontetrodepairs;
ratioB_t = ratios.ratioB.tetrodepairs;
ratioPC_t = ratios.ratioPC.tetrodepairs;
ratioTC_t = ratios.ratioTC.tetrodepairs;

ratioALL = ratios.ratioAll.all;
pairB_nt = ratios.pairB.nontetrodepairs;
pairPcont_nt = ratios.pairPC.nontetrodepairs;
pairTcont_nt = ratios.pairTC.nontetrodepairs;
pairB_t = ratios.pairB.tetrodepairs;
pairPcont_t = ratios.pairPC.tetrodepairs;
pairTcont_t = ratios.pairTC.tetrodepairs;

tInx = ratios.ttInx.tetrodepairs;
ntInx = ratios.ttInx.nontetrodepairs;

pairBT = tInx' & ratios.pairB.all;
pairBNT = ntInx' & ratios.pairB.all;
pairPCT = tInx' & ratios.pairPC.all;
pairPCNT = ntInx' & ratios.pairPC.all;
pairTCT = tInx' & ratios.pairTC.all;
pairTCNT = ntInx' & ratios.pairTC.all;

coords = [1, 2, 3, 4, 5.5, 6.5];

% Significance test
% [groups, pAll, pVals] = sigTest(ratioALL, pairBT, pairBNT, pairPCT,...
%     pairPCNT, pairTCT, pairTCNT, coords);

% x coordinate generation to scatter the plot
x1 = linspace(coords(1)-0.2, coords(1)+0.2, sum(pairBT));
x2 = linspace(coords(2)-0.2, coords(2)+0.2, sum(pairBNT));
x3 = linspace(coords(3)-0.2, coords(3)+0.2, sum(pairPCT));
x4 = linspace(coords(4)-0.2, coords(4)+0.2, sum(pairPCNT));
x5 = linspace(coords(5)-0.2, coords(5)+0.2, sum(pairTCT));
x6 = linspace(coords(6)-0.2, coords(6)+0.2, sum(pairTCNT));



% Bargraph of EF and Tonic containing group - MEDIAN
H8 = figure;
hold on;
bar(coords, [median(ratioALL(pairBT)) median(ratioALL(pairBNT))...
    median(ratioALL(pairPCT)) median(ratioALL(pairPCNT))...
    median(ratioALL(pairTCT)) median(ratioALL(pairTCNT))], 'FaceColor', 'w', 'LineWidth',2)
plot(x1, ratioALL(pairBT), 'o', 'MarkerFaceColor',cCode(1,:),'MarkerEdgeColor',cCode(1,:), 'LineWidth',3)
plot(x2, ratioALL(pairBNT), 'o', 'MarkerFaceColor',cCode(1,:),'MarkerEdgeColor',cCode(1,:), 'LineWidth',3)
plot(x3, ratioALL(pairPCT), 'o', 'MarkerFaceColor',cCode(2,:),'MarkerEdgeColor',cCode(2,:), 'LineWidth',3)
plot(x4, ratioALL(pairPCNT), 'o', 'MarkerFaceColor',cCode(2,:),'MarkerEdgeColor',cCode(2,:), 'LineWidth',3)
plot(x5, ratioALL(pairTCT), 'o', 'MarkerFaceColor',cCode(3,:),'MarkerEdgeColor',cCode(3,:), 'LineWidth',3)
plot(x6, ratioALL(pairTCNT), 'o', 'MarkerFaceColor',cCode(3,:),'MarkerEdgeColor',cCode(3,:), 'LineWidth',3)

x_lim = xlim;
y_lim = ylim;
% text((x_lim(2)*0.5), y_lim(2)*0.95, ['p B vs PC: ' num2str(pVals(1))]);
% text((x_lim(2)*0.5), y_lim(2)*0.85, ['p B vs TC: ' num2str(pVals(2))]);
% text((x_lim(2)*0.5), y_lim(2)*0.75, ['p PC vs TC: ' num2str(pVals(3))]);
% sigstar(groups,pAll);
% plotStr= {['BURST-SB n=',num2str(numB)];
%     ['Pcont n=',num2str(numPcont)];
%     ['Tcont n=',num2str(numTcont)]
%     };
% annotation('textbox',[.54 .55 .3 .3],'String',plotStr,'FitBoxToText','on',...
%     'EdgeColor', 'none');
setmyplot_tamas;
axis square;
title('Tetrode vs non-tetrodepairs')
ylabel('Synchrony index')
xticks([])
fnm = ['TETRODEPAIRS_' burstfilter '_' cb '_' tetTypes{iP} 'barplot_median.fig'];
fnm2 = ['TETRODEPAIRS_' burstfilter '_' cb '_' tetTypes{iP} 'barplot_median.jpeg'];
saveas(H8,fullfile([resdir1  PATH pairPath  '\'],fnm))   % save plot
saveas(H8,fullfile([resdir1  PATH pairPath  '\'],fnm2))   % save jpeg
close(H8);


%--------------------------------------------------------------------------
function [groupID, yLabel, pairPath, chInx] = categorizeTP(PairOfCells, group1, group2, pairTypes)
% Categorizes the cellpairs based on their TP identity
%OUTPUT: groupID of the pair, and the corresponding Label for plotting

% Subselect pair types inside the current cellbase, celltypes
switch pairTypes
    case 'normal' % Keeps all the cellids selected earlier for CCG
        % All cellids kept
        pairPath = '';
        chInx = 1:length(group1);
    case 'tagged' % only tagged
        ChAT1 = getvalue('ChAT+',PairOfCells(:,1));
        ChAT2 = getvalue('ChAT+',PairOfCells(:,2));
        bothCh = (ChAT1>0 & ChAT2>0);
        PairOfCells = PairOfCells(bothCh,:);
        chInx = find(bothCh);
        pairPath = '\tagged';
    case 'putative' % only putative
        pChAT1 = getvalue('pChAT+',PairOfCells(:,1));
        pChAT2 = getvalue('pChAT+',PairOfCells(:,2));
        bothpCh = (pChAT1>0 & pChAT2>0);
        chInx = find(bothpCh);
        PairOfCells = PairOfCells(bothpCh,:);
        pairPath = '\putative';
    case 'untagged' % only unidentified cells
        ChAT1 = getvalue('ChAT+',PairOfCells(:,1));
        ChAT2 = getvalue('ChAT+',PairOfCells(:,2));
        pChAT1 = getvalue('pChAT+',PairOfCells(:,1));
        pChAT2 = getvalue('pChAT+',PairOfCells(:,2));
        notCh = (ChAT1==0 & ChAT2==0 & pChAT1==0 & pChAT2==0);
        chInx = find(notCh);
        PairOfCells = PairOfCells(notCh,:);
        pairPath = '\untagged';
    case 'untaggedMix' % Pairs of an identified and a cholinergic cell
        ChAT1 = getvalue('ChAT+',PairOfCells(:,1));
        ChAT2 = getvalue('ChAT+',PairOfCells(:,2));
        pChAT1 = getvalue('pChAT+',PairOfCells(:,1));
        pChAT2 = getvalue('pChAT+',PairOfCells(:,2));
        notCh = (ChAT1==0 & ChAT2==0 & pChAT1==0 & pChAT2==0);
        bothCh = ((ChAT1>0 | pChAT1>0) & (ChAT2>0 | pChAT2>0));
        mixCh = sort([find(bothCh); find(notCh)]);
        chInx = find(mixCh);
        PairOfCells(mixCh,:) = [];
        pairPath = '\mixACh';
end

% Preallocate
groupID = zeros(1,size(PairOfCells,1));
yLabel=cell(1, size(PairOfCells,1));

% Categorize Pairs to TP groups
for numPair = 1:size(PairOfCells,1)
    currgroup1 = group1{numPair};
    currgroup2 = group2{numPair};
    [grouping, Label] = groupCellPair(currgroup1, currgroup2);
    groupID(numPair) = grouping;
    yLabel{1,numPair} = Label;
end

% -------------------------------------------------------------------------
function nM = normalizem(M,cp)

sM = smooth(M,'linear',15);
bsl = mean(sM(cp+100:end));
nM = sM - bsl;

%--------------------------------------------------------------------------
function [groups, pAll, p] = sigTest(currData, currData2, group1, group2,... 
        group3, group4, coords)

% Significance test
% Boxplots
[H1, p1] = boxstat(currData(group1),currData(group2),'','',0.05,'nonpaired');
close(H1);
[H2, p2] = boxstat(currData(group1),currData(group3),'','',0.05,'nonpaired');
close(H2);
[H3, p3] = boxstat(currData(group2),currData(group3),'','',0.05,'nonpaired');
close(H3);
[H4, p4] = boxstat(currData(group1),currData2(group4),'','',0.05,'nonpaired');
close(H4);
[H5, p5] = boxstat(currData(group2),currData2(group4),'','',0.05,'nonpaired');
close(H5);
[H6, p6] = boxstat(currData(group3),currData2(group4),'','',0.05,'nonpaired');
close(H6);

p(1) = round(p1,3);
p(2) = round(p2,3);
p(3) = round(p3,3);
p(4) = round(p4,3);
p(5) = round(p5,3);
p(6) = round(p6,3);

groupInc = [p1<=0.05, p2<=0.05, p3<=0.05, p4<=0.05, p5<=0.05, p6<=0.05];

groups={[coords(1),coords(2)],[coords(1),coords(3)],...
    [coords(2),coords(3)],[coords(1),coords(4)],...
    [coords(2),coords(4)],[coords(3),coords(4)]};
groups = groups(groupInc);
pAll = p(groupInc);


%--------------------------------------------------------------------------
function [groups, pAll, p] = sigTest2(currData, group1, group2, group3,...
    group4, group5, group6, coords)

% Significance test
% Boxplots
[H1, p1] = boxstat(currData(group1),currData(group2),'','',0.05,'nonpaired');
close(H1);
[H2, p2] = boxstat(currData(group1),currData(group3),'','',0.05,'nonpaired');
close(H2);
[H3, p3] = boxstat(currData(group2),currData(group3),'','',0.05,'nonpaired');
close(H3);
[H4, p4] = boxstat(currData(group2),currData(group3),'','',0.05,'nonpaired');
close(H4);
[H5, p5] = boxstat(currData(group2),currData(group3),'','',0.05,'nonpaired');
close(H5);
[H6, p6] = boxstat(currData(group2),currData(group3),'','',0.05,'nonpaired');
close(H6);

p(1) = round(p1,3);
p(2) = round(p2,3);
p(3) = round(p3,3);
p(4) = round(p4,3);
p(5) = round(p5,3);
p(6) = round(p6,3);

groupInc = [p1<=0.05, p2<=0.05, p3<=0.05, p4<=0.05, p5<=0.05,...
    p6<=0.05];

groups={[coords(1),coords(2)],[coords(1),coords(3)],...
    [coords(2),coords(3)]};
groups = groups(groupInc);
pAll = p(groupInc);

