function eventDataPlot_All(cb, actBurst, t_Type, c_Type, f_Type, doraster)
%PSTHEXTRACT   Plot function for psthExtract data
%   EVENTDATAPLO_ALL(CB, actBurst) Plotter function for psthExtract
%
%   See also PSTHEXTRACT.

%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Initilize
dbstop if error;
global RESDIR;
global cCode;
global PATH;
fs = filesep;

% Filters
adapFilt = 'NotAdaptive';
psthFilt = 'spsth';

if actBurst
    outLType = '_removed';
    outText = 'removed';
else
    outLType = [];
    outText = 'original';
end

% Define directories
resdir = [RESDIR cb fs];   % results directory
psthDir = [resdir 'psth' cb fs 'figures' PATH fs];
dataDir = [RESDIR cb fs 'PSTH_newdata' cb PATH fs];
listDir = dir([dataDir 'EventData*']);

% Load all psth the data for every cellbase
% Find latest EventData file
edFiles = {listDir.name};

if ~isempty(edFiles)
    if actBurst
        remInx = find(~cellfun('isempty',strfind(edFiles,'removed')));
    else
        remInx = find(cellfun('isempty',strfind(edFiles,'removed')));
    end
    currFiles = {edFiles{remInx}};
    dates = [listDir(remInx).datenum];
    latestD = find(dates==max(dates));
    allData = load([dataDir currFiles{latestD}]);
    
else
    % Pool data from separate cb source files
    allData = loadSeparate(actBurst);
end


if strcmp(PATH(2:end), 'CPu')
    pLimits = [-250 750];
    xtick_labels = {'-250', '0', '250', '750'};
    x_ticks = [-250 0 250 750];
else
    pLimits = [-500 500];
    xtick_labels = {'0', '40', '80'};
    x_ticks = [0 40 80];
end

fieldNames = fieldnames(allData);
if contains(fieldNames, 'EventData')
    EventData = allData.EventData;
    indData = true; % individual cellbase
else
    EventData = allData;
    indData = false; % pooled cellbase
end

% Define plotting parameters
if indData
    centerP = round(size(EventData.FB.FA.NotAdaptive.none.(psthFilt),2)/2);
else
    centerP = round(size(EventData.FB.FA.none.(psthFilt),2)/2);
end
time = pLimits(1):1:pLimits(2); % time vector for plotting
plotWindow = (centerP+time(1)):(centerP+time(end));
actWindow = centerP:(centerP+100); % activation win for burstratio
baseWindow = (centerP+100):plotWindow(end); % base win for burstratio
LWidth = 3;
plotTime = [-20 80];

dataNames = fieldnames(EventData);

% Clear empty cells and Normalize
for iD = 1:length(dataNames)
    tName = fieldnames(EventData.(dataNames{iD}));
    for iT = 1:numel(tName)
        switch tName{iT}
            case 'Miss'
                currData = EventData.(dataNames{iD}).Miss;
            case 'CorrectRejection'
                currData = EventData.(dataNames{iD}).CorrectRejection;
            otherwise
                currData = EventData.(dataNames{iD}).(tName{iT});
        end
        
        currFields = fieldnames(currData);
        if contains(currFields, 'NotAdaptive')
            currData = currData.NotAdaptive;
        elseif contains(currFields, 'fa')
            currData = currData.fa.NotAdaptive;
        elseif contains(currFields, 'hit')
            currData = currData.hit.NotAdaptive;
        end

        [currData, ~] = clearEmpty(currData); % Removes empty cells
        normData.(dataNames{iD}).(tName{iT}) = select_normalization(currData);
        allChAT.(dataNames{iD}) = currData.none.cellid;
    end
end

% Store cellids and groups
for iCB=1:numel(dataNames)
    numChAT = numel(allChAT.(dataNames{iCB}));
    fieldName = fieldnames(EventData.(dataNames{iCB}));
    groups = [];
    groupID = [];
    groups = normData.(dataNames{iCB}).(fieldName{1}).none.stats;
    
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
    % TP groups
    phasicB = groupID==1;
    poissonL = groupID==2;
    tonic = groupID==3;
    
    % Preallocate
    nbNum = 0;
    hdbNum = 0;
    phdbNum = 0;
    cbUsed = [0 0 0];
    ChAT = [];
    pChAT = [];
    % Define ChAT, pChAT properties for all CBs
    for iCC = 1:numChAT
        cbSwitcher(allChAT.(dataNames{iCB})(iCC))
        currCB = whichcb;
        
        switch currCB
            case 'NB'
                ChAT(iCC) = logical(getvalue('ChAT+',allChAT.(dataNames{iCB})(iCC)));
                pChAT(iCC) = logical(getvalue('pChAT+',allChAT.(dataNames{iCB})(iCC)));
                nbNum = iCC;
                cbUsed(1) = 1;
            case 'HDB'
                ChAT(iCC) = true;
                pChAT(iCC) = false;
                hdbNum = iCC;
                cbUsed(2) = 1;
            case 'PannaHDB'
                ChAT(iCC) = true;
                pChAT(iCC) = false;
                phdbNum = iCC;
                cbUsed(3) = 1;
        end
    end
    
    
    if sum(cbUsed)>1 % Multiple cb-s pooled together
        cellNums = sort([nbNum, hdbNum, phdbNum], 'ascend');
        typeNums = [cellNums(1) diff(cellNums)];
        numNB = typeNums(1);
        numHDB = typeNums(2);
        numPHDB = typeNums(3);
        % ChAT
        phasiC.(dataNames{iCB}).ChAT = [ChAT(1:numNB) & phasicB(1:numNB),...
            phasicB((numNB+1):end)]';
        poissoN.(dataNames{iCB}).ChAT = [ChAT(1:numNB) & poissonL(1:numNB),...
            poissonL((numNB+1):end)]';
        toniC.(dataNames{iCB}).ChAT = [ChAT(1:numNB) & tonic(1:numNB),...
            tonic((numNB+1):end)]';
        % pChAT
        phasiC.(dataNames{iCB}).pChAT = [pChAT(1:numNB) & phasicB(1:numNB),...
            phasicB((numNB+1):end)==2]';
        poissoN.(dataNames{iCB}).pChAT = [pChAT(1:numNB) & poissonL(1:numNB),...
            poissonL((numNB+1):end)==2]';
        toniC.(dataNames{iCB}).pChAT = [pChAT(1:numNB) & tonic(1:numNB),...
            tonic((numNB+1):end)==2]';
        % all
        phasiC.(dataNames{iCB}).allC = phasicB';
        poissoN.(dataNames{iCB}).allC = poissonL';
        toniC.(dataNames{iCB}).allC = tonic';
        
    else % individual cb
        phasiC.(dataNames{iCB}).ChAT = ChAT & phasicB;
        poissoN.(dataNames{iCB}).ChAT = ChAT & poissonL;
        toniC.(dataNames{iCB}).ChAT = ChAT & tonic;
        phasiC.(dataNames{iCB}).pChAT = pChAT & phasicB;
        poissoN.(dataNames{iCB}).pChAT = pChAT & poissonL;
        toniC.(dataNames{iCB}).pChAT = pChAT & tonic;
        phasiC.(dataNames{iCB}).allC = phasicB';
        poissoN.(dataNames{iCB}).allC = poissonL';
        toniC.(dataNames{iCB}).allC = tonic';
    end
end
% Putative containing cbs
putativeCont = {'NB', 'POOLED'};
if ismember(cb,putativeCont)
    cTypes = {'ChAT', 'pChAT', 'allC'};
else
    cTypes = {'ChAT', 'allC'};
end

% Filter types
bFilt = {'none', 'burst1', 'single', 'burstall', 'burst1Single'};
fType = {'FA', 'HIT'};
tTypes = fieldnames(normData);

%% Structure normalized data
for iC = 1:length(cTypes)
    for iT=1:length(tTypes)
        switch tTypes{iT}
            case 'MISS'
                fType = {'Miss'};
            case 'CR'
                fType = {'CorrectRejection'};
            otherwise
                fType = {'FA', 'HIT'};
        end
        for iF = 1:length(fType)
            for iB = 1:length(bFilt)
                plotData = normData.(tTypes{iT}).(fType{iF}).(bFilt{iB}).(psthFilt);
                % average ACG, phasic, bursting cells
                mn_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iC}).(bFilt{iB}) =...
                    nanmean(plotData(phasiC.(tTypes{iT}).(cTypes{iC}), plotWindow),1);   
                % SE, phasic, bursting cells
                se_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iC}).(bFilt{iB}) =...
                    std(plotData(phasiC.(tTypes{iT}).(cTypes{iC}),...
                    plotWindow)) / sqrt(sum(phasicB));   
                % average ACG, poisson-like cells
                mn_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iC}).(bFilt{iB}) =...
                    nanmean(plotData(poissoN.(tTypes{iT}).(cTypes{iC}),...
                    plotWindow),1);   
                % SE, poisson-like cells
                se_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iC}).(bFilt{iB}) =...
                    std(plotData(poissoN.(tTypes{iT}).(cTypes{iC}),...
                    plotWindow)) / sqrt(sum(poissonL));   
                % average ACG, tonic cells
                mn_tonic.(tTypes{iT}).(fType{iF}).(cTypes{iC}).(bFilt{iB}) =...
                    nanmean(plotData(toniC.(tTypes{iT}).(cTypes{iC}),...
                    plotWindow),1);   
                % SE, tonic cells
                se_tonic.(tTypes{iT}).(fType{iF}).(cTypes{iC}).(bFilt{iB}) =...
                    std(plotData(toniC.(tTypes{iT}).(cTypes{iC}),...
                    plotWindow)) / sqrt(sum(tonic));   
            end
        end
    end
end

% Subselect based on input arguments
if ~isempty(t_Type)
    tTypes = {t_Type};
end
if ~isempty(f_Type)
    fType = {f_Type};
end
if ~isempty(c_Type)
    cTypes = {c_Type};
end

% FIGURES
for iT = 1:length(tTypes)
    for iF = 1:length(fType)
        for iP = 1:length(cTypes)
            %% calculate selectivity data
            phasicC = [];
            poissonN = [];
            burst1Data = [];
            singleData = [];
            burstInd = [];
            singleInd = [];
            
            burst1Data = normData.(tTypes{iT}).(fType{iF}).burst1.(psthFilt);
            singleData = normData.(tTypes{iT}).(fType{iF}).single.(psthFilt);
            
            % For phasicB cells
            phasicC.ChAT = phasiC.(tTypes{iT}).ChAT;
            phasicC.pChAT = phasiC.(tTypes{iT}).pChAT;
            phasicC.allC = phasiC.(tTypes{iT}).allC;
            % For poisson cells
            poissonN.ChAT = poissoN.(tTypes{iT}).ChAT;
            poissonN.pChAT = poissoN.(tTypes{iT}).pChAT;
            poissonN.allC = poissoN.(tTypes{iT}).allC;
            
            %PhasicB
            %BurstIndex (Act. window mean - base win mean / act+base win means)
            burstInd = (mean(burst1Data(phasicC.(cTypes{iP}),actWindow),2)-...
                mean(burst1Data(phasicC.(cTypes{iP}),baseWindow),2))./...
                (mean(burst1Data(phasicC.(cTypes{iP}),actWindow),2)+...
                mean(burst1Data(phasicC.(cTypes{iP}),baseWindow),2));
            
            singleInd = (mean(singleData(phasicC.(cTypes{iP}),actWindow),2)-...
                mean(singleData(phasicC.(cTypes{iP}),baseWindow),2))./...
                (mean(singleData(phasicC.(cTypes{iP}),actWindow),2)+...
                mean(singleData(phasicC.(cTypes{iP}),baseWindow),2));
            
            
            %PoissonL
            %BurstIndex (Act. window mean - base win mean / act+base win means)
            burstIndP = (mean(burst1Data(poissonN.(cTypes{iP}),actWindow),2)-...
                mean(burst1Data(poissonN.(cTypes{iP}),baseWindow),2))./...
                (mean(burst1Data(poissonN.(cTypes{iP}),actWindow),2)+...
                mean(burst1Data(poissonN.(cTypes{iP}),baseWindow),2));
            
            singleIndP = (mean(singleData(poissonN.(cTypes{iP}),actWindow),2)-...
                mean(singleData(poissonN.(cTypes{iP}),baseWindow),2))./...
                (mean(singleData(poissonN.(cTypes{iP}),actWindow),2)+...
                mean(singleData(poissonN.(cTypes{iP}),baseWindow),2));

            
            bInx.(tTypes{iT}).(fType{iF}).(cTypes{iP}) = burstInd;
            sInx.(tTypes{iT}).(fType{iF}).(cTypes{iP}) = singleInd;
            
            % Eliminate bursting cells which has 0 bursts in the activation window
            nanInx = [];
            outlierInx = [];
            
            % PhasicB
            nanInx = ~isnan(burstInd);
            outlierInx = find(~nanInx);
            burstInd = burstInd(nanInx);
            singleInd = singleInd(nanInx);
            
            phasicCells = [];
            nonActCells = [];
            phasicCells = normData.(tTypes{iT}).(fType{iF}).burst1...
                .cellid(phasicC.(cTypes{iP}));
            nonActCells = phasicCells(~nanInx);
            
            % PoissonL
            nanInxP = ~isnan(burstIndP);
            outlierInxP = find(~nanInxP);
            burstIndP = burstIndP(nanInxP);
            singleIndP = singleIndP(nanInxP);
            
            poissonCells = [];
            nonActCellsP = [];
            poissonCells = normData.(tTypes{iT}).(fType{iF}).burst1...
                .cellid(poissonN.(cTypes{iP}));
            nonActCellsP = poissonCells(~nanInxP);
            
            if ~isempty(nonActCells)
                % Save it and re-run psthExtract
                removeCellDir = [resdir 'psth' cb fs 'removeCells' fs];
                save([removeCellDir 'removeCells_' tTypes{iT} '_' cTypes{iP}...
                    '_' fType{iF} '.mat'],...
                    'nonActCells', 'outlierInx');
            end
            
            % FIGURE3 PANEL D
            % SUPPLEMENTARY FIGURE1 PANEL B, C
            % Figure for Raw
            H1 = figure;
            set(gcf, 'Renderer', 'painters');
            hold on;
            if ~strcmp(PATH(2:end), 'CPu')
                errorshade(time, mn_phasicB.(tTypes{iT}).(fType{iF})...
                    .(cTypes{iP}).none-min(mn_phasicB.(tTypes{iT}).(fType{iF})...
                    .(cTypes{iP}).none),...
                    se_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iP}).none,...
                    'LineColor', cCode(1,:), 'ShadeColor', cCode(1,:),...
                    'LineWidth', LWidth);
                errorshade(time, mn_poissonL.(tTypes{iT}).(fType{iF})...
                    .(cTypes{iP}).none-min(mn_poissonL.(tTypes{iT}).(fType{iF})...
                    .(cTypes{iP}).none), se_poissonL.(tTypes{iT}).(fType{iF})...
                    .(cTypes{iP}).none, 'LineColor', cCode(2,:),...
                    'ShadeColor', cCode(2,:), 'LineWidth', LWidth);
                legText = {'BURST-SB','BURST-PL','REG'};
            else
                legText = {'REG'};
            end
            errorshade(time, mn_tonic.(tTypes{iT}).(fType{iF}).(cTypes{iP})...
                .none-min(mn_tonic.(tTypes{iT}).(fType{iF}).(cTypes{iP}).none),...
                se_tonic.(tTypes{iT}).(fType{iF}).(cTypes{iP}).none,...
                'LineColor', cCode(3,:), 'ShadeColor', cCode(3,:), 'LineWidth',...
                LWidth);
            title(['All AP ' tTypes{iT} ' ' cTypes{iP} ' ' fType{iF} ' ' outText])
            legend(legText, 'FontSize', 12, 'Location','northeast')
            legend('boxoff');
            setmyplot_tamas;
            xlim(plotTime);
            xticks(x_ticks);
            xticklabels(xtick_labels);
            switch fType{iF}
                case 'FA'
                    alignN = 'punishment';
                case 'HIT'
                    alignN = 'reward';
                case 'Miss'
                    alignN = 'feedback';
                case 'CorrectRejection'
                    alignN = 'feedback';
            end
            xlabel(['Time from ' alignN ' (ms)']);
            ylabel('Firing rate (Hz)');
            legend({['n=' num2str(sum(phasiC.(tTypes{iT}).(cTypes{iP})))],['n='...
                num2str(sum(poissoN.(tTypes{iT}).(cTypes{iP})))],...
                ['n=' num2str(sum(toniC.(tTypes{iT}).(cTypes{iP})))]},...
                'FontSize', 12, 'Location','northeast');
            axis square;
            fName = [psthDir 'raw_psth_' tTypes{iT} '_' cTypes{iP} '_'...
                fType{iF} '_' outText '.fig'];   % save PSTH figure
            fNameJ = [psthDir 'raw_psth_' tTypes{iT} '_' cTypes{iP}  '_'...
                fType{iF} '_' outText '.jpeg'];   % save STA jpeg
            saveas(H1,fName);
            saveas(H1,fNameJ);
            close(H1);
            
            % Burst selectivity
            % PhasicB
            H3 = figure;
            set(gcf, 'Renderer', 'painters');
            hold on;
            errorshade(time, mn_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iP})...
                .burst1-min(mn_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iP}).burst1),...
                se_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iP}).burst1,...
                'LineColor', 'm', 'ShadeColor', 'm', 'LineWidth', LWidth);
            errorshade(time, mn_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iP})...
                .single-min(mn_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iP}).single),...
                se_phasicB.(tTypes{iT}).(fType{iF}).(cTypes{iP}).single,...
                'LineColor', [0 0 0], 'ShadeColor', [0 0 0], 'LineWidth', LWidth);
            title(['Selectivity ' tTypes{iT} ' ' cTypes{iP} ' '...
                fType{iF} ' phasicB'])
            setmyplot_tamas;
            legend({'Bursts', 'Single spikes'}, 'FontSize', 12, 'Location',...
                'northeast')
            legend('boxoff');
            xticks(x_ticks);
            xticklabels(xtick_labels);
            xlim(plotTime);
            xlabel('Time from punishment (ms)');
            ylabel('Firing rate (Hz)');
            legend({['n=' num2str(sum(phasiC.(tTypes{iT}).(cTypes{iP})))]},...
                'FontSize', 12, 'Location','northeast');
            axis square;
            fName = [psthDir 'burstSelectivity_' tTypes{iT} '_' cTypes{iP} '_'...
                fType{iF} '_' outText '_phasicB.fig'];   % save PSTH figure
            fNameJ = [psthDir 'burstSelectivity_' tTypes{iT} '_' cTypes{iP} '_'...
                fType{iF} '_' outText '_phasicB.jpeg'];   % save STA jpeg
            saveas(H3,fName);
            saveas(H3,fNameJ);
            close(H3);
            
            
            % Burst selectivity
            % PoissonL
            H32 = figure;
            set(gcf, 'Renderer', 'painters');
            hold on;
            errorshade(time, mn_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iP})...
                .burst1-min(mn_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iP}).burst1),...
                se_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iP}).burst1,...
                'LineColor', 'm', 'ShadeColor', 'm', 'LineWidth', LWidth);
            errorshade(time, mn_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iP})...
                .single-min(mn_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iP}).single),...
                se_poissonL.(tTypes{iT}).(fType{iF}).(cTypes{iP}).single,...
                'LineColor', [0 0 0], 'ShadeColor', [0 0 0], 'LineWidth', LWidth);
            title(['Selectivity ' tTypes{iT} ' ' cTypes{iP} ' '...
                fType{iF} ' poissonL'])
            setmyplot_tamas;
            legend({['n=' num2str(sum(poissoN.(tTypes{iT}).(cTypes{iP})))]},...
                'FontSize', 12, 'Location','northeast');
            legend('boxoff');
            xticks(x_ticks);
            xticklabels(xtick_labels);
            xlim(plotTime);
            xlabel('Time from punishment (ms)');
            ylabel('Firing rate (Hz)');
            axis square;
            fName = [psthDir 'burstSelectivity_' tTypes{iT} '_' cTypes{iP} '_'...
                fType{iF} '_' outText '_poissonL.fig'];   % save PSTH figure
            fNameJ = [psthDir 'burstSelectivity_' tTypes{iT} '_' cTypes{iP} '_'...
                fType{iF} '_' outText '_poissonL.jpeg'];   % save STA jpeg
            saveas(H32,fName);
            saveas(H32,fNameJ);
            close(H32);
            
            
            % Boxplot for burstIndexes
            if ~isempty(burstInd)
                % PhasicB
                [H4, Wp] = boxstat(burstInd,singleInd,'Burst1','Single AP',...
                    0.05,'paired');
                title(['Burst1 vs Single ' cTypes{iP} ' ' fType{iF} ' ' outText])
                setmyplot_tamas;
                storeP = round(Wp,4);
                close(H4);
                
                % Barplot for burstIndexes median
                H5 = figure;
                hold on;
                line([1, 2], [burstInd, singleInd], 'Color', [0.6 0.6 0.6],...
                    'LineWidth',3)
                bar(1, median(burstInd), 'FaceColor','none','EdgeColor','m',...
                    'LineWidth',3)
                bar(2, median(singleInd), 'FaceColor','none','EdgeColor',...
                    [0 0 0],'LineWidth',3)
                xticks([0.9 2.1]);
                xticklabels({'Bursts', 'Single spikes'});
                ylabel('Selectivity Index')
                title(['Burst1 vs Single ' cTypes{iP} ' ' fType{iF} ' phasicB'])
                axis square;
                setmyplot_tamas;
                x_lim = xlim;
                y_lim = ylim;
                text((x_lim(2)*0.7), y_lim(2)*0.95, ['p=' num2str(storeP)]);
                setmyplot_tamas;
                fName = [psthDir 'burstIndex_' tTypes{iT} '_' cTypes{iP} '_'...
                    fType{iF} '_' outText '_median_phasicB.fig'];   % save PSTH figure
                fNameJ = [psthDir 'burstIndex_' tTypes{iT} '_' cTypes{iP} '_'...
                    fType{iF} '_' outText '_median_phasicB.jpeg'];   % save STA jpeg
                saveas(H5,fName);
                saveas(H5,fNameJ);
                close(H5);
                
                % PoissonL
                [H42, Wp2] = boxstat(burstIndP,singleIndP,'Burst1','Single AP',...
                    0.05,'paired');
                title(['Burst1 vs Single ' cTypes{iP} ' ' fType{iF} ' ' outText])
                setmyplot_tamas;
                storeP2 = round(Wp2,4);
                close(H42);
                
                % Barplot for burstIndexes median
                H52 = figure;
                hold on;
                line([1, 2], [burstIndP, singleIndP], 'Color', [0.6 0.6 0.6],...
                    'LineWidth',3)
                bar(1, median(burstIndP), 'FaceColor','none','EdgeColor','m',...
                    'LineWidth',3)
                bar(2, median(singleIndP), 'FaceColor','none','EdgeColor',...
                    [0 0 0],'LineWidth',3)
                xticks([0.9 2.1]);
                xticklabels({'Bursts', 'Single spikes'});
                ylabel('Selectivity Index')
                title(['Burst1 vs Single ' cTypes{iP} ' ' fType{iF} ' poissonL'])
                axis square;
                setmyplot_tamas;
                x_lim = xlim;
                y_lim = ylim;
                text((x_lim(2)*0.7), y_lim(2)*0.95, ['p=' num2str(storeP2)]);
                setmyplot_tamas;
                fName = [psthDir 'burstIndex_' tTypes{iT} '_' cTypes{iP} '_'...
                    fType{iF} '_' outText '_median_poissonL.fig'];   % save PSTH figure
                fNameJ = [psthDir 'burstIndex_' tTypes{iT} '_' cTypes{iP} '_'...
                    fType{iF} '_' outText '_median_poissonL.jpeg'];   % save STA jpeg
                saveas(H52,fName);
                saveas(H52,fNameJ);
                close(H52);
            end
            
            % FIGURE3 PANEL C
            % Individual raster+psth plots
            if doraster
                if strcmp((cTypes{iP}), 'allC')
                    for aI = 1:length(allChAT.(tTypes{iT}))
                        exampleCells = {'n046_130108a_8.1', 'n029_120215a_3.4',...
                            'n046_121230a_1.2'};
                        
                        if any(contains(exampleCells, allChAT.(tTypes{iT}){aI}))
                            %Store raster data
                            rasterNone = normData.(tTypes{iT}).(fType{iF})...
                                .none.raster{aI};
                            rasterBurst1 = normData.(tTypes{iT}).(fType{iF})...
                                .burst1.raster{aI};
                            cellidt = regexprep(allChAT.(tTypes{iT}){aI},'\.','_');
                            noneData = normData.(tTypes{iT}).(fType{iF}).none...
                                .(psthFilt);
                            
                            burInx = bInx.(tTypes{iT}).(fType{iF}).(cTypes{iP});
                            sinInx = sInx.(tTypes{iT}).(fType{iF}).(cTypes{iP});
                            
                            % Psth
                            psthDat = normData.(tTypes{iT}).(fType{iF}).none...
                                .spsth(aI,plotWindow);
                            
                            % Creates plot matrix (1 - Single spikes, 2 - Burst1)
                            plotImage = rasterNone(:,plotWindow);
                            middle = round(size(plotImage,2)/2);
                            
                            cbSwitcher(allChAT.(tTypes{iT}){aI});
                            cbName = whichcb;
                            
                            switch groupID(aI)
                                case 1
                                    group = 'PhasicB';
                                case 2
                                    group = 'PoissonL';
                                case 3
                                    group = 'Tonic';
                            end
                            
                            fnmS1 = [resdir 'raster' cb '\individuals' PATH '\'...
                                group '_' cellidt...
                                '_' tTypes{iT} '_' (fType{iF}) '_' psthFilt...
                                '_marged.fig'];
                            fnmSE1 = [resdir 'raster' cb '\individuals' PATH...
                                '\'  group '_' cellidt...
                                '_' tTypes{iT} '_' (fType{iF}) '_' psthFilt...
                                '_merged.jpeg'];
                            
                            
                            H8 = figure;
                            hold on;
                            yyaxis left
                            imagesc(plotImage);
                            map = [1 1 1; cCode(groupID(aI),:)];
                            colormap(map);
                            xlim([middle middle]+plotTime);
                            yticks(ylim)
                            xticks(x_ticks);
                            xticklabels(xtick_labels);
                            xticklabels({'0' '500'});
                            ylabel('Trials#');
                            xlabel('Time from punishment (ms)')
                            title([group ' example'])
                            y_lim2 = ylim;
                            ylim(y_lim2)
                            
                            yyaxis right
                            plot(psthDat,'Color', [0.6 0.6 0.6], 'LineWidth',...
                                LWidth);
                            xlim([middle middle]+plotTime);
                            line([middle; middle], ylim, 'Color', [0 0 1 0.5],...
                                'LineWidth', LWidth);
                            y_lim = ylim;
                            y_ticks = [0 y_lim(2)/3 y_lim(2)/3*2 y_lim(2)];
                            yticks(y_ticks);
                            yticklabels({'0', ' ', ' ', num2str(y_lim(2))});
                            ylabel('Firing rate (Hz)');
                            axis square;
                            setmyplot_tamas;
                            ylim(y_lim);
                            ax1 = gca;
                            ax1.YAxis(2).Color = [0.6 0.6 0.6];
                            ax1.YAxis(1).Color = [0 0 0];
                            saveas(H8, fnmSE1);
                            saveas(H8, fnmS1);
                            close(H8);
                        end
                    end
                end
            end
        end
    end
end
keyboard;

pLimits = [-500 500];
xtick_labels = {'-500', '0', '500'};
x_ticks = [-500 0 500];
time = pLimits(1):1:pLimits(2); % time vector for plotting
plotWindow = (centerP+time(1)):(centerP+time(end));
plotTime = pLimits;



% REVIEWER FIGURE6
% Figure for allTtypes
legText = {'Miss','CR','Hit', 'FA'};
alignN = 'cue';

% Burst-SB
H1 = figure;
set(gcf, 'Renderer', 'painters');
hold on;
errorshade(time, mn_phasicB.MISS.Miss.allC.none-min(mn_phasicB.MISS.Miss...
    .allC.none),...
    se_phasicB.MISS.Miss.allC.none, 'LineColor', [0 0 0], 'ShadeColor',...
    [0 0 0], 'LineWidth', LWidth);
errorshade(time, mn_phasicB.CR.CorrectRejection.allC.none-min(mn_phasicB.CR...
    .CorrectRejection.allC.none),...
    se_phasicB.CR.CorrectRejection.allC.none, 'LineColor', [0 0 0.75],...
    'ShadeColor', [0 0 0.75], 'LineWidth', LWidth);
errorshade(time, mn_phasicB.TONE.HIT.allC.none-min(mn_phasicB.TONE.HIT...
    .allC.none),...
    se_phasicB.TONE.HIT.allC.none, 'LineColor', cCode(3,:), 'ShadeColor',...
    cCode(3,:), 'LineWidth', LWidth);
errorshade(time, mn_phasicB.TONE.FA.allC.none-min(mn_phasicB.TONE.FA.allC...
    .none),...
    se_phasicB.TONE.FA.allC.none, 'LineColor', cCode(1,:), 'ShadeColor',...
    cCode(1,:), 'LineWidth', LWidth);
legend(legText, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
title('Burst-SB')
setmyplot_tamas;
xlim(plotTime);
ylim([-2 20])
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from ' alignN ' (ms)']);
ylabel('Firing rate (Hz)');
axis square;
fName = [psthDir '\allPlot\phasicB.fig'];   % save PSTH figure
fNameJ = [psthDir '\allPlot\phasicB.jpeg'];   % save STA jpeg
saveas(H1,fName);
saveas(H1,fNameJ);
close(H1);

% Burst-PL
H1 = figure;
set(gcf, 'Renderer', 'painters');
hold on;
errorshade(time, mn_poissonL.MISS.Miss.allC.none-min(mn_poissonL.MISS...
    .Miss.allC.none),...
    se_phasicB.MISS.Miss.allC.none, 'LineColor', [0 0 0], 'ShadeColor',...
    [0 0 0], 'LineWidth', LWidth);
errorshade(time, mn_poissonL.CR.CorrectRejection.allC.none-min(mn_poissonL...
    .CR.CorrectRejection.allC.none),...
    se_phasicB.CR.CorrectRejection.allC.none, 'LineColor', [0 0 0.75],...
    'ShadeColor', [0 0 0.75], 'LineWidth', LWidth);
errorshade(time, mn_poissonL.TONE.HIT.allC.none-min(mn_poissonL.TONE.HIT...
    .allC.none),...
    se_phasicB.TONE.HIT.allC.none, 'LineColor', cCode(3,:), 'ShadeColor',...
    cCode(3,:), 'LineWidth', LWidth);
errorshade(time, mn_poissonL.TONE.FA.allC.none-min(mn_poissonL.TONE.FA...
    .allC.none),...
    se_phasicB.TONE.FA.allC.none, 'LineColor', cCode(1,:), 'ShadeColor',...
    cCode(1,:), 'LineWidth', LWidth);
legend(legText, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
title('Burst-PL')
setmyplot_tamas;
xlim(plotTime);
ylim([-2 20])
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from ' alignN ' (ms)']);
ylabel('Firing rate (Hz)');
axis square;
fName = [psthDir '\allPlot\poissonL.fig'];   % save PSTH figure
fNameJ = [psthDir '\allPlot\poissonL.jpeg'];   % save STA jpeg
saveas(H1,fName);
saveas(H1,fNameJ);
close(H1);

% Regular
H1 = figure;
set(gcf, 'Renderer', 'painters');
hold on;
errorshade(time, mn_tonic.MISS.Miss.allC.none-min(mn_tonic.MISS.Miss.allC.none),...
    se_phasicB.MISS.Miss.allC.none, 'LineColor', [0 0 0], 'ShadeColor',...
    [0 0 0], 'LineWidth', LWidth);
errorshade(time, mn_tonic.CR.CorrectRejection.allC.none-min(mn_tonic.CR...
    .CorrectRejection.allC.none),...
    se_phasicB.CR.CorrectRejection.allC.none, 'LineColor', [0 0 0.75],...
    'ShadeColor', [0 0 0.75], 'LineWidth', LWidth);
errorshade(time, mn_tonic.TONE.HIT.allC.none-min(mn_tonic.TONE.HIT.allC...
    .none),...
    se_phasicB.TONE.HIT.allC.none, 'LineColor', cCode(3,:), 'ShadeColor',...
    cCode(3,:), 'LineWidth', LWidth);
errorshade(time, mn_tonic.TONE.FA.allC.none-min(mn_tonic.TONE.FA.allC.none),...
    se_phasicB.TONE.FA.allC.none, 'LineColor', cCode(1,:), 'ShadeColor',...
    cCode(1,:), 'LineWidth', LWidth);
legend(legText, 'FontSize', 12, 'Location','northeast')
legend('boxoff');
title('Regular')
setmyplot_tamas;
xlim(plotTime);
ylim([-2 20])
xticks(x_ticks);
xticklabels(xtick_labels);
xlabel(['Time from ' alignN ' (ms)']);
ylabel('Firing rate (Hz)');
axis square;
fName = [psthDir '\allPlot\tonic.fig'];   % save PSTH figure
fNameJ = [psthDir '\allPlot\tonic.jpeg'];   % save STA jpeg
saveas(H1,fName);
saveas(H1,fNameJ);
close(H1);







% -------------------------------------------------------------------------
function [currData, indexes] = clearEmpty(currData)

indexes = [];
filters = fieldnames(currData);

for i=1:length(filters)
    fieldCheck = fieldnames(currData.(filters{i}));
    if contains(fieldCheck, 'NotAdaptive')
        currData.(filters{i}) = currData.(filters{i}).NotAdaptive;
    end
    indexes = find(cellfun(@isempty,currData.(filters{i}).stats));
    currData.(filters{i}).stats(indexes)=[];
    currData.(filters{i}).psth(indexes,:) = [];
    currData.(filters{i}).spsth(indexes,:) = [];
    currData.(filters{i}).spsth_se(indexes,:) = [];
    currData.(filters{i}).raster(indexes) = [];
    if isfield(currData.(filters{i}), 'cellid')
        currData.(filters{i}).cellid(indexes) = [];
    end
end

% -------------------------------------------------------------------------
function currDataNorm = select_normalization(currData)

fld = fieldnames(currData);
NumFields = length(fld);
currDataNorm = currData;
for iF = 1:NumFields
    nfcn = @normfun0;  % no normalization
    if ismember(iF,[2 3])
        nfcn = @normfun3;   % divide by baseline (relative increase)
    end
    fld2 = fieldnames(currData.(fld{iF}));
    for iF2 = 1:3
        lpsth = currData.(fld{iF}).(fld2{iF2});
        if ~ismember(iF,[2 3])
            nlpsth = feval(nfcn,lpsth);  % call normalization function
        else
            % normalization: divide by baseline (relative increase) for both 
            % burst1 and single for proper comparison
            nlpsth = feval(nfcn,lpsth);  
        end
        currDataNorm.(fld{iF}).(fld2{iF2}) = nlpsth;
    end
end

% -------------------------------------------------------------------------
function N = normfun0(M)
% no normalization - note that Z-scoring will also abolish all vs burs1+single differences
N = M;  

% -------------------------------------------------------------------------
function [N, mn, sd] = normfun1(M)

s = size(M);
mn = repmat(mean(M,2),1,s(2));
sd = repmat(std(M,[],2),1,s(2));
N = (M - mn) ./ sd;   % Z-score
N = nan2zero(N);

% -------------------------------------------------------------------------
function [N, mn, sd] = normfun2(M,mn,sd)

N = (M - mn) ./ sd;   % Z-score with input mean and SD
N = nan2zero(N);

% -------------------------------------------------------------------------
function [N, mn] = normfun3(M)

s = size(M);
baseWindow = 600:750;
mn = repmat(mean(M(baseWindow),2),1,s(2));
N = M ./ (1 + mn);   % devide by baseline mean
N = nan2zero(N);

% -------------------------------------------------------------------------
function [allData] = loadSeparate(actBurst)

if actBurst
    nb = load('D:\_MATLAB_DATA\NB\PSTH_newdataNB\EventData_review_removed_22-Jan-2020.mat');
    hdb = load('D:\_MATLAB_DATA\HDB\PSTH_newdataHDB\EventData_review_removed_22-Jan-2020.mat');
    p_hdb = load('D:\_MATLAB_DATA\PannaHDB\PSTH_newdataPannaHDB\EventData_review_removed_22-Jan-2020.mat');
else
    nb = load('D:\_MATLAB_DATA\NB\PSTH_newdataNB\EventData_review__21-Jan-2020.mat');
    hdb = load('D:\_MATLAB_DATA\HDB\PSTH_newdataHDB\EventData_review__21-Jan-2020.mat');
    p_hdb = load('D:\_MATLAB_DATA\PannaHDB\PSTH_newdataPannaHDB\EventData_review__21-Jan-2020.mat');
end
bFilt = {'none', 'burst1', 'single', 'burstall', 'burst1Single'};
Type = fieldnames(nb.EventData);

% Extract data for all cb-s, trialtypes, and burstfilters
for iT = 1:length(Type)
    fType = fieldnames(nb.EventData.(Type{iT}));
for iF = 1:length(fType)
    if strcmp(fType,'Miss')
        nb.EventData.(Type{iT}).(fType{iF}) = nb.EventData.(Type{iT}).(fType{iF}).fa;
        hdb.EventData.(Type{iT}).(fType{iF}) = hdb.EventData.(Type{iT}).(fType{iF}).fa;
        p_hdb.EventData.(Type{iT}).(fType{iF}) = p_hdb.EventData.(Type{iT}).(fType{iF}).fa;
    elseif strcmp(fType, 'CorrectRejection')
        nb.EventData.(Type{iT}).(fType{iF}) = nb.EventData.(Type{iT}).(fType{iF}).hit;
        hdb.EventData.(Type{iT}).(fType{iF}) = hdb.EventData.(Type{iT}).(fType{iF}).hit;
        p_hdb.EventData.(Type{iT}).(fType{iF}) = p_hdb.EventData.(Type{iT}).(fType{iF}).hit;
    end

    for iB = 1:length(bFilt)
        % spsth
        allData.(Type{iT}).(fType{iF}).(bFilt{iB}).spsth = [nb.EventData.(Type{iT}).(fType{iF})...
            .NotAdaptive.(bFilt{iB}).spsth;...
            hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).spsth;...
            p_hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).spsth];
        % psth
        allData.(Type{iT}).(fType{iF}).(bFilt{iB}).psth = [nb.EventData.(Type{iT}).(fType{iF})...
            .NotAdaptive.(bFilt{iB}).psth;...
            hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).psth;...
            p_hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).psth];
        % spsth_se
        allData.(Type{iT}).(fType{iF}).(bFilt{iB}).spsth_se = [nb.EventData.(Type{iT}).(fType{iF})...
            .NotAdaptive.(bFilt{iB}).spsth_se;...
            hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).spsth_se;...
            p_hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).spsth_se];
        % stats
        allData.(Type{iT}).(fType{iF}).(bFilt{iB}).stats = [nb.EventData.(Type{iT}).(fType{iF})...
            .NotAdaptive.(bFilt{iB}).stats';...
            hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).stats';...
            p_hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).stats'];
        % cellid
        allData.(Type{iT}).(fType{iF}).(bFilt{iB}).cellid = [nb.EventData.(Type{iT}).(fType{iF})...
            .NotAdaptive.(bFilt{iB}).cellid';...
            hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).cellid';...
            p_hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).cellid'];
        % raster
        allData.(Type{iT}).(fType{iF}).(bFilt{iB}).raster = [nb.EventData.(Type{iT}).(fType{iF})...
            .NotAdaptive.(bFilt{iB}).raster';...
            hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).raster';...
            p_hdb.EventData.(Type{iT}).(fType{iF}).NotAdaptive.(bFilt{iB}).raster'];
    end
end
end
