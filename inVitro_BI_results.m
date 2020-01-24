function inVitro_BI_results
% INVITRO_BI_RESULTS
% Plotting function for resampled in vivo BurstIndex data
% It compares the original and the downsampled BurstIndex changes
% and plots their distribution. It also plots the difference between
% the in vivo and in vitro segment lengths.
%
% Check also the parent function where the data comes from:
% See also INVITRO_BI_CHECK

global cCode;

% Load in vitro and in vivo acg data
fNameVitro = which('ACG_matrices_POOLED_inVitro.mat');
fNameVivo = which('ACG_matrices_POOLED.mat');
inVitro = load(fNameVitro);
inVivo = load(fNameVivo);
[resDir,~,~] = fileparts(fNameVitro);

% Match cellids and structure data
for iC=1:length(inVitro.cellids)
    
    currCell = inVitro.cellids{iC};
    cellPos = find(contains(inVivo.cellids, currCell));
    if ~isempty(cellPos)
        
        % Load segment length difference for current cell
        diffDat = load(which([regexprep(currCell,'\.','_')...
            '_lengthDiff.mat']));
        biData.minDiff(iC) = diffDat.minDiff;
        
        % Exclude outliers for minDiff larger than 100s
        errInx(iC) = diffDat.minDiff < 100;
        
        % Structure data inVitro
        biData.inVitro.BurstIndex(iC) = inVitro.BurstIndex(iC);
        biData.inVitro.groupID(iC) = inVitro.groupID(iC);
        biData.inVitro.cellids(iC) = inVitro.cellids(iC);
        
        % Structure data inVivo
        biData.inVivo.BurstIndex(iC) = inVivo.BurstIndex(cellPos);
        biData.inVivo.groupID(iC) = inVivo.groupID(cellPos);
        biData.inVivo.cellids(iC) = inVitro.cellids(iC);
    end
    
end

%% Plot BI differences

% Boxplot for burstIndexes

origiData = biData.inVivo.BurstIndex;
resampData = biData.inVitro.BurstIndex;

% Cells which change their group identity
groupChange = biData.inVivo.groupID ~= biData.inVitro.groupID;

upInx = origiData<resampData;
downInx = ~upInx;
upGroupID = biData.inVivo.groupID(upInx);
downGroupID = biData.inVivo.groupID(downInx);

filts.bfcnBurst = biData.inVivo.groupID ~= 3;
filts.bfcnReg = ~filts.bfcnBurst;

plotFilts = {'bfcnBurst', 'bfcnReg'};

errInx_Ori = errInx;
for iF = 1:numel(plotFilts)
    
    if strcmp(plotFilts{iF}, 'bfcnBurst')
       plotClr = cCode(1,:);
    else
        plotClr = cCode(3,:);
    end
    
    
    currFilt = filts.(plotFilts{iF});
    % BI plots
    [H1, Wp] = boxstat(origiData(currFilt),resampData(currFilt),'Original','Downsampled',0.05,'paired');
    title([plotFilts{iF} ' Original BI vs Downsampled'])
    setmyplot_tamas;
    storeP = round(Wp,4);
    fName = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI.fig'];
    fNameJ = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI.jpeg'];
    saveas(H1,fName);
    saveas(H1,fNameJ);
    close(H1);
    
    % Barplot for burstIndexes median
    H2 = figure;
    hold on;
    line([1, 2], [origiData(currFilt)', resampData(currFilt)'], 'Color', [0.6 0.6 0.6], 'LineWidth',3)
    bar(1, median(origiData(currFilt)), 'FaceColor','none','EdgeColor',plotClr,'LineWidth',3)
    bar(2, median(resampData(currFilt)), 'FaceColor','none','EdgeColor',[0 0 0],'LineWidth',3)
    xticks([0.9 2.1]);
    xticklabels({'Original', 'Downsampled'});
    ylabel('BurstIndex')
    axis square;
    title([plotFilts{iF} ' Original BI vs Downsampled'])
    setmyplot_tamas;
    x_lim = xlim;
    y_lim = ylim;
    if storeP<0.01
        text((x_lim(2)*0.7), y_lim(2)*0.95, 'p<0.01');
    else
        text((x_lim(2)*0.7), y_lim(2)*0.95, ['p=' num2str(storeP)]);
    end
    setmyplot_tamas;
    fName = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI_dist.fig'];
    fNameJ = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI_dist.jpeg'];
    saveas(H2,fName);
    saveas(H2,fNameJ);
    close(H2);
    
    % Barplot for minDiff median
    
    errInx = currFilt & errInx_Ori;
    
    x1 = rand(length(biData.minDiff(errInx)),1)+0.5;
    
    H3 = figure;
    hold on;
    bar(1, median(biData.minDiff(errInx)), 'FaceColor','none','EdgeColor', [0.6 0.6 0.6], 'LineWidth',3)
    plot(x1, biData.minDiff(errInx), 'o', 'Color', [0 0.6 0], 'LineWidth',3)
    xticks([1]);
    xticklabels('Length difference (s)');
    ylabel('Count')
    axis square;
    title([plotFilts{iF} ' Length Diff'])
    setmyplot_tamas;
    x_lim = xlim;
    y_lim = ylim;
    text((x_lim(2)*0.7), y_lim(2)*0.95,...
        ['median: ' num2str(median(biData.minDiff(errInx)))]);
    setmyplot_tamas;
    fName = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI_lengthDiff.fig'];
    fNameJ = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI_lengthDiff.jpeg'];
    saveas(H3,fName);
    saveas(H3,fNameJ);
    close(H3);
    
    % minDiff histogram
    H4 = figure;
    hold on;
    hist(biData.minDiff(errInx), 50)
    xlabel('Length difference (s)');
    ylabel('Count')
    axis square;
    title([plotFilts{iF} ' Length Diff'])
    setmyplot_tamas;
    x_lim = xlim;
    y_lim = ylim;
    text((x_lim(2)*0.7), y_lim(2)*0.95,...
        ['median: ' num2str(median(biData.minDiff(errInx)))]);
    setmyplot_tamas;
    fName = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI_lengthDiff_hist.fig'];
    fNameJ = [resDir filesep plotFilts{iF} '_original_vs_downsampled_BI_lengthDiff_hist.jpeg'];
    saveas(H4,fName);
    saveas(H4,fNameJ);
    close(H4);
end

