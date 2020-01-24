function ccg_NBsync(cellids,varargin)
%CCG   Cross-correlation.
%   CCG calculates cross-correlations. Window size is set to +-1000 ms with
%   a 1 ms resolution. Maximum 50000 spikes are included to avoid memory
%   problems. CCG is not calculated if one of the cells has less than 100
%   spikes. Segments are filtered with FINDSEGS3. Minimal shift for
%   shuffled CCGs is set to 1100 ms. For details on the algorithm, see
%   SOMCCG_CONF_FILTER.
%
%   CCG(I) calls CCG for pairs of cells defined by the cell ID list
%   (or index set to CELLIDLIST, see CellBase documentation) I. By default,
%   non-tetrode pairs within the given list are selected for analysis.
%   Optional input arguments (parameter-value pairs with default values):
%       'issave', false - controls saving behavior; plots and
%           cross-correlation matrices with confidence intervals are saved
%           only if 'issave' is set to true
%       'whichcells', 'nontetrodepairs' - method of pair selection;
%           'nontetrodepairs' selects cells from other tetrodes,
%           'tetrodepairs' selects cells from the same tetrode and
%           'allpairs' selects all cells from the session
%       'include', 'list' - by default, only pairs for which both cells are
%           included in I are analyzed; if 'include' is set to 'cellbase',
%           all cells in CellBase that are paired with the ones in I
%           according to 'whichcells' are analyzed
%
%   See also ACG, SOMCCG_CONF_FILTER and XCORR_WRAND_FILTER.

%   Balazs Hangya, Cold Spring Harbor Laboratory
%   1 Bungtown Road, Cold Spring Harbor
%   balazs.cshl@gmail.com
%   14-June-2013

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
addParameter(prs,'issave',false,@islogical)   % control saving behavior
addParameter(prs,'whichcells','nontetrodepairs',...
    @(s)ismember(s,{'nonterodepairs','tetrodepairs','allpairs'}))   % which cells to include
addParameter(prs,'include','list',...
    @(s)ismember(s,{'list','cellbase'}))   % cell pair selection behavior: only from input list or full CellBase
addParameter(prs,'burstfilter','none',@(s)ischar(s)&ismember(s,{'none',...
    'burst1', 'burstall', 'single', 'burst1Single'}))   % burstfilter selector
addParameter(prs,'runtime','first',@(s)ischar(s)&ismember(s,{'first',...
    'repeated', 'burst1Corr'}))   % select between first run or repeated when the pairs are already found
addParameter(prs,'cells','NB',@(s)ischar(s)&ismember(s,{'NB',...
    'HDB', 'CPu', 'ACb'}))   % selecting celltype
addParameter(prs,'sessiontype','behavior',@(s)ischar(s)&ismember(s,{'behavior',...
    'ITI', 'sleep', 'freely moving', 'quiet wakefulness'}))   % session type selector
parse(prs,cellids,varargin{:})
g = prs.Results;

% Pass the control to the user in case of error
dbstop if error
cellbase = whichcb;
cellbase = 'POOLED';
% Directories
global RESDIR;
global PATH;
fs = filesep;
resdir = [RESDIR cellbase fs 'ccg' cellbase fs regexprep(g.sessiontype,' ','_') fs g.burstfilter cellbase PATH fs];  % results directory
fnmm = ['CCG_matrices_' regexprep(g.sessiontype,' ','_') '_' g.burstfilter '_' cellbase '.mat'];   % filename for saving the result matrices
acgDir = [RESDIR cellbase fs 'acg' cellbase PATH fs];
acgPath = [acgDir 'ACG_matrices_' cellbase '.mat'];

% Include only long enough segments
longsegments = false;  % control whether to use this option

% Determine time window
sr = 1000;      % sampling rate
wn = 250 * sr / 1000;    % 2 * 250 ms window

% Cell pairs
PairOfCells = cell(0,2);
numCells = length(cellids);  % number of cells
switch g.runtime
    case 'first'
        for iC = 1:numCells
            cellid = cellids{iC};
            [animalID, sessionID, tetrode1, unit1] = cellid2tags(cellid);
            switch g.whichcells
                case 'nontetrodepairs'
                    [nm, ps] = nontetrodepairs(cellid);   % cells on other tetrodes
                case 'tetrodepairs'
                    [nm, ps] = tetrodepairs(cellid);   % cells on the same tetrode
                    ps = setdiff(ps,cellid);
                    nm = length(ps);
                case 'allpairs'
                    ps = findcell('rat',animalID,'session',sessionID);   % all concurrently recorded cells
                    ps = setdiff(ps,cellid);
                    nm = length(ps);
            end
            if isequal(g.include,'list')
                ps = intersect(cellids,ps);   % include only those that are in the input list
                nm = length(ps);
            end
            for k = 1:nm
                [j1, j2, tetrode2, unit2] = cellid2tags(ps(k));
                if (tetrode2*10+unit2) > (tetrode1*10+unit1)   % prevent duplicates
                    PairOfCells(end+1,1:2) = {cellid ps{k}}; %#ok<AGROW>
                end
            end
        end
        numPairs = size(PairOfCells,1);
    case 'repeated'
        load([RESDIR cellbase '\ccg' cellbase '\behavior\none' cellbase...
            PATH '\CCG_matrices_behavior_none_' cellbase '.mat']);   % Load CCG matrix
        numPairs = size(PairOfCells,1);
        disp('Pair of cells already calculated. Load old CCG matrix')
    case 'burst1Corr'
        PairOfCells = [cellids', cellids'];
        numPairs = size(PairOfCells,1);
end

% Cell pair loop for CCG
wb = waitbar(0,'Please wait...','Name','Running CCG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [0 50000];   % include max 50000 spikes; calculate only if min 100 spikes
[CCR, LCCR, UCCR] = deal(zeros(numPairs,2*wn));
[MeanH0 SDH0] = deal(nan(numPairs,1));
SDH0 = nan(numPairs,1);
SegmentLength = nan(numPairs,1);
chatIndex = [];
delInx = [];
repInx = [];
for iP = 1:numPairs   % loop through pairs of cells
    cell1 = PairOfCells{iP,1};
    cell2 = PairOfCells{iP,2};
    try
        % Set cellbase
        cbSwitcher(cell1);
        
        switch g.sessiontype
            case 'behavior'
                tseg = findSegs3(cell1,'segfilter','stimfb_excl_nb',...
                    'light_activation_duration',[-5 5],'margins',[0 0]);  % find time segments
            case 'ITI'
                tseg = findSegs3(cell1,'segfilter','ITI');  % find time segments
            otherwise
                [xnum, xstr, xcell] = xlsread([RESDIR cellbase '\sta' cellbase '\' 'sleep_segments2_refilled.xlsx'], 1);
                sleepsessinx = find(strcmp(g.sessiontype,xcell(1,:)));
                if isempty(sleepsessinx)
                    disp('The current cell is not recorded in a freely moving session')
                end
                [rowinx, ~] = find(strcmp(cell1,xcell(:,2)),1);   % find cell in table
                SegNum = length(rowinx);   % number of segments
                if SegNum==0
                    disp('The current cell is not recorded in a freely moving session')
                end
                starttime = xcell{rowinx, sleepsessinx};
                endtime = xcell{rowinx, sleepsessinx+1};
                if isnan(endtime)
                    disp([g.sessiontype ' data is missing for the current cell'])
                end
                %         tseg = findSegs3(cellid,'segfilter','prestim');  % find time segments
                SE = loadcb(cell1,'StimEvents');
                tseg = [SE.PulseOn - 1; SE.PulseOn + 5];
                [rowN, columnN] = find(tseg>endtime,1);
                if isempty(rowN) % All pulses are in the defined segment (prestim)
                    FMseg = [starttime tseg(2,end); tseg(1,1) endtime];
                elseif rowN==1 && columnN==1 % Freely moving recording instead of a pretrial tagging
                    FMseg = [xcell{rowinx, sleepsessinx}; endtime];
                else
                    if rowN==1
                        FMseg = [starttime tseg(2,columnN-1);  tseg(1,1), endtime];
                    elseif rowN==2
                        FMseg = [starttime tseg(1,columnN);  tseg(1,1) endtime];
                    else
                        keyboard;
                    end
                end
                tseg = FMseg;
        end
        %         tseg = findSegs3(cell1,'segfilter','prestim3');  % find time segments
        %         tseg = findSegs3(cell1,'segfilter','fb_incl_nb',...
        %             'feedback_duration',[-0.5 0.5],'margins',[0 0],'min_int',0);  % find time segments
        %         tseg = findSegs3(cell1,'segfilter','cue_incl_nb',...
        %             'feedback_duration',[-1.5 0],'margins',[0 0],'min_int',0);  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if longsegments   % later to be implemented as an input option
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < seglim   % use the longest segment if it's longer than the threshold
                continue
            end
        end
        SegmentLength(iP) = sum(ltseg);  % cumulative length of the segments
        if strcmp('burst1Corr', g.runtime)
            [ncc1,seltsind,selisi] = extractSegSpikes_NBsync(cell1,tseg, 'burstfilter', 'burst1');   % find spikes in the time segments
            [ncc2,seltsind,selisi] = extractSegSpikes_NBsync(cell2,tseg, 'burstfilter', 'none');
        else
            [ncc1,seltsind,selisi] = extractSegSpikes_NBsync(cell1,tseg, 'burstfilter', g.burstfilter);   % find spikes in the time segments
            [ncc2,seltsind,selisi] = extractSegSpikes_NBsync(cell2,tseg, 'burstfilter', g.burstfilter);
        end
    catch ME
        disp(ME.message)
        disp('Could not extract the segment.');
    end
    
    %     ncc1 = loadcb(cell1,'SPIKES');   % use all spikes
    %     ncc2 = loadcb(cell2,'SPIKES');
    
    if length(ncc1) > limit_spikes(2);      % crop if too long to avoid out of memory
        ncc1 = ncc1(1:limit_spikes(2));
    end
    if length(ncc2) > limit_spikes(2);
        ncc2 = ncc2(1:limit_spikes(2));
    end
    
    if length(ncc1) > limit_spikes(1) && length(ncc2) > limit_spikes(1)     % minimum 100 spikes
        [H1, ccr, lwr, upr, rccg] = somccg_conf_filter_NBsync(ncc1,ncc2,wn,1100);    % 1->2
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        %         text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.9,regexprep(cell1,'_',' '))
        %         text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,regexprep(cell2,'_',' '))
        try
            [group1{iP}, ~] = cellid2TPgroup(cell1, 'acgPath', acgPath);
        catch
            group1{iP} = {'NaN'};
            delInx = [delInx, iP];
        end
        try
            [group2{iP}, ~] = cellid2TPgroup(cell2, 'acgPath', acgPath);
        catch
            group2{iP} = {'NaN'};
            if ~isempty(delInx)
                if delInx(end)~=iP
                    delInx = [delInx, iP];
                end
            else
                delInx = [delInx, iP];
            end
        end
        
        % Chat groups based on path
        if strcmp(PATH(2:end), 'other') 
            chatID(1,1) = getvalue('ChAT+', cell1);
            pchatID(1,1) = getvalue('pChAT+', cell1);
            chatID(1,2) = getvalue('ChAT+', cell2);
            pchatID(1,2) = getvalue('pChAT+', cell2);
            if sum(chatID(1,1))>0 || sum(pchatID(1,1))>0
                if sum(chatID(1,2))>0 || sum(pchatID(1,2))>0
                    chatProp = 'bothChAT/';
                    ChGroup1 = '';
                    ChGroup2 = '';
                else
                    chatProp = 'oneChAT/';
                    ChGroup1 = 'Ch: ';
                    ChGroup2 = '';
                    chatIndex = [chatIndex; iP];
                end
            else
                if sum(chatID(1,2))>0 || sum(pchatID(1,2))>0
                    chatProp = 'oneChAT/';
                    ChGroup1 = '';
                    ChGroup2 = 'Ch: ';
                    chatIndex = [chatIndex; iP];
                else
                    chatProp = 'noChAT/';
                    ChGroup1 = '';
                    ChGroup2 = '';
                end
            end
        elseif strcmp(PATH(2:end), 'ACx')
            chatID1 = getvalue('Area1', cell1);
            chatID2 = getvalue('Area1', cell2);
            if strcmp(chatID1{1}, 'ACx') && strcmp(chatID2{1}, 'ACx')
                chatProp = 'bothACx/';
                ChGroup1 = 'ACx ';
                ChGroup2 = 'ACx ';
            elseif ~strcmp(chatID1{1}, 'ACx') && ~strcmp(chatID2{1}, 'ACx')
                chatProp = 'bothNB/';
                ChGroup1 = 'NB ';
                ChGroup2 = 'NB ';
                chatIndex = [chatIndex; iP];
            elseif strcmp(chatID1{1}, 'ACx') && ~strcmp(chatID2{1}, 'ACx')
                chatProp = 'mixNBACx/';
                ChGroup1 = 'ACx ';
                ChGroup2 = 'NB ';
            elseif ~strcmp(chatID1{1}, 'ACx') && strcmp(chatID2{1}, 'ACx')
                chatProp = 'mixNBACx/';
                ChGroup1 = 'NB ';
                ChGroup2 = 'ACx ';
            end
        else
            chatProp = '';
            ChGroup1 = '';
            ChGroup2 = '';
        end
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.8, [ChGroup1 group1{iP}{1} ' ' ChGroup2 group2{iP}{1}])
%         title([g.burstfilter ' cb: ' cellbase])
        [~, Label] = groupCellPair(group1{iP}, group2{iP});
        axis square;
        setmyplot_tamas;
        x_lim = xlim;
        xticks([x_lim(1) 0 x_lim(2)]);
        xlabel('Lag (ms)');
        title('   ');
        if isempty(delInx)
            if g.issave   % save figure
                ncl1 = regexprep(cell1,'\.','_');
                ncl2 = regexprep(cell2,'\.','_');
                fnm = [chatProp Label '_CCG_' ncl1 '_' ncl2 '_' g.burstfilter '_' g.cells '_new.fig'];
                fnm2 = [chatProp Label '_CCG_' ncl1 '_' ncl2 '_' g.burstfilter '_' g.cells '_new.pdf'];
                fnm3 = [chatProp Label '_CCG_' ncl1 '_' ncl2 '_' g.burstfilter '_' g.cells '_new.jpeg'];
                saveas(H1,fullfile(resdir,fnm))   % save CCG plot
                saveas(H1,fullfile(resdir,fnm2))   % save CCG plot
                saveas(H1,fullfile(resdir,fnm3))   % save CCG plot
                close(H1)
            end
        elseif delInx(end)~=iP
            if g.issave   % save figure
                ncl1 = regexprep(cell1,'\.','_');
                ncl2 = regexprep(cell2,'\.','_');
                fnm = [chatProp Label '_CCG_' ncl1 '_' ncl2 '_' g.burstfilter '_' g.cells '_new.fig'];
                fnm2 = [chatProp Label '_CCG_' ncl1 '_' ncl2 '_' g.burstfilter '_' g.cells '_new.pdf'];
                saveas(H1,fullfile(resdir,fnm))   % save CCG plot
                saveas(H1,fullfile(resdir,fnm2))   % save CCG plot
                saveas(H1,fullfile(resdir,fnm3))   % save CCG plot
                close(H1)
            end
        end
        CCR(iP,:) = ccr;   % cross-correlogram
        LCCR(iP,:) = lwr;  % lower significance limit
        UCCR(iP,:) = upr;  % upper significance limit
        MeanH0(iP) = mean(mean(rccg,2),1);   % surrogate mean
        SDH0(iP) = mean(std(rccg,[],2),1);   % surrogate SD
        disp(['Pair #' num2str(iP) ' / ' num2str(numPairs) ' done......'])
    else
        disp('Not enough spike for CCG');
        delInx = [delInx, iP];
    end
    
    % Save
    if isempty(delInx)
        if g.issave
            if isequal(mod(iP,50),0)   % save after every 50 pairs to prevent data loss
                save(fullfile(resdir,fnm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength', 'group1', 'group2')
                if ~isempty(chatProp)
                    save([resdir chatProp 'chatIndex.mat'],'chatIndex');
                end
                disp('Autosave done.')
            end
        end
    elseif delInx(end)~=iP
        if g.issave
            if isequal(mod(iP,50),0)   % save after every 50 pairs to prevent data loss
                save(fullfile(resdir,fnm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength', 'group1', 'group2')
                if ~isempty(chatProp)
                    save([resdir chatProp 'chatIndex.mat'],'chatIndex');
                end
                disp('Autosave done.')
            end
        end
    end
    waitbar(iP/numPairs)
end
close(wb)   % eliminate progress indicator
if exist(chatProp, 'var')
    save([resdir chatProp 'chatIndex.mat'],'chatIndex');
end

if ~isempty(repInx) % delete duplicated cellpairs
    repInx = unique(repInx);
    CCR(repInx, :) = [];
    group1(repInx) = [];
    group2(repInx) = [];
    LCCR(repInx,:) = [];
    MeanH0(repInx,:) = [];
    PairOfCells(repInx,:) = [];
    SDH0(repInx) = [];
    SegmentLength(repInx) = [];
    UCCR(repInx,:) = [];
end
if ~isempty(delInx) % delete empty datarows
    CCR(delInx, :) = [];
    group1(delInx) = [];
    group2(delInx) = [];
    LCCR(delInx,:) = [];
    MeanH0(delInx,:) = [];
    PairOfCells(delInx,:) = [];
    SDH0(delInx) = [];
    SegmentLength(delInx) = [];
    UCCR(delInx,:) = [];
end

% Save
if g.issave
    save(fullfile(resdir,fnmm),'PairOfCells','CCR','LCCR','UCCR','MeanH0','SDH0','SegmentLength', 'group1', 'group2')
end

