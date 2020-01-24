function acg_NBsync(cellids,varargin)
%ACG   Auto-correlation.
%   ACG calculates auto-correlations. Window size is set to +-500 ms with a
%   0.5 ms resolution. Maximum 50000 spikes are included to avoid memory
%   problems. Segments are filtered with FINDSEGS3. For details on the
%   algorithm, see XCORR.
%
%   Three indeces are derived from the auto-correlogram (ACG). Burst index
%   is calculated as the normalized difference between maximum ACG for lags
%   0-10 ms and mean ACG for lags 180-200 ms. The normalizing factor is the
%   greater of the two numbers, yielding and index between -1 and 1 (see
%   Royer et al., 2012). Theta index is calculated as the normalized
%   difference between mean ACG for a +-25 ms window around the peak within
%   lags 100 and 200 ms and the mean ACG for lags 180-200 and 65-85 ms.
%   Normalization is performed similarly as for the burst index. Refractory
%   is calculated as full width half hight of the central gap in the
%   smoothed ACG. Smoothing is performed by a 10 ms moving average.
%
%   ACG(I) calls ACG for cells defined by the cell ID list (or index set
%   to CELLIDLIST, see CellBase documentation) I.
%   Optional input arguments (parameter-value pairs with default values):
%       'issave', false - controls saving behavior; plots and
%           auto-correlation matrices are saved only if 'issave' is set to
%           true
%
%   Reference:
%   Royer S, Zemelmen BV, Losonczy A, Kim J, Chance F, Magee JC, Buzsaki G
%   (2012) Control of timing, rate and bursts of hippocampal place cells by
%   dendritic and somatic inhibition. Nat Neurosci 15:769-775
%
%   See also CCG and XCORR.

%   Balazs Hangya, Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   hangya.balazs@koki.mta.hu

% Input arguments
prs = inputParser;
addRequired(prs,'cellids',@(s)iscell(s)|iscellstr(s)|ischar(s))
addParameter(prs,'issave',false,@islogical)   % control saving behavior
addParameter(prs,'altSess',false,@islogical)   % specific parameter for alternative sessions
addParameter(prs,'original','',@(s)iscell(s)|iscellstr(s)|ischar(s)) % original cellid in case of altSession
addParameter(prs, 'inVitroCut', false, @(s)islogical(s)|ismember(s,[0 1])) % downsample based on in-vitro lengths
addParameter(prs, 'pooled', false, @(s)islogical(s)|ismember(s,[0 1])) % in case of pooled cellbases
addParameter(prs,'cellbase','',@(s)iscell(s)|iscellstr(s)|ischar(s)) % original cellid in case of altSession
addParameter(prs, 'longSeg', false, @(s)islogical(s)|ismember(s,[0 1])) % downsample based on in-vitro lengths

parse(prs,cellids,varargin{:})
g = prs.Results;
if ischar(cellids)
    cellids = {cellids};  % one cell ID
end

% Pass the control to the user in case of error
dbstop if error

% Choose cb
if g.pooled
    cellbase = 'POOLED';
else
    if ~isempty(g.cellbase)
        cellbase = g.cellbase;
    else
        cellbase = whichcb;   % Check the current cellbase
    end
end

% Select result directory
global RESDIR;
global PATH;
fs = filesep;

% Set inVitro directory
if g.inVitroCut
    vitroDir = ['inVitro' fs];
else
    vitroDir = [];
end

% Alternative sessions
if g.altSess
    resdir = [RESDIR cellbase fs 'altSession' cellbase fs vitroDir];
else
    resdir = [RESDIR cellbase fs 'acg' cellbase PATH fs vitroDir];
end

if isempty(resdir)
    error('Did not select any folder. Program aborts')
else
    cd(resdir);   % Set current folder
end


% Directories
if g.inVitroCut
    fnmm = ['ACG_matrices_' cellbase '_inVitro.mat'];   % filename for saving the result matrices
else
    fnmm = ['ACG_matrices_' cellbase '.mat'];   % filename for saving the result matrices
end

% Include only long enough segments
longsegments = false;  % control whether to use this option
seglim = 300;
nameSave = []; % Creating variable for storing cellids
delInx = [];
% Determine time window
sr = 1000;      % sampling rate
if g.longSeg
    wn = 500 * sr / 1000;    % 2 * 500 ms window
        resdir = [resdir 'longSeg' fs];
else
    wn = 250 * sr / 1000;    % 2 * 500 ms window
end
res = 0.5;   % resolution for ACG in ms

% Cell loop for ACG
wb = waitbar(0,'Please wait...','Name','Running ACG...');  % progress indicator
global WB
WB(end+1) = wb;
limit_spikes = [1 50000];   % include max 50000 spikes; no lower limit
numCells = length(cellids);
[CCR, SCCR] = deal(zeros(numCells,2*wn/res));
[SegmentLength, BurstIndex, Refractory, ThetaIndex] = deal(nan(numCells,1));
for iC = 1:numCells   % loop through the cells
    currCell = cellids{iC};
    cbSwitcher(currCell);
    try
        tseg = findSegs3(currCell,'segfilter','stim_excl_nb',...
            'light_activation_duration',[-5 5],'margins',[0 0]);  % find time segments
        ltseg = tseg(2,:) - tseg(1,:);  % length of the segments
        if longsegments   % later to be implemented as an input option
            seginx = find(ltseg==max(ltseg));
            tseg = tseg(:,seginx(1));   % find the longest segment
            ltseg = ltseg(seginx(1));
            if tseg < seglim   % use the longest segment if it's longer than the threshold
                continue
            end
        end
        SegmentLength(iC) = sum(ltseg);  % cumulative length of the segments
        [ncc,seltsind,selisi] = extractSegSpikes_NBsync(currCell,tseg,'inVitroCut',g.inVitroCut);   % find spikes in the time segments
    catch ME
        disp(ME.message)
        disp('Could not extract the segment.');
    end
    
    %     ncc = loadcb(cell,'SPIKES');   % use all spikes
    
    if length(ncc) > limit_spikes(2)      % crop if too long to avoid out of memory
        ncc = ncc(1:limit_spikes(2));
    end
    
    if length(ncc) > limit_spikes(1)     % minimum criterion
        [H1, ccr, lags] = acorr(ncc,wn,res);
        sccr = smooth(ccr,'linear',21);    % smoothed ACR
        %         nqf = 1 / res * 1000 / 2;   % Nyquist freq.
        %         flt = fir1(32,[4 10]/nqf,'bandpass');
        %         sccr2 = filtfilt(flt,1,sccr);   % high-pass filter > 1 Hz
        hold on
        %         plot(lags,sccr,'Color',[0.7 0.7 0.7])
        %         plot(lags,sccr2,'c')
        %         bar(lags(lags>-10&lags<10),ccr(lags>-10&lags<10), 1, 'FaceColor',[0 0.8 0],'EdgeColor',[0 0.8 0])
        %         bar(lags(lags>-200&lags<-180),ccr(lags>-200&lags<-180), 1, 'FaceColor',[0.8 0 0],'EdgeColor',[0.8 0 0])
        %         bar(lags(lags>180&lags<200),ccr(lags>180&lags<200), 1, 'FaceColor',[0.8 0 0],'EdgeColor',[0.8 0 0])
        %         bar(lags((lags>=180&lags<=200)|(lags>=65&lags<=85)),...
        %             ccr((lags>=180&lags<=200)|(lags>=65&lags<=85)),...
        %             1, 'FaceColor','c','EdgeColor','c');
        BurstIndex(iC) = burstinx(ccr,lags);   % burst index
        Refractory(iC) = refract(ccr,sccr,lags);   % refractory
        ThetaIndex(iC) = thetainx4(sccr,lags);   % theta index
        ThetaIndex_IN(iC) = thetainx4_2(sccr,lags);   % theta index
             
        xl = xlim;   % put the cell IDs in the plot
        yl = ylim;
        
        % Defines Phasic-Tonic group identity
        if Refractory(iC) < 40 && BurstIndex(iC) > 0.2
            group = 'Bursting';
            groupID(iC) = 1;
        elseif Refractory(iC) < 40 && BurstIndex(iC) <= 0.2
            group = 'Poisson-like';
            groupID(iC) = 2;
        else
            group = 'Tonic';
            groupID(iC) = 3;
        end
        
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1)),regexprep(currCell,'_',' '))
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.95,['Burst index: ' num2str(BurstIndex(iC))])
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.90,['Refractory: ' num2str(Refractory(iC))])
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.85,['Theta index: ' num2str(ThetaIndex(iC))])
        text(xl(1)+(xl(2)-xl(1))*0.7,yl(1)+(yl(2)-yl(1))*0.80,group)
        xticks([xl(1) 0 xl(2)]);
        
        if g.altSess % Makes title for defining session type
            if ~isempty(g.original)
                if strcmp(g.original, currCell)
                    title('Original Session')
                else
                    title('Alternative Session')
                end
            end
        end
        setmyplot_tamas;
        %         if max(ccr) >2
        if g.issave   % save figure
            ncl = regexprep(currCell,'\.','_');
            nameSave = [nameSave '_' ncl];
            fnm = [regexprep(group,'-','_') '_ACG_' ncl '.fig'];
            saveas(H1,fullfile(resdir,fnm))   % save ACG plot
            close(H1)
        end
        CCR(iC,:) = ccr;   % cross-correlogram
        SCCR(iC,:) = sccr;   % cross-correlogram
        disp(['Cell #' num2str(iC) ' / ' num2str(numCells) ' done......'])
    end
    
    % Save
    if g.issave
        if isequal(mod(iC,50),0)   % save after every 50 pairs to prevent data loss
            save(fullfile(resdir,fnm),'cellids','CCR','SCCR','lags',...
                'SegmentLength','BurstIndex','Refractory','ThetaIndex',...
                'ThetaIndex_IN', 'groupID')
            disp('Autosave done.')
        end
    end
    
    waitbar(iC/numCells)
end
close(wb)   % eliminate progress indicator

% Delete empty data
if ~isempty(delInx)
    cellids(delInx) = [];
    CCR(delInx,:) = [];
    SCCR(delInx,:) = [];
    SegmentLength(delInx) = [];
    BurstIndex(delInx) = [];
    Refractory(delInx) = [];
    ThetaIndex(delInx) = [];
    ThetaIndexN(delInx) = [];
    ThetaIndexN2(delInx) = [];
    groupID(delInx) = [];
end

% Save
if g.issave
    if g.altSess % Saves cellids to filename in case of altSession
        fnmm = ['ACG_matrices_' cellbase '_' nameSave '.mat'];
    end
    save(fullfile(resdir,fnmm),'cellids','CCR','SCCR','lags',...
        'SegmentLength','BurstIndex','Refractory','ThetaIndex',...
        'ThetaIndex_IN', 'groupID')
end

% -------------------------------------------------------------------------
function [H1 ccr lags] = acorr(ncc,wn,res)

% Calculate spike times in milliseconds
sr = 1000;
nc = ncc * sr;
mn = nc(1);  % only relative spike times count; avoid out of memory
nc = nc - mn;
nc(nc<0.5) = [];  % drop spikes before 0.5 ms(!) to avoid error in line 39
wn2 = wn / 1000;    % window size in seconds

% Auto-correlogram
zunit = zeros(1,round(nc(end)/res)+5);
zunit(round(nc/res)) = 1;
[ccr lags] = xcorr(zunit,zunit,wn2*sr/res);     % 1->2; window: -wn ms - wn ms
ccr(length(ccr)/2+0.5) = [];    % auto-correlation: drop middle bin
lags(length(lags)/2+0.5) = [];
lags = lags * res;   % in ms

% Plot
H1 = figure;
bar(lags,ccr, 2, 'FaceColor','black')
set(gca,'XLim',[-wn wn])

% -------------------------------------------------------------------------
function bi = burstinx(ccr,lags)

a = max(ccr(lags>0&lags<10));
b = mean(ccr(lags>180&lags<200));
bi = (a - b) / max(a,b);   % burst index

% -------------------------------------------------------------------------
function tau = refract(ccr,sccr,lags)

mx = max(sccr);   % peak: from smoothed ACG
halfmx = mx / 2;   % half-max
ccr2 = ccr(lags>0);   % half ACG
lags2 = lags(lags>0);   % corresponding lag values
tau = lags2(find(ccr2>halfmx,1,'first'));   % refractory

% -------------------------------------------------------------------------
function thinx = thetainx4(ccr,lags)

% Theta index
thpeak = ccr(lags>=100&lags<=200);   % ACG first theta peak
[jnk pl] = max(thpeak);
ploc = lags(find(lags>=100,1,'first')+pl-1);   % peak location
mthp = mean(ccr(lags>=ploc-25&lags<ploc+25));

thtrough = ccr(lags>=225&lags<=275);   % ACG first theta trough
mtht = mean(thtrough);
thinx = (mthp - mtht) / max(mthp,mtht);   % theta index


% -------------------------------------------------------------------------
function thinx = thetainx4_2(ccr,lags)

% Theta index
thpeak = ccr(lags>=100&lags<=200);   % ACG first theta peak
[jnk pl] = max(thpeak);
ploc = lags(find(lags>=100,1,'first')+pl-1);   % peak location
mthp = mean(ccr(lags>=ploc-25&lags<ploc+25));

thtrough = ccr((lags>=225&lags<=275)|(lags>=65&lags<=85));   % ACG first theta trough
mtht = mean(thtrough);
thinx = (mthp - mtht) / max(mthp,mtht);   % theta index