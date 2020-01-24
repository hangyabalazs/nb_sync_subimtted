function groupStaErs(varargin)
%GROUPSTAERS   Sta, spectrum and phase averages by TP groups.
%   GROUPSTAERS Loads all data file from input folder to a struct.
%
%   See also STA_ERS, FIVEBAR
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Input arguments
prs = inputParser;
addParameter(prs,'sessiontype','behavior',@(s)ischar(s)&ismember(s,{'behavior',...
    'burstOn', 'stimHIT', 'stimFA', 'CR', 'Miss', 'CueHIT', 'CueFA', 'CRhit',...
    'Missfa', 'fb_exclude', 'fb_exclude_new'}))   % Plot CCG and corresponding ACG's together
addParameter(prs,'burstfilter','none',@(s)ischar(s)&ismember(s,{'none',...
    'burst1', 'burstall', 'single' 'burst1Single', 'all'}))   % Plot CCG and corresponding ACG's together
parse(prs,varargin{:})
g = prs.Results;

% Initialize
dbstop if error;
global RESDIR;
global cCode;
choosecb('NB')
cb = whichcb;
fs = filesep;
resdir = [RESDIR cb '\sta' cb '\_ALLDATA'];
cd(resdir);   % Set current folder
if iscell(g.sessiontype)
    sessiontype = g.sessiontype;
else
    sessiontype = {g.sessiontype};
end
if strcmp(g.burstfilter, 'all')
    bFilt = {'none', 'burst1', 'single', 'burst1Single'};
else
    bFilt = {g.burstfilter};
end


for sI = 1:length(sessiontype) % Session cycle
    switch sessiontype{sI}
        case 'behavior'
            burstfilter = bFilt;
            groupNames = {'Bursting', 'PoissonL', 'Tonic'};
        case 'burstOn'
            burstfilter = {'none'};
            groupNames = {'All'};
        case 'fb_exclude_new'
            burstfilter = bFilt;
            groupNames = {'Bursting', 'PoissonL', 'Tonic'};
        otherwise
            burstfilter = {'none'};
            groupNames = {'Bursting', 'PoissonL', 'Tonic'};
    end
    for bI = 1:length(burstfilter) % Burstfilter cycle
        if ~strcmp(burstfilter{bI}, 'none')
            if strcmp(burstfilter{bI}, 'burst1Single')
                sync = {'sync', 'nonsync'};
            else
                sync = {'original'};
            end
            groupNames = {'Bursting'};
        else
            sync = {'original'};
        end
        for sY = 1:length(sync) % Sync cycle
            filedir = [RESDIR cb fs 'sta' cb fs burstfilter{bI} ...
                cb fs 'SpectData' fs regexprep(sessiontype{sI},' ','_')...
                fs sync{sY} fs]; % Load corresponding path
            listing = dir(filedir); % List folders data
           
            if length(listing)>2 % if folder is not empty
                %Preallocate
                cellids = [];
                Groups = [];
                Filters = [];
                STA = [];
                SE = [];
                groupDATA = [];
                phase = [];
                spect = [];
                erspLog = [];
                invertList = [];
                
         
                % Load all sta data for every cellbase
                % Find latest sta file
                numVer = '20';
                
                jF=0; % index for file number
                for iF = 1:length(listing) % Read in data from folder
                    if strcmp(listing(iF).name(1), 'n') &&...
                            contains(listing(iF).name, ['_new' numVer]) &&...
                            ~contains(listing(iF).name, 'stAllMatrix')
                        jF = jF+1;
                        fname = [listing(iF).folder '\' listing(iF).name];
                        load(fname); % Load current file
                        fName = listing(iF).name;
                        underS = strfind(fName, '_');
                        dataType = fName(underS(end-1)+2:underS(end)-7); % Define datatype
                        cellids{jF}=cellid;
                        if exist('inverted')
                            invertList(jF) = inverted;
                            clear inverted;
                        end
                        if strcmp(cellid(3), '7')
                            fakeIndex(jF) = 1; % cells without cortical lfp
                        else
                            fakeIndex(jF) = 0;
                        end
                        Groups(jF) = groupID;
                        switch dataType % sta, spect or phase data
                            case 'sta'
                                groupDATA.se{jF} = SE;
                                indexes = find(cellfun(@isempty,groupDATA.se));
                                groupDATA.se(indexes)=[];
                                currData = 'sta';
                            case 'phase'
                                currData = 'ersph';
                            case 'spect'
                                currData = 'I';
                        end
                        % Load current data to struct
                        groupDATA.(dataType){jF} = eval(currData);
                        indexes = find(cellfun(@isempty,groupDATA.(dataType)));
                        groupDATA.(dataType)(indexes)=[]; % Delete empty data
                    end
                end
                
                % Delete redundant information and empty cells
                % (cellids, Groups stored 3 time -because 3 datatype used
                cellids = cellids(1:3:end);
                invertList = invertList(3:3:end);
                Groups = Groups(1:3:end);
                fakeIndex = find(fakeIndex(1:3:end));
                groupDATA.phase(fakeIndex) = [];
                groupDATA.spect(fakeIndex) = [];
                groupDATA.sta(fakeIndex) = [];
                groupDATA.se(fakeIndex) = [];
                cellids(fakeIndex) = [];
                invertList(fakeIndex) = [];
                Groups(fakeIndex) = [];
                STA = cell2mat(groupDATA.sta');
                SE = cell2mat(groupDATA.se');
                
                %Tonic-Phasic groups
                phasicB = Groups == 1;
                poissonL = Groups == 2;
                tonic = Groups == 3;
                BPL = Groups == 4;
                BT = Groups == 5;
                PLT = Groups == 6;
                phasicBPL = Groups == 1 | Groups == 4;
                
                % STA
                wn = size(STA,2)-1; %Plot window size
                time = linspace(-wn/2,wn/2,length(STA));
                if strcmp(sessiontype{sI}, 'burstOn')
                    plotType = 'burstOn';
                    plotGroups = {phasicB, poissonL, tonic};
                    H1 = staPlotter(STA, SE, time, plotGroups, cCode, wn, plotType);
                else
                    plotType = 'groupSTA';
                    plotGroups = {phasicB, poissonL, tonic};
                    H1 = staPlotter(STA, SE, time, plotGroups, cCode, wn, plotType);
                end
                title([sessiontype{sI}])
                saveas(H1, [resdir '\' regexprep(sessiontype{sI},' ','_')...
                    '\' burstfilter{bI} '\' sync{sY} '\' plotType '_new' numVer '.fig']);
                saveas(H1, [resdir '\' regexprep(sessiontype{sI},' ','_')...
                    '\' burstfilter{bI} '\' sync{sY} '\' plotType '_new' numVer '.pdf']);
                                saveas(H1, [resdir '\' regexprep(sessiontype{sI},' ','_')...
                    '\' burstfilter{bI} '\' sync{sY} '\' plotType '_new' numVer '.jpeg']);
                close(H1);
                
                % ERS
                load([RESDIR cb fs 'sta' cb fs 'scaleVector67.mat']);
                scaleVector = f;
                fRange = size(groupDATA.spect{1},1);
                [~,uppLim] = min(abs(f-100));
                f = f(uppLim:fRange); % cuts frequency range at 100Hz
                fValues = round(f, 1); % Frequency scale
                [~,uppGamma] = min(abs(fValues-40));
                [~,uppTheta] = min(abs(fValues-12));
                [~,uppDelta] = min(abs(fValues-4));
                fPos = [1 uppGamma uppTheta uppDelta length(f)];
                fTicks = [fValues(1) fValues(uppGamma) fValues(uppTheta)...
                    fValues(uppDelta) fValues(end)];
                
                % Loop through cells
                deltainx = find(f>0.5&f<4);
                thetainx = find(f>4&f<13);
                gammainx = find(f>13&f<40);
                for aP = 1:length(groupDATA.phase)
                    groupDATA.phase{aP} = groupDATA.phase{aP}(uppLim:end,:);
                    groupDATA.spect{aP} = groupDATA.spect{aP}(uppLim:end,:);
                    currHz = groupDATA.spect{aP}; % spect data
                    centerP = round(size(currHz, 2)/2); % center for plotting
                    deltaS{aP} = mean(currHz(deltainx,centerP:centerP+50)); % delta range <4Hz
                    thetaS{aP} = mean(currHz(thetainx,centerP:centerP+50)); % theta range ~4Hz-12Hz
                    gammaS{aP} = mean(currHz(gammainx,centerP:centerP+50)); % gamma range ~12Hz-41Hz
                end
                erspLog = groupDATA.spect;
                cLim = [-0.4 0.4];
                xrange = size(erspLog{1},2)-1;
                
                % Spectrum
                spectB = erspLog(phasicB);
                spectP = erspLog(poissonL);
                spectT = erspLog(tonic);
                
                % Phase
                phase = groupDATA.phase;
                phaseB = phase(phasicB);
                phaseP = phase(poissonL);
                phaseT = phase(tonic);
                
                % Sync change of phasicB
                if ~strcmp(sync{sY}, 'original')
                    phasicB = Groups == 1 | Groups == 4;
                end
                
                % Data for Bar plots
                deltaAll = cell2mat(deltaS');
                thetaAll = cell2mat(thetaS');
                gammaAll = cell2mat(gammaS');
                deltaB = deltaAll(phasicB,:);
                deltaP = deltaAll(poissonL,:);
                deltaT = deltaAll(tonic,:);
                thetaB = thetaAll(phasicB,:);
                thetaP = thetaAll(poissonL,:);
                thetaT = thetaAll(tonic,:);
                gammaB = gammaAll(phasicB,:);
                gammaP = gammaAll(poissonL,:);
                gammaT = gammaAll(tonic,:);
                meanDB = mean(deltaB,2)';
                meanTB = mean(thetaB,2)';
                meanGB = mean(gammaB,2)';
                meanDP = mean(deltaP,2)';
                meanTP = mean(thetaP,2)';
                meanGP = mean(gammaP,2)';
                meanDT = mean(deltaT,2)';
                meanTT = mean(thetaT,2)';
                meanGT = mean(gammaT,2)';
                stackedB = [meanDB;  meanTB;  meanGB];
                stackedP = [meanDP; meanTP;  meanGP];
                stackedT = [meanDT;  meanTT;  meanGT];
                % Save bar data
                barData.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).Bursting = stackedB;
                barData.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).Poisson = stackedP;
                barData.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).Tonic = stackedT;
            end
            
            % Spect and Phase plots
            for gN = 1:length(groupNames)
                Bursting = mean(cat(3,spectB{:}),3);
                PoissonL = mean(cat(3,spectP{:}),3);
                Tonic = mean(cat(3,spectT{:}),3);
                All = mean(cat(3,erspLog{:}),3);
                currDataS = eval(groupNames{gN});
                if ~isempty(currDataS)
                    
                    % SPECT
                    H2 = figure;
                    imagesc(currDataS);
                    yticks(fPos);
                    b_rescaleaxis('Y',round(fValues))
                    setappdata(gca,'scaley',round(fValues))
                    b_zoomset_for_wavelet
                    xlim([0 xrange]);
                    xticks([0 xrange/2 xrange]);
                    xticklabels({num2str(-xrange/2), '0', num2str(xrange/2)});
                    colorbar;
                    caxis(cLim);
                    colormap jet;
                    axis square;
                    ylabel('Frequency (Hz)');
                    xlabel('Lag (ms)');
                    title([sessiontype{sI},' ',groupNames{gN}]);
                    setmyplot_tamas;
                    saveas(H2, [resdir '\' regexprep(sessiontype{sI},' ','_')...
                        '\' burstfilter{bI} '\' sync{sY} '\_new\' groupNames{gN} 'CellsSpect_new' numVer '.fig']);
                    saveas(H2, [resdir '\' regexprep(sessiontype{sI},' ','_')...
                        '\' burstfilter{bI} '\' sync{sY} '\_new\' groupNames{gN} 'CellsSpect_new' numVer '.pdf']);
                    saveas(H2, [resdir '\' regexprep(sessiontype{sI},' ','_')...
                        '\' burstfilter{bI} '\' sync{sY} '\_new\' groupNames{gN} 'CellsSpect_new' numVer '.jpeg']);
                    close(H2);
                else
                    disp(['No data for the current group: ' groupNames{gN}])
                end
                
            end
            
            % Store data in struct
            master.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).sta = STA;
            master.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).se = SE;
            master.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).phase = phase;
            master.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).spect = erspLog;
            master.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).groups = Groups;
            master.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).cellids = cellids;
            master.(sessiontype{sI}).(burstfilter{bI}).(sync{sY}).invertList = invertList;
        end
    end
end
keyboard;
invertListB = master.(sessiontype{sI}).none.original.invertList;
% Save sta master file and bardata
save([resdir fs 'masterData_new' numVer '.mat'], 'master', '-v7.3');
save([resdir fs 'barData' numVer '.mat'], 'barData', '-v7.3');
save([resdir fs 'invertList' numVer '.mat'], 'cellids', 'invertListB', '-v7.3');


%--------------------------------------------------------------------------

function H = staPlotter(STA, SE, time, plotGroups, cCode, wn, plotType)

H = figure;
hold on;

switch plotType
    case 'groupSTA' % STA plot for TP groups
        phasicB = plotGroups{1};
        poissonL = plotGroups{2};
        tonic = plotGroups{3};
        errorshade(time,mean(STA(phasicB,:), 1), std(STA(phasicB,:), [], 1)/sqrt(sum(phasicB)),'LineColor',cCode(1,:),'ShadeColor',cCode(1,:),'LineWidth',1.5)
        errorshade(time,mean(STA(poissonL,:), 1),std(STA(poissonL,:), [], 1)/sqrt(sum(poissonL)),'LineColor',cCode(2,:),'ShadeColor',cCode(2,:),'LineWidth',1.5)
        errorshade(time,mean(STA(tonic,:), 1),std(STA(tonic,:), [], 1)/sqrt(sum(tonic)),'LineColor',cCode(3,:),'ShadeColor',cCode(3,:),'LineWidth',1.5)
        xlim([-wn/2 wn/2]);
        xticks([-wn/2 0 wn/2]);
        xlabel('Lag (ms)');
        y_lim=ylim;
        y_lim = round(y_lim,1,'significant');
        ylim(y_lim);
        yticks([y_lim(1) 0 y_lim(2)]);
        ylabel('uV');
        title('   ');
        axis square;
        setmyplot_tamas;
    case 'burstOn' % STA plot for allcells
        errorshade(time,mean(STA, 1), std(STA, [], 1)/sqrt(size(STA,1)),'LineColor',[0 0 0],'ShadeColor',[0 0 0],'LineWidth',1.5)
        xlim([-wn/2 wn/2]);
        xticks([-wn/2 0 wn/2]);
        xlabel('Lag (ms)');
        y_lim=ylim;
        y_lim = round(y_lim,1,'significant');
        ylim(y_lim);
        yticks([y_lim(1) 0 y_lim(2)]);
        ylabel('uV');
        title('   ');
        axis square;
        setmyplot_tamas;
end