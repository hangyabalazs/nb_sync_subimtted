function [group, groupIDS] = cellid2TPgroup(currID, varargin)
% CELLID2TPGROUP Returns the Tonic/Phasic group identifier for a given
%   cellid.
%   [GROUP, GROUPIDS] = CELLID2TPGROUP(CELLID, VARARGIN) Defines TP group
%   for input cellids. ACGpath can be defined for current cells.
%
%   Output arguments:
%       GROUP - TP name stored as string
%       GROUPIDS - groupID numeric codes
%
%       See also: ALTSESSION2CELLID, CCG

%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

prs = inputParser;
addRequired(prs,'cellid',@(s)iscell(s)|iscellstr(s)|ischar(s))
addParameter(prs,'acgPath','',@(s)isempty(s)|ischar(s))   % ACG matrix path
parse(prs,currID,varargin{:})
g = prs.Results;
dbstop if error;

% Input arguments
if ischar(currID) % convert to cell
    currID = {currID};
end
cellbase = whichcb;
global PATH;
global DATAPATH;
acgPath = g.acgPath;

if ~isempty(acgPath)
    load(acgPath);   % Load ACG matrix
else % Try to load ACG matrix for current cellbase
    try
        fnm_acg = [DATAPATH cellbase '\acg' cellbase PATH '\ACG_matrices_'...
            cellbase '.mat'];
        load(fnm_acg);
    catch
        [acg_file, acg_filedir] = uigetfile([DATAPATH cellbase '\acg'...
            cellbase PATH '\ACG_matrices_' cellbase '.mat'],...
            'Choose ACG matrix!');
        if acg_file < 1
            error('No such file. The program aborted.')
        end
        fnm_acg = fullfile(acg_filedir, acg_file);
        load(fnm_acg);   % Load ACG matrix
    end
end

for iC = 1:length(currID) % Loop through input cells
    for iCS = 1:length(cellids) % Cellids for the current cellbase
        if strcmp(cellids{iCS}, currID{iC})
            index = iCS;
        end
    end
    % If cellid not found check alternative sessions
    if ~exist('index', 'var')
        % Find original cellid for the alternative session
        [~, index] = altsession2cellid(currID(iC), cellids);
    end
    % Define TP group
    if exist('index', 'var')
        [phasicB, poissonL, tonic] = groupsTP(cellids{index});
        % Define tonic/phasic groups
        if phasicB
            group{iC} = 'phasicB';
            groupIDS(iC) = 1;
        elseif poissonL
            group{iC} = 'poissonL';
            groupIDS(iC) = 2;
        elseif tonic
            group{iC} = 'tonic';
            groupIDS(iC) = 3;
        end
    else
        group{iC} = '';
        groupIDS(iC) = NaN;
    end
    clearvars index;
end
