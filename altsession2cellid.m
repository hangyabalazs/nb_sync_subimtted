function [alterID, indexOut] = altsession2cellid(currID, cellids)
% ALTSESSION2CELLID Returns the original cellid for the input alternative
%   sessions cellid.
%   [alterID, index] = ALTSESSION2CELLID(CURRID, CELLIDS) Defines the
%   original cellid of currID. 
%
%   Output arguments:
%       ALTERID - Original cellid of the alternative session
%       INDEX - Index of the original cellid in the current cell of IDs
%       See also: CELLID2TPGROUP, CCG


% Defining variables
global PATH;
if isempty(cellids) % If no input list given load cells for current CB
    if isempty(PATH)
        cells = 'NB'; 
    else
        cells = PATH(2:end);
    end
    [~, ~, allChAT] = selectChAT(cells);   % Select cholinergic cells
    NumChAT = length(allChAT);
else
    NumChAT = length(cellids);
    allChAT = cellids;
end
% Load alternative sessions of the current cellid
alternates = getvalue('alternative_sessions', currID);
index = [];
% Loop through cellids to find the alternative session
if iscell(alternates)
    for iA = 1:length(alternates{1})
        if isempty(index)
            for allcells = 1:NumChAT
                if strcmp(alternates{1}{iA}, allChAT{allcells})
                    index = iA; % Store original cellids index
                    indexOut = allcells;
                    break;
                end
            end
        else
            break
        end
    end
end

if isempty(index)
    disp('No cellid found for the given alternative session')
    alterID = {};
    alterID{1} = [];
else
    alterID = alternates{1}(index);
end


