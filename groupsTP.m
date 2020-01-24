function [phasicB, poissonL, tonic] = groupsTP(cellIDs)
%GROUPSTP   Defines TP group identity for the input cellids.
%   [PHASICB, POISSONL, TONIC] = GROUPSTP(CELLIDS) Defines TP identity
%   from cellids.
%
%   Output arguments:
%       PHASICB - indexes of cells with Ref < 40 & BI > 0.2
%       POISSONL - indexes of cells with Ref < 40 & BI <= 0.2
%       TONIC - indexes of cells with refractory >= 40 ms
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

cb = whichcb;
global RESDIR;
global PATH;
fs = filesep;

if ~iscell(cellIDs)
   cellIDs = {cellIDs}; 
end
% Load BurstIndex and Refractory for the current cells
BI = getvalue('BurstIndex',cellIDs);
Ref = getvalue('Refractory',cellIDs);
if iscell(BI)
    BI = cell2mat(BI);
    Ref = cell2mat(Ref);
end

if isempty(BI) || isempty(Ref)
    % If data has not been loaded to cellbase load from acg matrix
    acgDir = [RESDIR cb fs 'acg' cb PATH fs];
    acgPath = [acgDir 'ACG_matrices_' cb '.mat'];
    load(acgPath);
    for iC = 1:length(cellIDs)
        cellIDS = cellIDs(iC);
        for iCS = 1:length(cellids) % Cellids for the current cellbase
            if strcmp(cellids{iCS}, cellIDS)
                index = iCS;
            end
        end
        BI(iC) = BurstIndex(index);
        Ref(iC) = Refractory(index);
    end
end

% TP group indexes
phasicB = Ref < 40 & BI > 0.2;   % phasic, bursting
poissonL = Ref < 40 & BI <= 0.2;   % poisson-like
tonic = Ref >= 40;   % refractory above 40 ms
