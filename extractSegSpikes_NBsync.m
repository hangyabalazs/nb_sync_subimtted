function [selts, seltsind, selisi] = extractSegSpikes_NBsync(cellid,segs, varargin)
%EXTRACTSEGSPIKES   Extract spikes between specific timestamps.
%   [SELTS, SELTSIND, SELISI] = EXTRACTSEGSPIKES(CELLID,SEGS) extracts
%   timestamps for spikes of a cell (CELLID) that lie between start (first
%   row of SEGS) and end (second row of SEGS) points of given segments.
%   Output arguments:
%       SELTS: timestamps of the extracted spikes
%       SELTSIND: indices of the timestamps with respect to all spikes in
%           the cluster
%       SELISI: interspike intervals where more than one spike was found 
%           within a time segment; NaN is returned when there is only one
%           spike
%
%   EXTRACTSEGSPIKES(SPT,SEGS) overloads the function with numeric input in
%   the first argument, implemented as spike times.
%
%   See also FINDSEGS3 and ABS2RELTIMES.

%   Edit log: SPR 5/10, BH 4/25/12, 5/8/12

% Input argument check
prs = inputParser;
addRequired(prs,'cellid',@(s)iscellid(s)|isnumeric(s))
addRequired(prs,'segs',@isnumeric)
addParameter(prs,'burstfilter','none',@(s)ischar(s)&ismember(s,{'none',...
    'burst1', 'burstall', 'single', 'burst1Single'}))   % filter on burst first spikes or single spikes
addParameter(prs, 'inVitroCut', false, @(s)islogical(s)|ismember(s,[0 1])) 
parse(prs,cellid,segs,varargin{:})
g = prs.Results;


% burstfilter = 'single';   % POSSIBLE VALUES: none, burst1, burstall, single
% Load spikes
if iscellid(cellid)
    switch g.burstfilter
        case 'none'
            spk = loadcb(cellid);   % RAW DATA
        case 'burst1'
            spk = loadcb(cellid);
            [spk, ~, ~] = burstDetect(spk);
            spk = cell2mat(spk);
%             spk = burstStimes_CCG(spk, g.burstfilter);   % BURST'S FIRST SPIKES
        case 'burstall'
            spk = loadcb(cellid);
            [~, ~, spk] = burstDetect(spk);
            spk = cell2mat(spk);
%             spk = burstStimes_CCG(spk, g.burstfilter);   % BURST'S ALL SPIKES
        case 'single'
            spk = loadcb(cellid);
            [~, spk, ~] = burstDetect(spk);
            spk = cell2mat(spk);
%             spk = burstStimes_CCG(spk, g.burstfilter);   % SINGLE SPIKES
        case 'burst1Single'
            spk = loadcb(cellid);
            [burst1, single, ~] = burstDetect(spk);
            burst1 = cell2mat(burst1);
            single = cell2mat(single);
            spk = sort([burst1; single]);
%             spk = burstStimes_CCG(spk, g.burstfilter);   % SINGLE SPIKES AND BURST1 SPIKE
    end
else
    spk = cellid;   % overload extractSegSpikes: implement first argument as spike times
end
if iscell(spk)
    spk = spk{:};
else
    spk = spk(:);
end

% Extract spikes for each segment
seltsind = [];
selisi = [];
for is = 1:size(segs,2)   % segment loop
    ladle = find(spk>segs(1,is)&spk<segs(2,is));
    seltsind = [seltsind; ladle]; %#ok<AGROW>
    if size(ladle,1)>1
        selisi = [selisi; diff(spk(ladle))]; %#ok<AGROW>
    else
        selisi = [selisi; NaN]; %#ok<AGROW>
    end
end
selts = spk(seltsind);

% Downsample spikes according to inVitro numbers
if g.inVitroCut
    selts = inVitro_BI_check(cellid, selts);
end
