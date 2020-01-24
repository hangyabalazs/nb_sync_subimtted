function [groupID, Label] = groupCellPair(group1, group2)
% GROUPCELLPAIR Defines TP group identity for cellpair
%   [GROUPID, LABEL] = GROUPCELLPAIR(GROUP1, GROUP2) Defines TP group
%   for input cellpair.
%
%   Output arguments:
%       GROUPID - groupID numeric codes
%       LABEL - TP name stored as string
%
%       See also: GROUPCCG, STA_ERS

%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

if iscell(group1)
   group1 = group1{1}; 
end
if iscell(group2)
    group2 = group2{1};
end

if strcmp(group1, '') || strcmp(group2, '')
    groupID = 7;
    Label = 'unknown';
end



if strcmp(group1, group2) % In case of the pairs group ID is the same
    switch group1
        case 'phasicB'
            groupID = 1;
            Label = 'phasicB';
        case 'poissonL'
            groupID = 2;
            Label = 'poissonL';
        case 'tonic'
            groupID = 3;
            Label = 'tonic';
    end
else % Mixed pair
    if strcmp(group1, 'phasicB') & strcmp(group2, 'poissonL') | strcmp(group1, 'poissonL') & strcmp(group2, 'phasicB')
        groupID = 4;
        Label = 'BNB';
    elseif strcmp(group1, 'phasicB') & strcmp(group2, 'tonic') | strcmp(group1, 'tonic') & strcmp(group2, 'phasicB')
        groupID = 5;
        Label = 'BT';
    elseif strcmp(group1, 'tonic') & strcmp(group2, 'poissonL') | strcmp(group1, 'poissonL') & strcmp(group2, 'tonic')
        groupID = 6;
        Label = 'NBT';
    end
end