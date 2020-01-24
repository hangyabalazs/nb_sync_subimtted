function [allChAT, NumChAT] = selectPCells(cells, allChAT)
%SELECTPCELLS   Returns the specified original and putative cellids from 
%   the selected cellbase.
%   [ALLCHAT, NUMCHAT] = SELECTPCELLS(CELLS, ALLCHAT) returns the original 
%   and the putative cellids and their number from the selected cellbase 
%   based on the properties defined below.
%
%   See also SELECTCELLS
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

%   Edit log: LT 9/28/2017

switch cells
    case 'NB'
        ChAT2 = selectcell(['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
        pChAT2 = selectcell(['"pChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
        allChAT = [allChAT ChAT2 pChAT2];   % identified and putative
        NumChAT = length(allChAT);   % number of cholinergic cells
    case 'HDB'
        ChAT2 = selectcell(['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified
        allChAT = [allChAT ChAT2];   % identified
        NumChAT = length(allChAT);   % number of cholinergic cells
    case 'CPu'
        ChAT2 = selectcell(['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''Cpu'',''CPu''})']);  % identified
        allChAT = [allChAT ChAT2];   % identified
        NumChAT = length(allChAT);   % number of cholinergic cells
    case 'ACb'
        ChAT2 = selectcell(['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''Acb''})']);  % identified
        allChAT = [allChAT ChAT2];   % identified
        NumChAT = length(allChAT);   % number of cholinergic cells
    case 'PV'
        selstr = ['"ChAT+"==0&"pChAT+"==0&"PV+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
        pChAT = selectcell(selstr);  % identified
        allChAT = [allChAT pChAT];   % identified and putative
        NumChAT = length(allChAT);   % number of PV cells
    case 'other'
        selstr = ['"ChAT+"==0&"pChAT+"==0&"PV+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
        ChAT = selectcell(selstr);  % unidentified
        allChAT = ChAT;
        NumChAT = length(allChAT);   % number of unidentified cells
    case 'NB-ACx' % NB cholinergic cells with ACx cells
        ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
        ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
        pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
        % Alternative sessions
        ChAT2 = selectcell(['"ChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
        pChAT2 = selectcell(['"pChAT+"==2&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
        % A1 Non-Cholinergic neurons from A1
        selstr = ['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''ACx''})'];
        ACx = selectcell(selstr);
        allChAT = [ChAT pChAT ChAT2 pChAT2 ACx];   % NB and A1
        NumChAT = length(allChAT);   % number of cholinergic cells
end