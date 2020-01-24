function [ChAT, pChAT, allChAT] = selectChAT(area)
%SELECTCHAT   Select cholinergic neurons.
%   [CHAT, PCHAT, ALLCHAT] = SELECTCHAT(AREA) selects cholinergic neurons
%   from NB (AREA = 'NB', default) of HDB (AREA = 'HDB').
%
%   Output arguments:
%       CHAT - cellIDs of tagged cholinergic neurons
%       PCHAT - cellIDs of putative cholinergic neurons
%       ALLCHAT - cellIDs of both putative and tagged cholinergic neurons
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Cholinergic cells
switch area
    case 'NB'
        ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % identified
        ChAT = [ChAT 'n045_121217x_4.6'];   % clustered based on light-evoked spikes
        pChAT = selectcell(['"pChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})']);  % putative
        
    case 'HDB'
        ChAT = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified ChAT+ cells
        ChAT = [ChAT 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
        ChAT = [ChAT 'n067_141019a_5.2'];   % light spike assisted clustering
        pChAT = [];
        
    case 'other'
        % NB cellbase
        
        % Non-Cholinergic neurons from gonogo, and feedbackdelay
        selstr = ['"pChAT+"==0&"ChAT+"==0&"PV+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''GP'',''GP/SI'',''SI'',''IC'',''RT/IC'',''EP'',''EA'',''EAC''})'];
        ChAT = selectcell(selstr);  % identified
        pChAT = [];
        
    case 'otherHDB'
        
        % HDB Non-Cholinergic neurons
        selstr = ['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''HDB'',''SI'',''VP''})'];  % identified (n = 12)
        ChAT = selectcell(selstr);
        pChAT = [];
    case 'ACx'
        
        % A1 Non-Cholinergic neurons
        selstr = ['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''ACx''})'];  % identified (n = 12)
        ChAT = selectcell(selstr);
        pChAT = [];
    case 'CPu'
        
        % A1 Non-Cholinergic neurons
        selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''CPu'',''Cpu''})'];  % identified (n = 12)
        ChAT = selectcell(selstr);
        pChAT = [];
        
    case 'Acb'
        
        % A1 Non-Cholinergic neurons
        selstr = ['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''Acb'',''ACb''})'];  % identified (n = 12)
        ChAT = selectcell(selstr);
        pChAT = [];
        
    case 'PannaHDB'
        ChAT = selectcell(['"ChAT+"==1&' ...
            'ismember("Area1",{''HDB'',''SI''})']);  % identified ChAT+ cells

        pChAT = [];
        
        
    case 'bothHDB'
        % Select HDB cells from Cell paper
        choosecb('HDB')
        ChAT_HDB = selectcell(['"ChAT+"==1&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified ChAT+ cells
        ChAT_HDB = [ChAT_HDB 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
        ChAT_HDB = [ChAT_HDB 'n067_141019a_5.2'];   % light spike assisted clustering
        
        choosecb('PannaHDB')
        ChAT_PHDB = selectcell(['"ChAT+"==1&' ...
            'ismember("Area1",{''HDB'',''SI''})']);  % identified ChAT+ cells
        
        ChAT = [ChAT_HDB ChAT_PHDB];
        pChAT = [];
    case 'otherBothHDB'
                % Select HDB cells from Cell paper
        choosecb('HDB')
        ChAT_HDB = selectcell(['"ChAT+"==0&"ID_PC">20&"Lr_PC"<0.15&"validity"==1&' ...
            'ismember("session_type",{''gonogo'',''feedbackdelay''})&' ...
            'ismember("Area1",{''HDB'',''SI'',''VP''})']);  % identified ChAT+ cells
        ChAT_HDB = [ChAT_HDB 'n067_141017a_1.3'];   % without miss-triggering it passes the cluster criteria
        ChAT_HDB = [ChAT_HDB 'n067_141019a_5.2'];   % light spike assisted clustering
        
        choosecb('PannaHDB')
        ChAT_PHDB = selectcell(['"ChAT+"==0&' ...
            'ismember("Area1",{''HDB'',''SI''})']);  % identified ChAT+ cells
        
        ChAT = [ChAT_HDB ChAT_PHDB];
        pChAT = [];

end

allChAT = [ChAT pChAT];   % identified and putative