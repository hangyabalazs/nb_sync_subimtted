function cbSwitcher(cellid)
% CBSWITCHER It switches between cellbases based on the input cellid

% Get current cellbase
currCB = whichcb;

% Store HDB cells experiment codes
hdbCells = {'n067'; 'n070'; 'n078'};

% 
% Convert cellid from cell to string
if iscell(cellid)
    cellid = cellid{1};
end

% Try to figure out current cells cellbase and switch to it if necessary
try
    % HDB
    if any(contains(hdbCells, cellid(1:4)))
        cellCB = 'HDB';
    % PannaHDB
    elseif strcmp(cellid(1:3),'HDB')
        cellCB = 'PannaHDB';
    % NB
    else
        if ~strcmp(currCB,'NB')
            choosecb('NB');
        end
        Area = getvalue('Area1', cellid);
        nbAreas = {'GP','GP/SI','SI','IC','RT/IC','EP','EA','EAC'};
        if ismember(Area,nbAreas)
            cellCB = 'NB';
        elseif ismember(Area,{'CPu', 'Cpu'})
            cellCB = 'NB';
        elseif ismember(Area,{'Acb', 'ACb'})
            cellCB = 'HDB';
        else % In case of unknown cellid terminate
            disp('CellID comes from unknown CB')
            keyboard;
            cellCB = currCB;
        end
    end
catch % In case of unknown cellid terminate
    disp('CellID comes from unknown CB')
    keyboard;
end

% Switch cellbase if necessary
if ~strcmp(currCB, cellCB)
    choosecb(cellCB);
end



% NBcells =  {'n028_120211a_8.1';'n028_120211a_8.2';'n029_120215a_2.2';...
%     'n029_120215a_3.4';'n046_130102b_4.1';'n046_130102b_6.1';...
%     'n046_130102x_4.1';'n046_130102x_4.3';'n046_130103a_6.1';...
%     'n046_130103a_6.2';'n046_130104a_4.1';'n046_130104a_6.1';...
%     'n046_130104a_6.2';'n046_130107a_4.2';'n046_130107a_8.1';...
%     'n046_130108a_4.1';'n046_130108a_4.2';'n046_130108a_8.1';...
%     'n046_130108d_4.1';'n046_130108d_8.2';'n046_130108e_4.1';...
%     'n046_130108e_8.2';'n046_130108f_4.1';'n046_130108f_4.2';...
%     'n071_141201a_2.3';'n071_141201a_7.2'};
%
% HDBcells = {'n067_141017a_1.3';'n067_141017a_4.1'};
%
% PannaHDBcells = {'HDB25_180526a_3.2';'HDB25_180526a_5.1';...
%     'HDB36_190426a_2.1';'HDB36_190426a_3.1';'HDB36_190426a_4.2';...
%     'HDB36_190504a_3.1';'HDB36_190504a_5.1';'HDB36_190511a_3.3';...
%     'HDB36_190511a_5.1';'HDB36_190515a_3.1';'HDB36_190515a_4.1'};

% This list is good for ACG
% allCells = {'n023_111218a_1.2';'n029_120202a_3.5';'n029_120203a_3.1';...
%     'n029_120214b_2.2';'n029_120301a_3.3';'n046_121229a_1.1';...
%     'n046_121230a_1.2';'n046_121231a_1.2';'n046_130101a_6.1';...
%     'n046_130102a_6.1';'n046_130102x_4.1';'n046_130103a_6.2';...
%     'n046_130104a_6.2';'n046_130108a_4.1';'n046_130108a_4.2';...
%     'n071_141129a_2.3';'n071_141201a_2.3';'n071_141201a_7.2';...
%     'n071_141229a_3.2';'n071_141218a_6.1';'n072_141223a_4.2';...
%     'n072_141222a_4.3';'n045_121217x_4.6';'n023_111220a_1.2';...
%     'n029_120207b_1.1';'n029_120210a_3.3';'n028_120211a_8.1';...
%     'n028_120211a_8.2';'n029_120215a_3.4';'n029_120220a_3.1';...
%     'n029_120221b_6.1';'n029_120222b_4.1';'n029_120302a_3.3';...
%     'n029_120313a_1.1';'n029_120314a_3.1';'n037_121006a_4.1';...
%     'n046_121213a_3.1';'n046_121218a_2.2';'n046_121219a_8.1';...
%     'n045_121231a_8.1';'n046_130102x_4.3';'n046_130104a_4.1';...
%     'n046_130104a_6.1';'n046_130108a_8.1';'n077_141218a_5.1';...
%     'n067_141009b_2.2';'n067_141012a_3.2';'n067_141017a_4.1';...
%     'n067_141018a_5.1';'n067_141011x_1.1';'n070_141112a_4.1';...
%     'n070_150104a_5.1';'n078_150104a_1.1';'n078_150110a_3.1';...
%     'n070_150108x_6.1';'n067_141017a_1.3';'n067_141019a_5.2';...
%     'HDB18_170813a_4.2';'HDB25_180420a_2.2';'HDB25_180421a_2.1';...
%     'HDB25_180526a_3.2';'HDB25_180526a_5.1';'HDB25_180611a_4.1';...
%     'HDB36_190426a_2.1';'HDB36_190426a_3.1';'HDB36_190426a_4.2';...
%     'HDB36_190504a_3.1';'HDB36_190504a_5.1';'HDB36_190505a_3.1';...
%     'HDB36_190507a_3.3';'HDB36_190508a_3.2';'HDB36_190510a_3.2';...
%     'HDB36_190511a_3.3';'HDB36_190511a_5.1';'HDB36_190512a_3.2';...
%     'HDB36_190513a_3.3';'HDB36_190515a_3.1';'HDB36_190515a_4.1';...
%     'HDB36_190517a_3.3';'HDB36_190518a_3.3';'HDB38_190522a_4.1'};


% cellPos = find(contains(allCells, cellid));
%
% if cellPos < 46
%     cellCB = 'NB';
% elseif cellPos > 45 && cellPos < 58
%     cellCB = 'HDB';
% elseif cellPos > 57
%     cellCB = 'PannaHDB';
% else
%     disp('CellID comes from unknown CB')
%     keyboard;
% end

% nbCont = find(contains(NBcells, cellid));
% hdbCont = find(contains(HDBcells, cellid));
% PannahdbCont = find(contains(PannaHDBcells, cellid));
%
% if ~isempty(nbCont)
%     cellCB = 'NB';
% elseif ~isempty(hdbCont)
%     cellCB = 'HDB';
% elseif ~isempty(PannahdbCont)
%     cellCB = 'PannaHDB';
% else
%     disp('CellID comes from unknown CB')
%     keyboard;
% end
%




