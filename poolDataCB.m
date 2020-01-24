function poolDataCB

onlyHDB = false;
global RESDIR;

acgNB = load('D:\_MATLAB_DATA\NB\acgNB\ACG_matrices_NB.mat');
acgHDB = load('D:\_MATLAB_DATA\HDB\acgHDB\ACG_matrices_HDB.mat');
acgPannaHDB  = load('D:\_MATLAB_DATA\PannaHDB\acgPannaHDB\ACG_matrices_PannaHDB.mat');
choosecb('PannaHDB');
Pareas = getvalue('Area1', acgPannaHDB.cellids);
allInx = ~strcmp(Pareas, 'HDB');
vpInx = strcmp(Pareas, 'VP');
accInx = strcmp(Pareas, 'AcbC');
allInx(vpInx) = 0;
allInx(accInx) = 0;
pRev = acgPannaHDB.cellids(allInx);

if onlyHDB
    acgNB.cellids = [];
    acgNB.CCR = [];
    acgNB.SCCR = [];
    acgNB.SegmentLength = [];
    acgNB.BurstIndex = [];
    acgNB.Refractory = [];
    acgNB.ThetaIndex = [];
    acgNB.groupID = [];
end


% Concatenate Data
cellids = [acgNB.cellids, acgHDB.cellids, acgPannaHDB.cellids];
CCR = [acgNB.CCR; acgHDB.CCR; acgPannaHDB.CCR(:,501:1500)];
SCCR = [acgNB.SCCR; acgHDB.SCCR; acgPannaHDB.SCCR(:,501:1500)];
SegmentLength = [acgNB.SegmentLength; acgHDB.SegmentLength; acgPannaHDB.SegmentLength];
BurstIndex = [acgNB.BurstIndex; acgHDB.BurstIndex; acgPannaHDB.BurstIndex];
Refractory = [acgNB.Refractory; acgHDB.Refractory; acgPannaHDB.Refractory];
ThetaIndex = [acgNB.ThetaIndex; acgHDB.ThetaIndex; acgPannaHDB.ThetaIndex];
lags = acgHDB.lags;

lengthCB = length(acgPannaHDB.cellids)-length(pRev);
% Remove Unnecessary VP cells
for iR = 1:length(pRev)
    revCell(iR) = find(strcmp(cellids, pRev{iR}));
end
cellids(revCell) = [];
CCR(revCell,:) = [];
SCCR(revCell,:) = [];
SegmentLength(revCell) = [];
BurstIndex(revCell) = [];
Refractory(revCell) = [];
ThetaIndex(revCell) = [];

lastCB = (length(cellids)-(lengthCB-1)):length(cellids);
% TP group indexes
phasicB = Refractory < 40 & BurstIndex > 0.2;   % phasic, bursting
poissonL = Refractory < 40 & BurstIndex <= 0.2;   % poisson-like
tonic = Refractory >= 40;   % refractory above 40 ms
groupID_HP = phasicB(lastCB) + (poissonL(lastCB)*2) + (tonic(lastCB)*3);
groupID = [acgNB.groupID'; acgHDB.groupID'; groupID_HP];

% Define results directory
resdir = [RESDIR 'POOLED\acg\'];
fnmm = 'ACG_matrices_POOLED.mat';
save(fullfile(resdir,fnmm),'cellids','CCR','SCCR','lags',...
    'SegmentLength','BurstIndex','Refractory','ThetaIndex', 'groupID')
