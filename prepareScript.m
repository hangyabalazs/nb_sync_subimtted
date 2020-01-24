function [all_ChAT, all_pChAT, allChAT] = prepareScript
%PREPARESCRIPT
% [ALL_CHAT, ALL_PCHAT, ALLCHAT] = PREPARESCRIPT
% Initialize necessary variables for master wrapper
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu
%% Prepare
% Global variables
cb = 'NB'; % select cellbase for nucleus basalis
currCB = whichcb;
if ~strcmp(cb, currCB)
    choosecb(cb);
end

% Initialize global variables
initGlobals(cb); 

% Select cholinergic cells for ALL CELLBASES
% NB
[nbChAT, nbpChAT, nballChAT] = selectChAT(cb);

% HDB
choosecb('HDB');
[hdbChAT, hdbpChAT, hdballChAT] = selectChAT('HDB');

% PannaHDB
% Currently HDB, SI, VP, ACB cells are also added
choosecb('PannaHDB');
[PhdbChAT, PhdbpChAT, PhdballChAT] = selectChAT('PannaHDB');

% Pool all the cellids from the cellbases
all_ChAT = [nbChAT'; hdbChAT'; PhdbChAT'];
all_pChAT = [nbpChAT'; hdbpChAT'; PhdbpChAT'];
allChAT = [nballChAT'; hdballChAT'; PhdballChAT'];

