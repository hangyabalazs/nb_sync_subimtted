function pathSelect(cells)
%PATHSELECT   Selects a path for saving for the current cells. More cell
%   types are used in one cellbase so path is a subcategory inside the 
%   current cellbase.
%   RAWTRACEPLOT(CELLS) Creates global variable PATH based on the input
%   cells.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

global PATH
switch cells % options: NB, HDB, ACb, CPu, PV, other
    case 'NB' % nucleus basalis
        PATH = '';
    case 'HDB' % horizontal limb of the diagonal band of Broca
        PATH = '';
    case 'ACb' % Accumbens
        PATH = '\ACb';
    case 'CPu' % Caudate putamen
        PATH = '\CPu';
    case 'PV' % PV cells in the NB
        PATH = '\PV';
    case 'other' % unidentified cells in the NB
        PATH = '\other';
    case 'mixNB' % unidentified and cholinergic cells in the NB
        PATH = '\mixNB';
    case 'mixPV' % unidentified and PV cells
        PATH = '\mixPV';
    case 'otherHDB' % unidentified HDB cells
        PATH = '\otherHDB';
    case 'ACx' % unidentified cells from the auditory cortex
        PATH = '\ACx';
    case 'ChAT'
        PATH = '\ChAT';
    case 'pChAT'
        PATH = '\pChAT';
    case 'PannaHDB'
        PATH = '';
end