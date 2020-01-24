function initGlobals(cb)
%INITGLOBALS   Initialize global variables.
%   INITGLOBALS(CELLBASE) Current cellbase used as input argument.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neuroscience
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Directory for saving results
global RESDIR;
RESDIR = 'D:\_MATLAB_DATA\';

% Initialize random seed
rng('shuffle');

% Update current CellBase
pathSelect(cb);

% Color codes for Bursting, Poisson-like, Tonic cells
global cCode;
cCode = [0.6, 0, 0; 0.9, 0.6, 0.1; 0, 0.6, 0];

% Time interval used for plotting
global plotWN;
plotWN = [-250 250];