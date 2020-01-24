function [ersp pers] = ers(vdisc,pow,sr,wn)
%ERS   Event-related spectrogram.
%   [ERSP ERSPA] = ERS(EVENTS,POW,SR,WN) calculates event-related spectrogram.
%   Input arguments:
%       EVENTS: aligning events (in data points)
%       POW: sectrum
%       SR: sampling rate
%       WN (optional): window size (in data points; default: 5 sec)
%   Output arguments:
%       ERSP: event-related spectrogram
%       ERSPA: all aligned spectra in a 3D array (1st dim corresponds to events)
%
%   See also ASTANORM2

% Input argument check
narginchk(3,4)
if nargin < 4
    wn = 1 * sr;
end

% Calculate ERS
wn2 = round(wn/2);
lenv = length(vdisc);
pers = zeros(lenv,size(pow,1),2*wn2+1);
for t = 1:lenv
    pow2 = pow(:,vdisc(t)-wn2:vdisc(t)+wn2);
    pers(t,:,1:2*wn2+1) = pow2;
end
ersp = squeeze(mean(pers,1));