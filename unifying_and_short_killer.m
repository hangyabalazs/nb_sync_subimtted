function [theta, s1, s2] = unifying_and_short_killer(theta, windowS1, windowS2)
%UNIFYING_AND_SHORT_KILLER() finds and kills too short segments in logical 
%vector 'theta'.
%
%   Author: Barnabas Kocsis
%   Institute of Experimental Medicine, MTA
%   Date: 03/08/2018

if nargin < 3
    windowS1 = 5*nsr; %window size (neglect delta segments shorter than this value)
    windowS2 = windowS1; %window size (neglect theta segments shorter than this value)
end

%unifier (kills short delta segments):
dtheta = diff(theta);
s1 = find(dtheta==1);  % change from non-theta to theta.
s2 = find(dtheta==-1);  % change from theta to non-theta

for it = 1:length(s1)-1
    if s1(it+1)-s2(it) < windowS1
        theta(s2(it):s1(it+1)) = 1;
    end
end

%short_killer (kills short theta segments):
theta = [0 theta 0];
dtheta = diff(theta);
s1 = find(dtheta==1);  % change from non-theta to theta.
s2 = find(dtheta==-1);  % change from theta to non-theta

for it = 1:length(s1)-1
    if s2(it)-s1(it) < windowS2
        theta(s1(it):s2(it)) = 0;
    end
end

theta = theta(3:end-2);
dtheta = diff(theta);
s1 = find(dtheta==1);
s2 = find(dtheta==-1);
end