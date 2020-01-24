function [burst1Spikes, singleSpikes, burstAllSpikes] = burstDetect(spikeTimes)
%BURSTDETECT   Detects burst first spikes, single spikes, burstall spikes
%   [BURST1, SINGLE, BURSTALL] = BURSTDETECT(ALLSPIKES) Detects ISI-s
%   shorter than 10ms, and the consecutive ones below 15ms.
%
%   See also BURSTRATIO.
%
%   Tamas Laszlovszky
%   Laboratory of Systems Neurosciecnce
%   Hungarian Academy of Sciences
%   laszlovszky.tamas@koki.mta.hu

% Input argument check
if iscell(spikeTimes)
    numCell = length(spikeTimes);
else
    spikeTimes = {spikeTimes};
    numCell = 1;
end

%Preallocate
[burst1Spikes, singleSpikes, burstAllSpikes] = deal(cell(1,numCell));

for iC = 1:numCell   % loop through cells
    spike_times = spikeTimes{iC};
    
    % Inter-spike intervals
    isi = diff(spike_times);  % inter-spike interval
    burstinx = find(isi<0.01);   % ISI < 10 ms
    
    % Detect bursts
    % Burst: first ISI < 10 ms, subsequent ISIs < 15 ms
    bursts = {};
    used = [];
    burst1 = {};
    for k = 1:length(burstinx)
        if ismember(burstinx(k),used)
            continue   % already included in the previous burst
        end
        bursts{end+1} = spike_times([burstinx(k) burstinx(k)+1]); %#ok<AGROW>
        next = burstinx(k) + 1;
        if next > length(isi)
            burst1{end+1} = bursts{end}(1); %#ok<AGROW>
            continue   % last ISI of the cell
        end
        while isi(next) < 0.015  % if conseq. ISIs < 15 ms
            bursts{end} = [bursts{end}; spike_times(next+1)];
            used = [used next]; %#ok<AGROW>
            next = next + 1;
            if next > length(isi)
                break   % last ISI of the cell, quit while loop
            end
        end
        burst1{end+1} = bursts{end}(1); %#ok<AGROW>
    end
    bursts = vertcat(bursts{:});
    single = setdiff(spike_times, bursts);
    burst1 = cell2mat(burst1)';
    
    % Check output
    if ~isempty(spike_times)
        if ~isequal(spike_times, sort([bursts; single]))
            disp('Missing data during detection')
            keyboard;
        end
        if ~isequal(spike_times, union(bursts, single))
            disp('Missing data during detection')
            keyboard;
        end
    end
    
    % Output
    burst1Spikes{iC} = burst1;
    singleSpikes{iC} = single;
    burstAllSpikes{iC} = bursts;
end