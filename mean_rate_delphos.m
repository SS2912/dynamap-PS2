
function [mean_window, events] = mean_rate_delphos(results, windowsize, event)

%% Analyze DELPHOS results with sliding window
% 
% Syntax:  
%    [mean_window, events] = mean_rate_delphos(results, windowsize)
%
% Inputs:
%   results    - matlab file obtained with DELPHOS
%   windowsize - length of the window in seconds
%   event      - string of the event ('Spike' or 'Gamma' etc)
% 
% Outputs:
%   mean_window - [Nch x Nw] Mean rate at each window and for each channel
%   spikes      - All the detections (timepoints) for each channel
%

% Authors: Victor LM (original window_spikes.m, modified by Sara S. Feb2023)


overlap = 0; % no overlap when calculating delphos 
spikes=cell(length(results.labels),1);
events=[];

for iter=1:length(results.markers)
    if strcmp(results.markers(iter).label, event) == 1
        for iter2=1:length(results.labels)
            if strcmp(results.markers(iter).channels, results.labels(iter2)) == 1
                spikes{iter2} = [spikes{iter2} results.markers(iter).position];
                events = [events results.markers(iter).position];
            end
        end
    end
end

step = round(windowsize*(1-overlap));
init_time=results.cfg.start;
Nw = round((results.cfg.duration-windowsize)/step+1);

for w=1:Nw
    init = init_time + (w-1)*step;
    endt = init + windowsize;
    
    for ch=1:length(results.labels)
        vector_of_spikesperwindow = spikes{ch}>init & spikes{ch}<endt;
%         mean_window(ch,w) = sum(vector_of_spikesperwindow)/windowsize;
        mean_window(ch,w) = sum(vector_of_spikesperwindow);

    end
end

