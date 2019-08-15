% Need validation and used
% Need modify the following description and put it in the front of this script
% Need set an threshold for maxdis used for detecting touch, lower than the threshold be treated as failsure trial

% extract the trial
% 1. find the start point of kluver, using jointstartpoint_hist
% 
% 2. find the time point of target onset
%   1. find the peaks of distances between each sample and the start point
%   of kluver.
%   
%   2. the time point before the peak distances will be identifies as the
%   onset of target if the distance of this time point reaches the
%   dis_onset (dis_onset if a threshold should be decided ahead, now = 10)

%
% 3. find the time point of reach onset
%   3-1. the rough range of a trial: from target onset to the time point 1.5s before next peak.
%   3-2. 


% final step, check the segmented trial, criteria: (1) total trial length
% (2) reaction time (3)