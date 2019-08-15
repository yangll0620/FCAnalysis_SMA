function [data_filtered]= filter_lpbutter(data,fre,Fs)
%  band pass filtered using butterworth filter
%
%  Inputs:
%	data: a vector
%   fre:  a scalar value
%  
%  Output:
%      data_filtered: filtered data
%

% the approximate order n
% normalized frequency band edges Wn
PassbandFre = fre;
StopbandFre = fre * 1.1;
PassbandRipple = 1;
StopbandAttenuation = 40;
designMethod = 'butter'; % Design method
d = designfilt('lowpassiir', ...        % Response type
       'PassbandFrequency',PassbandFre, ...     % Frequency constraints
       'StopbandFrequency',StopbandFre, ...
       'PassbandRipple',PassbandRipple, ...          % Magnitude constraints
       'StopbandAttenuation',StopbandAttenuation, ...
       'DesignMethod',designMethod, ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',Fs);               % Sample rate

data_filtered    = filtfilt(d,data);


