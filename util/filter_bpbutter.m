function [data_filtered]= filter_bpbutter(data,fre,Fs)
%  band pass filtered using butterworth filter
%
%  Inputs:
%	data: a vector
%   fre:  a vector with two values
%   default Fs = 512
%  
%  Output:
%      data_filtered: filtered data
%
if nargin < 3
    Fs = 512;
end


% the approximate order n
% normalized frequency band edges Wn
data = double(data);
resp = 'bandpassiir'; % Response type
designMethod = 'butter'; % Design method
PassbandFrequency1 = fre(1);
StopbandFrequency1 = fre(1)*0.9;
PassbandFrequency2 = fre(2);
StopbandFrequency2 = fre(2)*1.1;
StopbandAttenuation1 = 40;
StopbandAttenuation2 = 40;
PassbandRipple = 1;
d = designfilt(resp, ...       % Response type
       'StopbandFrequency1',StopbandFrequency1, ...    % Frequency constraints
       'PassbandFrequency1',PassbandFrequency1, ...
       'PassbandFrequency2',PassbandFrequency2, ...
       'StopbandFrequency2',StopbandFrequency2, ...
       'StopbandAttenuation1',StopbandAttenuation1, ...   % Magnitude constraints
       'PassbandRipple',PassbandRipple, ...
       'StopbandAttenuation2',StopbandAttenuation2, ...
       'DesignMethod',designMethod, ...      % Design method
       'MatchExactly','passband', ...   % Design method options
       'SampleRate',Fs);

data_filtered    = filtfilt(d,data);


