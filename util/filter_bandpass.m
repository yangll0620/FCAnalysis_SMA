function [datalpfilter]= filter_bandpass(data,fre,Fs)
% 	[datalpfilter]= filter_bandpass(data,fre,Fs)
%	data is vector, fre is a vector with two values
%   default Fs = 512

if nargin < 3
    Fs = 512;
end

w_bpf = fre/(Fs/2);

% the approximate order n
% normalized frequency band edges Wn
% F is a vector of band edge frequencies in Hz
F = [fre(1)*0.9 fre fre(2)*1.05];
[n,Wn,beta,ftype] = kaiserord(F,[0 1 0],[0.01 0.05 0.01],Fs); 
b_bpf = fir1(n,w_bpf); 

data = double(data);

datalpfilter    = filtfilt(b_bpf,1,data);


