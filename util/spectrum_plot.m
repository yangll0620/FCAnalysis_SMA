function [power, f] = spectrum_plot(x, Fs)
% plot the spectrum of signal x
% 
% x: a vector
%
% Output:
%   power: single-sided spectrum P1 for each frequency, n_f * 1
%   f: each frequency, n_f *1



L = length(x1);
fft_x = fft(x);
P2 = abs(fft_x/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs *(0:(L/2))/L;
plot(f, P1)

power = P1;

