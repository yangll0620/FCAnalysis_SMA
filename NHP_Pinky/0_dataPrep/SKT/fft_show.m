function fft_show(lfp, fs)

figure;
freqs = [5 50];

yfft = fft(lfp);
nfft = length(lfp);
fs_fft = (0:nfft-1)*(fs/nfft);
power = abs(yfft).^2/nfft;
idx = find(fs_fft<=freqs(2) & fs_fft>=freqs(1));

f_show = fs_fft(idx);
power_show = power(idx);
power_show = smooth(power_show, 20);

plot(f_show, power_show);