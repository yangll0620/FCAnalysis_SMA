clear
folder = '/home/lingling/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/pipeline/NHP_Pinky/0_dataPrep/STKtrialextract_0/';







lfpM1 = lfptrial(1:96,:,:);
lfpM1_mean = squeeze(mean(lfpM1, 1));

ind = min(idxevent(:,5));
lfpM1_mean = mean(lfpM1_mean(254:ind,:),2);


%% spectrogram
data4spectrogram = lfpM1_mean;
[s,f,t,power] = spectrogram(data4spectrogram,256,200,[4:15],fs);
t = t - idxevent(1,1)/fs;
ind = find(t<0);
power_base = repmat(mean(power(:, ind),2), 1, length(t));
relpower = power - power_base;
imagesc(t, f, relpower')
set(gca,'YDir','normal') 

%% spectrum calculation part
data4fft = lfpM1_mean;

len = length(data4fft);

% Fourier transform of the signal. 
fft_lfpM1 = fft(data4fft);

% the two-sided spectrum P2
P2 = abs(fft_lfpM1/len);

% single-sided spectrum P1 
P1 = P2(1:round(len/2+1));
P1(2:end-1) = 2*P1(2:end-1);

% the frequency domain f
f = fs*(0:round((len/2)))/len;

%% plot
fshow_range = [4 40];
ind_show = find(f>fshow_range(1) & f< fshow_range(2));
f_show = f(ind_show);
P1_show = P1(ind_show);

plot(f_show, P1_show) 
title('Single-Sided Amplitude Spectrum of M1 lfp')
xlabel('f (Hz)')
ylabel('|P1(f)|')
