function m0_rest_spectrum_SMA()
%% extract the spectrum in area SMA during state

%% extract the corresponding pipeline folder for this code

% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

[~, ~, pipelinefolder, ~] = exp_subfolders();

correspipelinefolder = code_corresfolder(codefilepath, true, false);


%% global parameters
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);
% using the temporal data in [1 4]min
t_str = 60; t_dur = 60 * 4;

%% Input setup
% load data folder
restdatafolder = fullfile(pipelinefolder, ['/NHP_' animal '/0_dataPrep/m1_restDataSMA_extract']);

%% save setup
savefolder = correspipelinefolder;

%5 load the averaged rest lfp across channels and files
fs_new = 500;
[lfplSMA_normal, lfprSMA_normal] = avg_restlfpSMA_extract(restdatafolder, 'normal', fs_new, t_str, t_dur);
[lfplSMA_mild, lfprSMA_mild] = avg_restlfpSMA_extract(restdatafolder, 'mild', fs_new, t_str, t_dur);

[p_show_lnormal, f_show_lnormal]= calpower(lfplSMA_normal, fs_new);
[p_show_rnormal, f_show_rnormal]= calpower(lfprSMA_normal, fs_new);
[p_show_lmild, f_show_lmild]= calpower(lfplSMA_mild, fs_new);
[p_show_rmild, f_show_rmild]= calpower(lfprSMA_mild, fs_new);


%% plot 
figure
plot(f_show_lnormal, p_show_lnormal)
hold on
plot(f_show_rnormal, p_show_rnormal)
title('Spectrum of SMA lfp of Rest in normal state')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('leftSMA', 'rightSMA')
saveas(gcf,fullfile(savefolder, ['rest_spectrum_normal.png']))

figure
plot(f_show_lmild, p_show_lmild)
hold on
plot(f_show_rmild, p_show_rmild)
title('Spectrum of SMA lfp of Rest in mild state')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('leftSMA', 'rightSMA')

%% save
saveas(gcf,fullfile(savefolder, ['rest_spectrum_mild.png']))

function [p_show, f_show]= calpower(data, fs)

%% spectrum calculation part
len = length(data);

% Fourier transform of the signal. 
fft_lfp = fft(data);

% the two-sided spectrum P2
P2 = abs(fft_lfp/len);

% single-sided spectrum P1 
P1 = P2(1:round(len/2+1));
P1(2:end-1) = 2*P1(2:end-1);

% the frequency domain f
f = fs*(0:round((len/2)))/len;

% plot related
fshow_range = [8 40];
ind_show = find(f>fshow_range(1) & f< fshow_range(2));
f_show = f(ind_show);
P1_show = P1(ind_show);
p_show = smooth(P1_show, 100, 'moving');




function [lfp_lSMA, lfp_rSMA] = avg_restlfpSMA_extract(restdatafolder, pdCond, fs_new, t_str, t_dur)
%% extract averaged rest lfp in area SMA for pdCond across all the files
%
%      Process:
%           1. averaged across channels of one file
%           2. downsample to fs_new
%           3. averaged across all the files
%
%      return:
%           lfp_lSMA, lfp_rSMA: averaged lfp in lSMA and rSMA
%               

ntemporal = length(fs_new*t_str:fs_new*t_str + t_dur * fs_new -1);

%% load lfpdata and fs
files = dir(fullfile(restdatafolder, ['*' pdCond '*']));
lfplSMAs = zeros(ntemporal, length(files));
lfprSMAs = zeros(ntemporal, length(files));
for i = 1: length(files)
    filename = files(i).name;
    load(fullfile(restdatafolder, filename), 'lfpdata_lSMA', 'lfpdata_rSMA', 'fs');
        
    % average across channels
    lfp_lSMA = mean(lfpdata_lSMA,2);
    lfp_rSMA = mean(lfpdata_rSMA,2);
    
    
    %% downsample to 500Hz
    n = round(fs/fs_new);
    lfp_lSMA = downsample(lfp_lSMA,n);
    lfp_rSMA = downsample(lfp_rSMA,n);
    
    lfplSMAs(:,i) =  lfp_lSMA(fs_new*t_str:fs_new*t_str + t_dur * fs_new -1,1);
    lfprSMAs(:,i) =  lfp_rSMA(fs_new*t_str:fs_new*t_str + t_dur * fs_new -1,1);
    
    clear filename lfpM1 n fs
end

lfp_lSMA = mean(lfplSMAs,2);
lfp_rSMA = mean(lfprSMAs,2);
