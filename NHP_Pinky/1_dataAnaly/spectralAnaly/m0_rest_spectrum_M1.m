function rest_spectrum_0()

%% codecorresfolder
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

codefolder = codefilepath(1: strfind(codefilepath, 'code') + length('code')-1);

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

[datafolder, ~, pipelinefolder, ~] = exp_subfolders();

codecorresfolder = code_corresfolder(codefilepath, true, false);

%% save folder
savefolder = codecorresfolder;

%% load data folder
restdatafolder = fullfile(pipelinefolder, '/NHP_Pinky/0_dataPrep/restDataextract_1');


% load rest data lfp
fs_new = 500;
lfpM1_normal = restlfpM1_extract(restdatafolder, 'normal', fs_new);
lfpM1_mild = restlfpM1_extract(restdatafolder, 'mild', fs_new);

[p_show_normal, f_show_normal]= calpower(lfpM1_normal, fs_new);
[p_show_mild, f_show_mild]= calpower(lfpM1_mild, fs_new);

plot(f_show_normal, p_show_normal) 
hold on
plot(f_show_mild, p_show_mild) 
title('Spectrum of M1 lfp of Rest')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend('normal', 'mild')

%% save
saveas(gcf,fullfile(savefolder, ['rest_spectrum.png']))

function [p_show, f_show]= calpower(data, fs)

%% spectrum calculation part
len = length(data);

% Fourier transform of the signal. 
fft_lfpM1 = fft(data);

% the two-sided spectrum P2
P2 = abs(fft_lfpM1/len);

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




function lfpM1 = restlfpM1_extract(restdatafolder, pdCond, fs_new)
% 
t_str = 60;
t_dur = 60*4;

ntemporal = length(fs_new*t_str:fs_new*t_str + t_dur * fs_new -1);

%% load lfpdata and fs
files = dir(fullfile(restdatafolder, ['*' pdCond '*']));
lfpM1s = zeros(ntemporal, length(files));
for i = 1: length(files)
    filename = files(i).name;
    load(fullfile(restdatafolder, filename));
    
    
    % average across channels
    lfpM1 = mean(lfpdata,2);
    
    
    %% downsample to 500Hz
    n = round(fs/fs_new);
    lfpM1 = downsample(lfpM1,n);
    
    lfpM1s(:,i) =  lfpM1(fs_new*t_str:fs_new*t_str + t_dur * fs_new -1,1);
    
    clear filename lfpM1 n fs
end

lfpM1 = mean(lfpM1s,2);
