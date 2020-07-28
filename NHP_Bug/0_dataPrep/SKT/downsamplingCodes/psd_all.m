clear 
close all

%% global parameters
% input
inputfolder = '.';

% global parameters
reachIdx = [2 3];
returnIdx = [4 5];

brainarea = 'M1';


% phase parameters
tphase_bef = 0.2; % the time before event onset
tmin_phase = 0.5;
tmax_phase = 1;

% baseline parameters
tbase_bef = 0.5; % the time before target onset


% parameters for spectrogram
twin = 0.4;
toverlap = twin * 0.9;

% frequency range of intest
f_roi = [10 50];
t_roi = -0.1;



%%
%eventName = 'reach'; cond = 'mild'; psd_reach;
eventName = 'reach'; cond = 'normal';psd_reach;

%eventName = 'return'; cond = 'mild';psd_return;
eventName = 'return'; cond = 'normal';psd_return;


