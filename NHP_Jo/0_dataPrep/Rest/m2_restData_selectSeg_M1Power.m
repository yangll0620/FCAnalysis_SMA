function m2_restData_selectSeg_M1Power()
%   Manually marked the good and bad segments
% 
%   Processing steps as follows:
%       1. add variable segsRemain, for manually marking each segment with 1 (good) or 0 (not good)
%           show the averaged lfp_m1 and each bipolar lfp_diffstn and lfp_diffgp
%
%       2. bipolar for STN and GP channels
%
%       3. Down sample trials into fs_new = 500


%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');


% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add util path
addpath(genpath(fullfile(codefolder,'util')));


% the corresponding pipeline folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

[datafolder, ~, ~, ~] = exp_subfolders();

%% global parameters
animal = 'Jo';

%% input setup
inputfolder = fullfile(codecorresParentfolder, 'm1_restData_cleaned_extract');

twin = 2;
toverlap = twin * 0.9;
freqs_roi = [0 50];

fs_new = 500;

%% save setup
savefolder = codecorresfolder;


%% Start Here

% segment check
files = dir(fullfile(inputfolder, ['*.mat']));
nfiles = length(files);
for fi = 1 : nfiles
    file = fullfile(files(fi).folder, files(fi).name);
    
    load(file, 'fs', 'data_segments',  'segsIndex')
    
    disp([num2str(fi) '/' num2str(nfiles) ' ' files(fi).name]);
    
    % nwin and noverlap for pwelch calculation
    nwin = round(twin * fs);
    noverlap = round(toverlap * fs);
    
    nsegs = length(data_segments);
    segsRemain = ones(1, nsegs);
    % pwelch for each segi
    for segi = 1: nsegs
        
        % lfp_m1: ntemp * nchns
        lfp_m1 = data_segments(segi).lfp_m1;
        lfp_meanm1 = mean(lfp_m1,2);
        
        [pxx_m1, F_m1] = pwelch(lfp_meanm1, nwin, noverlap, [freqs_roi(1):1/twin:freqs_roi(2)], fs);

        
        % lfp_stn: ntemp * 8
        lfp_stn = data_segments(segi).lfp_stn;
        lfp_diffstn = diff(lfp_stn,[],2);
        
        [pxx_stn, F_stn] = pwelch(lfp_diffstn, nwin, noverlap, [freqs_roi(1):1/twin:freqs_roi(2)], fs);
        
        
        % lfp_gp: ntemp * 8; lfp_diffgp: ntemp * 7
        lfp_gp = data_segments(segi).lfp_gp;
        lfp_diffgp = diff(lfp_gp,[],2);
        
        [pxx_gp, F_gp] = pwelch(lfp_diffgp, nwin, noverlap, [freqs_roi(1):1/twin:freqs_roi(2)], fs);
            
        % plot pxx_m1, pxx_stn and pxx_gp
        figure('units', 'normalized', 'outerposition', [0.1 0. 0.75 0.75])
        subplot(4,4,1); plot(F_m1, pxx_m1); legend('m1'); xlim(freqs_roi); 
        
        for chi = 1: size(lfp_diffstn, 2)
            subplot(4,4,chi + 1); plot(F_stn, pxx_stn(:, chi)); legend(['stn - ' num2str(chi)]);
            xlim(freqs_roi);
        end

        for chi = 1: size(lfp_diffstn, 2)
            subplot(4,4,chi + 8); plot(F_gp, pxx_gp(:, chi)); legend(['gp - ' num2str(chi)]);
            xlim(freqs_roi);
        end
        
        
        % manualy check the good('y') and the bad (n) segment
        reply = input(['segi=' num2str(segi) '/' num2str(nsegs) ', remain this seg? y/n [y]:'], 's');
        if isempty(reply)
            reply = 'y';
        end
        reply = lower(reply);
        if reply == 'n'
            segsRemain(segi) = 0;           
        end



        %%% bipolar stn and gp %%%%
        lfp_stn = lfp_diffstn;
        lfp_gp = lfp_diffgp;
        

        %%% down sample %%%
        lfp_m1 = resample(lfp_m1, round(fs_new), round(fs));
        lfp_stn = resample(lfp_stn, round(fs_new), round(fs));
        lfp_gp = resample(lfp_gp, round(fs_new), round(fs));
        
       
        data_segments(segi).lfp_m1 = lfp_m1;
        data_segments(segi).lfp_stn = lfp_stn;
        data_segments(segi).lfp_gp = lfp_gp;



        close all
        clear lfp_m1 lfp_meanm1 pxx_m1 F_m1
        clear lfp_stn lfp_diffstn pxx_stn F_stn lfp_gp lfp_diffgp pxx_gp F_gp
        clear reply
        
    end
    
    
    savefile = fullfile(savefolder, files(fi).name);
    save(savefile,'fs', 'data_segments',  'segsIndex','segsRemain');
    
    
    clear file fs data_segments
    clear nwin noverlap segi nsegs segsRemain
end


