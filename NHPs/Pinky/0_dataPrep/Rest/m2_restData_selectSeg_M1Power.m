function m2_restData_selectSeg_M1Power()
%   Manually marked the good and bad segments
% 
% 1. add variable segsRemain, for marking each segment with 1 (good) or 0 (not good)
%
% 2. add variable GMChnAreas from GMChnsarea file
%
%

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


%% input setup
inputfolder = fullfile(codecorresParentfolder, 'm1_restData_cleaned_extract');

twin = 2;
toverlap = twin * 0.9;
freqs_roi = [5 40];


% GrayMatter chn-area information file
filename_GMChnsarea = ['Pinky_GMChnAreaInf.csv'];
file_GMChnsarea =  fullfile(datafolder, filename_GMChnsarea);


%% save setup
savefolder = codecorresfolder;

%% Start Here

% GM chn Areas 
GMChnAreas = GMChnArea_extract(file_GMChnsarea);

% segment check
files = dir(fullfile(inputfolder, ['*moderate*.mat']));
nfiles = length(files);
for fi = 1 : nfiles
    file = fullfile(files(fi).folder, files(fi).name);
    
    load(file, 'fs', 'data_segments', 'chans_m1', 'segsIndex')
    
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
        
        [pxx_m1, F] = pwelch(lfp_meanm1, nwin, noverlap, [freqs_roi(1):1/twin:freqs_roi(2)], fs);
        
        
        % lfp_stn: ntemp * 8
        lfp_stn = data_segments(segi).lfp_stn;
        lfp_diffstn = diff(lfp_stn,[],2);
        
        [pxx_stn, F] = pwelch(lfp_diffstn(:,6), nwin, noverlap, [freqs_roi(1):1/twin:freqs_roi(2)], fs);
        
        
        % lfp_gp: ntemp * 8
        lfp_gp = data_segments(segi).lfp_gp;
        lfp_diffgp = diff(lfp_gp,[],2);
        
        [pxx_gp, F] = pwelch(lfp_diffgp(:,7), nwin, noverlap, [freqs_roi(1):1/twin:freqs_roi(2)], fs);
            
        %
        figure
        subplot(3,1,1);plot(F, pxx_m1)
        title([files(fi).name ', seg' num2str(segi)])
        
        subplot(3,1,2); plot(F, pxx_stn)
        subplot(3,1,3); plot(F, pxx_gp)
        
        
        
        
        % manualy check the good('y') and the bad (n) segment
        reply = input(['segi=' num2str(segi) '/' num2str(nsegs) ', remain this seg? y/n [y]:'], 's');
        if isempty(reply)
            reply = 'y';
        end
        reply = lower(reply);
        if reply == 'n'
            segsRemain(segi) = 0;           
        end
        
        close all
        clear lfp_m1 lfp_mean pxx reply
        
    end
    
    
    savefile = fullfile(savefolder, files(fi).name);
    save(savefile,'fs', 'data_segments', 'chans_m1', 'segsIndex','segsRemain', 'GMChnAreas');
    
    
    clear file chans_m1 fs data_segments
    clear nwin noverlap segi nsegs segsRemain
end


function GMChnAreas = GMChnArea_extract(file_GMChnsarea)
T = readtable(file_GMChnsarea);
chi_firstGM = 101;
nareas = height(T);
GMChnAreas = cell(nareas,1);
for areai = 1 : nareas 
    area = T.brainarea{areai};
    tmpcell = split(T.channels{areai}, ',');
    
    for j = 1 : length(tmpcell)
        chn = str2num(char(tmpcell{j}));
        GMChnAreas(chn-chi_firstGM+1,1) = {area};
    end
end

