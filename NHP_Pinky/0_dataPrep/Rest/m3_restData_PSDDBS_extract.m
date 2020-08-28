function m3_restData_PSDDBS_extract()
%   PSD estimates for mild and normal individually


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


%%  input setup


% pwelch psd estimate variable
twin_pwelch = 2;


% variables for plotting
plotF_AOI = [5 50];


% input folder 
inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_M1Power');


%% save setup
savefolder = codecorresfolder;
savefilename_prefix = 'psd_DBS';

%% Codes Start Here
file_psdallsegs = fullfile(savefolder, [savefilename_prefix '_allsegs_normalmild.mat']);
% calculate/load psd for mild and normal
if ~exist(file_psdallsegs, 'file') % not exist

    files_mild = dir(fullfile(inputfolder, '*_mild_*.mat'));
    files_normal = dir(fullfile(inputfolder, '*_normal_*.mat'));
    
    [psd_stn_allsegs_normal, F_stn_pxx_normal] = psd_allsegs_fromfiles(files_normal, twin_pwelch, 'stn');
    [psd_stn_allsegs_mild, F_stn_pxx_mild] = psd_allsegs_fromfiles(files_mild, twin_pwelch, 'stn');


    
    [psd_gp_allsegs_normal, F_gp_pxx_normal] = psd_allsegs_fromfiles(files_normal, twin_pwelch, 'gp');
    [psd_gp_allsegs_mild, F_gp_pxx_mild] = psd_allsegs_fromfiles(files_mild, twin_pwelch, 'gp');

    if ~isequal(F_stn_pxx_normal, F_stn_pxx_mild, F_gp_pxx_normal, F_gp_pxx_mild)
        disp('F_pxx_normal not equal F_pxx_mild');
    else
        F_pxx = F_stn_pxx_normal;

        save(file_psdallsegs, 'psd_stn_allsegs_normal', 'psd_stn_allsegs_mild', 'psd_gp_allsegs_normal', 'psd_gp_allsegs_mild','F_pxx');
    end
else
    load(file_psdallsegs, 'psd_stn_allsegs_normal', 'psd_stn_allsegs_mild', 'psd_gp_allsegs_normal', 'psd_gp_allsegs_mild','F_pxx');
end


% plot PSD
plotPSD_comp(psd_stn_allsegs_normal, psd_stn_allsegs_mild, F_pxx, plotF_AOI, savefolder, 'stn');
plotPSD_comp(psd_stn_allsegs_normal, psd_stn_allsegs_mild, F_pxx, plotF_AOI, savefolder, 'gp');


function plotPSD_comp(psd_allsegschns_normal, psd_allsegschns_mild, F_all, plotF_AOI, savefolder, brainarea)
%%  plot the psd comparison of normal and mild
%
%   Inputs:
%       psd_allsegs_normal, psd_allsegs_mild: psd of all segments in normal
%       or mild, nfs * nchns * nsegs
%
%       F_all: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
%
%       plotF_AOI: the plot frequence range (in hertz)  (1 * 2)
%
%       savefile_prefix: the savefile prefix including folder


% find the idx for F_AOI
idx_AOI = find(F_all>= plotF_AOI(1) & F_all <= plotF_AOI(2));


% colors setup
color_normal_range = [224,255,255]/255;
color_normal_mean = [0,0,255]/255;
color_mild_range = [255,228,225]/255;
color_mild_mean = [255,00,0]/255;

% F_AOI
F_AOI = F_all(idx_AOI);
% reshape F_AOI
[m, n] = size(F_AOI); F_AOI = reshape(F_AOI, 1, m *n); clear m n

nchns = size(psd_allsegschns_normal, 2);
for chni = 1: nchns
    
    psd_allsegs_normal = squeeze(psd_allsegschns_normal(:, chni, :));
    psd_allsegs_mild = squeeze(psd_allsegschns_mild(:, chni, :));
       
    
    % extract the  psd of AOI
    psd_normal_FAOI = psd_allsegs_normal(idx_AOI, :);
    psd_mild_FAOI = psd_allsegs_mild(idx_AOI, :);
    
    
    psd_normal_high = max(psd_normal_FAOI, [], 2);
    psd_normal_low = min(psd_normal_FAOI, [], 2);
    psd_normal_mean = mean(psd_normal_FAOI, 2);

    
    psd_mild_high = max(psd_mild_FAOI, [], 2);
    psd_mild_low = min(psd_mild_FAOI, [], 2);
    psd_mild_mean = mean(psd_mild_FAOI, 2);
    
    
    
    
    % reshape, psd_*_high/low/mean into 1 * nfs
    
    [m, n] = size(psd_normal_high); psd_normal_high = reshape(psd_normal_high, 1, m *n);    clear m n
    [m, n] = size(psd_normal_low);  psd_normal_low = reshape(psd_normal_low, 1, m *n);      clear m n
    [m, n] = size(psd_normal_mean); psd_normal_mean = reshape(psd_normal_mean, 1, m *n);    clear m n
    
    
    [m, n] = size(psd_mild_high);   psd_mild_high = reshape(psd_mild_high, 1, m *n);    clear m n
    [m, n] = size(psd_mild_low);    psd_mild_low = reshape(psd_mild_low, 1, m *n);      clear m n
    [m, n] = size(psd_mild_mean);   psd_mild_mean = reshape(psd_mild_mean, 1, m *n);    clear m n
   
    
    
    % plot range
    figure
    fill([F_AOI flip(F_AOI)],[psd_normal_high flip(psd_normal_low)],color_normal_range,'LineStyle','none')
    hold all
    fill([F_AOI flip(F_AOI)],[psd_mild_high flip(psd_mild_low)],color_mild_range,'LineStyle','none')
    
    % plot mean
    h1 = plot(F_AOI, psd_normal_mean, 'Color',color_normal_mean);
    h2 = plot(F_AOI, psd_mild_mean, 'Color',color_mild_mean);
    
    
    % find the frequency with maximum density
    [maxPSD, idx_max]= max(psd_mild_mean);
    F_maxPSD = round(F_AOI(idx_max));
    plot([F_maxPSD F_maxPSD], [0 maxPSD + maxPSD * 0.2], 'k--')
    
    
    
    xlim([min(F_AOI) max(F_AOI)])
    
    
    % legend
    legend([h1, h2], {'normal', 'mild'})
    
    % title
    title(['PSD in ' upper(brainarea) ', chi = ' num2str(chni) ', high freq = ' num2str(F_maxPSD)])
    
    
    % save figure
    savename = fullfile(savefolder, ['psd_' brainarea '_ch' num2str(chni)]);
    saveas(gcf, savename, 'png')
    
        
    clear psd_allsegs_normal psd_allsegs_mild
    clear psd_normal_FAOI psd_mild_FAOI
    clear psd_normal_high psd_normal_low psd_normal_mean psd_mild_high psd_mild_low psd_mild_mean
    clear h1 h2 maxPSD F_maxPSD idx_max
    clear savename
end


function [psd_allsegs, F_pxx] = psd_allsegs_fromfiles(files, twin_pwelch, brainarea)
%% extract psd of all the segments from all files
%
% Outputs:
%   
%       psd_allsegs: the PSD estimate of all the segments from the files  (nfs * nsegs)
%
%       F_pxx: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
%
%       brainarea: the brain area to analysis ('stn', or 'gp')




nfiles = length(files);
psd_allsegs = [];
for filei = 1 : nfiles
    
    % load data
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'data_segments', 'fs', 'segsRemain');
    
    
    for segi = 1: length(segsRemain)
        if segsRemain(segi) == 0 % ignore the segment marked with 0
            continue;
        end
        
        
        % extract the lfp data of segi in brainarea
        eval(['lfp_oneseg = data_segments(segi).lfp_' brainarea ';']);
        
        
        %%% calcualte PSD throuch zscore and pwelch%%%
        
        % zscore of the lfp_oneseg along each column
        lfp_zscore = zscore(lfp_oneseg);
        
        
        % psd pwelch paramers
        nwins = round(twin_pwelch * fs);
        noverlap = round(nwins * 0.9);
        
        % PSD is computed independently for each column, Pxx: nfs * nchns, f: nfs * 1
        [Pxx, f]= pwelch(lfp_zscore, nwins, noverlap, nwins, fs);
        
        
        if ~exist('F_pxx', 'var')
            F_pxx = f;
        else
            if ~isequal(F_pxx, f)
                disp(['F_pxx not equal f for ' filename ', segi = ' num2str(segi)])
                continue;
            end
        end
        
        % psd_allsegs: nfs * nchns * nsegs
        psd_allsegs = cat(3, psd_allsegs, Pxx);
        
        clear lfp_oneseg lfp_zscore Pxx f
    end
    
    clear filename data_segments fs segsRemain
end

