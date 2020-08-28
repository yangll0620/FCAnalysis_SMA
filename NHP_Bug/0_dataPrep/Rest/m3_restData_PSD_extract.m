function m3_restData_PSD_extract()
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
brainarea = 'm1';


% pwelch psd estimate variable
twin_pwelch = 2;


% variables for plotting
plotF_AOI = [5 50];


% input folder 
inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_M1Power');


%% save setup
savefolder = codecorresfolder;
savefilename_prefix = 'psd_';

%% Codes Start Here
file_psdallsegs = fullfile(savefolder, [savefilename_prefix 'allsegs_normalmild.mat']);

% calculate/load psd for mild and normal
if ~exist(file_psdallsegs, 'file') % not exist

    files_mild = dir(fullfile(inputfolder, '*_mild_*.mat'));
    [psd_allsegs_mild, F_pxx_mild] = psd_allsegs_fromfiles(files_mild, twin_pwelch, brainarea);


    files_normal = dir(fullfile(inputfolder, '*_normal_*.mat'));
    [psd_allsegs_normal, F_pxx_normal] = psd_allsegs_fromfiles(files_normal, twin_pwelch, brainarea);

    if ~isequal(F_pxx_mild, F_pxx_normal)
        disp('F_pxx_normal not equal F_pxx_mild');
    else
        F_pxx = F_pxx_normal;

        save(file_psdallsegs, 'psd_allsegs_normal', 'psd_allsegs_mild', 'F_pxx');
    end
else
    load(file_psdallsegs, 'psd_allsegs_normal', 'psd_allsegs_mild', 'F_pxx')
end


% plot PSD
savefigure_prefix =  fullfile(savefolder, [savefilename_prefix brainarea]);
plotPSD_comp(psd_allsegs_normal, psd_allsegs_mild, F_pxx, plotF_AOI, savefigure_prefix, brainarea);


function plotPSD_comp(psd_allsegs_normal, psd_allsegs_mild, F_all, plotF_AOI, savefile_prefix, brainarea)
%%  plot the psd comparison of normal and mild
%
%   Inputs:
%       psd_allsegs_normal, psd_allsegs_mild: psd of all segments in normal or mild, nfs1 * nsegs
%
%       F_all: the vector of frequencies (in hertz) at which the PSD is estimated (nfs * 1)
%       
%       plotF_AOI: the plot frequence range (in hertz)  (1 * 2)
%
%       savefile_prefix: the savefile prefix including folder


% find the idx for F_AOI
idx = find(F_all>= plotF_AOI(1) & F_all <= plotF_AOI(2));

% extract the F_AOI and psd of AOI
F_AOI = F_all(idx);
psd_normal_FAOI = psd_allsegs_normal(idx, :);
psd_mild_FAOI = psd_allsegs_mild(idx, :);


psd_normal_high = max(psd_normal_FAOI, [], 2);
psd_normal_low = min(psd_normal_FAOI, [], 2);
psd_normal_mean = mean(psd_normal_FAOI, 2);


psd_mild_high = max(psd_mild_FAOI, [], 2);
psd_mild_low = min(psd_mild_FAOI, [], 2);
psd_mild_mean = mean(psd_mild_FAOI, 2);




% reshape F_all, psd_*_high/low/mean into 1 * nfs
[m, n] = size(F_AOI); F_AOI = reshape(F_AOI, 1, m *n); clear m n
    
[m, n] = size(psd_normal_high); psd_normal_high = reshape(psd_normal_high, 1, m *n);    clear m n
[m, n] = size(psd_normal_low);  psd_normal_low = reshape(psd_normal_low, 1, m *n);      clear m n
[m, n] = size(psd_normal_mean); psd_normal_mean = reshape(psd_normal_mean, 1, m *n);    clear m n


[m, n] = size(psd_mild_high);   psd_mild_high = reshape(psd_mild_high, 1, m *n);    clear m n
[m, n] = size(psd_mild_low);    psd_mild_low = reshape(psd_mild_low, 1, m *n);      clear m n
[m, n] = size(psd_mild_mean);   psd_mild_mean = reshape(psd_mild_mean, 1, m *n);    clear m n




% colors setup
color_normal_range = [224,255,255]/255;
color_normal_mean = [0,0,255]/255;
color_mild_range = [255,228,225]/255;
color_mild_mean = [255,00,0]/255;

% plot range 
figure
fill([F_AOI flip(F_AOI)],[psd_normal_high flip(psd_normal_low)],color_normal_range,'LineStyle','none')
hold all
fill([F_AOI flip(F_AOI)],[psd_mild_high flip(psd_mild_low)],color_mild_range,'LineStyle','none')

% plot mean
h1 = plot(F_AOI, psd_normal_mean, 'Color',color_normal_mean);
h2 = plot(F_AOI, psd_mild_mean, 'Color',color_mild_mean);


% find the frequency with maximum density
[maxPSD, idx]= max(psd_mild_mean);
F_maxPSD = round(F_AOI(idx));
plot([F_maxPSD F_maxPSD], [0 maxPSD + maxPSD * 0.2], 'k--')


% legend
legend([h1, h2], {'normal', 'mild'})

% title
title(['PSD in ' upper(brainarea) ', center freq = ' num2str(F_maxPSD)])


% save figure
saveas(gcf, [savefile_prefix '_compNormalMild'], 'png')





