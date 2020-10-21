function m2_SKTData_spectrogram()
%  extract lfp data respect to reachonset
% 
%  return:
%        lfptrials: nchns * ntemp * ntrials


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables
% animal
[fi, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
animal = codecorresfolder(fi + length('NHPs/'):j);


%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm1_SKTData_avgArea');

%% save setup
savefolder = codecorresfolder;

%% starting: narrow filter the lfp data of all the files
files_normal = dir(fullfile(inputfolder, '*_normal_*.mat'));
files_mild = dir(fullfile(inputfolder, '*_mild_*.mat'));
files_moderate = dir(fullfile(inputfolder, '*_moderate_*.mat'));


tmin_reach_normal = 0.4;
tmax_reach_normal = 1;
tdur_trial_normal = [-0.8 0.8];
[lfptrials_normal, fs, T_chnsarea] = lfp_align2_reachonset(files_normal, tdur_trial_normal, tmin_reach_normal, tmax_reach_normal);
plot_spectrogram_allfiles(lfptrials_normal, T_chnsarea, fs, savefolder, animal, 'normal')



tmin_reach_mild = 0.5;
tmax_reach_mild = 1;
tdur_trial_mild = [-0.8 1.5];
[lfptrials_mild, fs, T_chnsarea] = lfp_align2_reachonset(files_mild, tdur_trial_mild, tmin_reach_mild, tmax_reach_mild);
plot_spectrogram_allfiles(lfptrials_mild, T_chnsarea, fs, savefolder, animal, 'mild')


tmin_reach_moderate = 0.5;
tmax_reach_moderate = 1.2;
tdur_trial_moderate = [-0.8 1.5];
[lfptrials_moderate, fs, T_chnsarea] = lfp_align2_reachonset(files_moderate, tdur_trial_moderate, tmin_reach_moderate, tmax_reach_moderate);
plot_spectrogram_allfiles(lfptrials_moderate, T_chnsarea, fs, savefolder, animal, 'moderate')



combined_imgs(savefolder, animal)


function plot_spectrogram_allfiles(lfptrials, T_chnsarea, fs, savefolder, animal, pdcond)
%%

% global parameters
twin = 0.2;
toverlap = 0.15;
f_AOI = [8 50];


nwin = round(twin * fs);
noverlap = round(toverlap * fs);
nfft = noverlap;

idxSTN = find(strcmp(T_chnsarea.brainarea, 'STN'));
for i = 1 : length(idxSTN)
    T_chnsarea.brainarea{idxSTN(i)} = ['stn' num2str(i-1) '-' num2str(i)];
    
end

idxGP = find(strcmp(T_chnsarea.brainarea, 'GP'));
for i = 1 : length(idxGP)
    T_chnsarea.brainarea{idxGP(i)} = ['gp' num2str(i-1) '-' num2str(i)];
    
end

% plot for each brain area
for areai = 1 : height(T_chnsarea)
    
    brainarea = T_chnsarea.brainarea{areai};
    chnMask =  strcmp(T_chnsarea.brainarea, brainarea);
    
    
    % calculate psds: nf * nt * ntrials
    psds = [];
    for triali = 1: size(lfptrials, 3)
        x = lfptrials(chnMask, : , triali);
        
        [~, f, t, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    psd = mean(psds, 3);
    
    idx_f = (f>=f_AOI(1) &  f<=f_AOI(2));
    f_selected =  f(idx_f);
    psd_selected = psd(idx_f, :);
    t_selected = t - 1;
    
    
    psd_selected = 10 * log10(psd_selected);
    psd_selected = imgaussfilt(psd_selected,'FilterSize',5);
    
    figure
    set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
    imagesc(t_selected, f_selected, psd_selected);
    colorbar
    set(gca,'YDir','normal')
    
    xlabel('time/s')
    ylabel('psd')
    xticks([-0.5 0 0.5])
    xticklabels({-0.5, 'reachonset', 0.5})
    
    title([animal ' ' pdcond ': ' brainarea])
    
    
    savefile = fullfile(savefolder, [animal '_' pdcond '_'  brainarea]);
    saveas(gcf, savefile, 'tiff');
   
    
    clear brainarea chnMask psds f t idx_f *_selected 
    clear savefile
end



function [lfptrials, fs, T_chnsarea] = lfp_align2_reachonset(files, tdur_trial, tmin_reach, tmax_reach)

coli_reachonset = 2;
coli_reach = 3;
coli_mouth = 5;

load(fullfile(files(1).folder, files(1).name),  'fs', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent');
    
    if(height(T_idxevent) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    
    ntrials = size(lfpdata, 3);
    for tri = 1: ntrials
        
        t_reach = (T_idxevent{tri, coli_reach} - T_idxevent{tri, coli_reachonset}) / fs;
        t_dur = (T_idxevent{tri, coli_mouth} - T_idxevent{tri, coli_reachonset}) / fs;
        if t_reach < tmin_reach || t_reach > tmax_reach || t_dur < tdur_trial(2)
            continue
        end
            
        idxdur = round(tdur_trial * fs) + T_idxevent{tri, coli_reachonset};
        lfp_phase_1trial = lfpdata(:, idxdur(1):idxdur(2), tri);
           
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
    end
end


function combined_imgs(folder, animal)

% --- M1 ---%
brainarea = 'M1';

img_normal = imread(fullfile(folder,[animal '_normal_' brainarea '.tif']));
img_mild = imread(fullfile(folder,[animal '_mild_' brainarea '.tif']));
img_moderate = imread(fullfile(folder,[animal '_moderate_' brainarea '.tif']));
imgs_M1 = cat(1, img_normal, img_mild, img_moderate);
imwrite(imgs_M1,  fullfile(folder, ['combined_' animal brainarea '.png']));
clear img_normal img_mild img_moderate



% --- STN ---%
brainarea = 'stn';
for i = 1: 7
    
    img_normal = imread(fullfile(folder,[animal '_normal_' brainarea num2str(i-1) '-' num2str(i) '.tif']));
    img_mild = imread(fullfile(folder,[animal '_mild_' brainarea num2str(i-1) '-' num2str(i) '.tif']));
    img_moderate = imread(fullfile(folder,[animal '_moderate_' brainarea num2str(i-1) '-' num2str(i) '.tif']));
    imgs_3 = cat(1, img_normal, img_mild, img_moderate);
    
    if mod(i, 2) == 1
        imgs = imgs_3;

    else
        imgs = cat(2, imgs, imgs_3);
        
        imwrite(imgs,  fullfile(folder, ['combined_' animal upper(brainarea)  num2str(i/2) '.png']));
    end
    
    if i == 7
        
        imwrite(imgs,  fullfile(folder, ['combined_' animal upper(brainarea)  num2str(i/2 + 1) '.png']));
    end
    
    clear img_normal img_mild img_moderate
end



% --- GP ---%
brainarea = 'gp';
for i = 1: 7
    
    img_normal = imread(fullfile(folder,[animal '_normal_' brainarea num2str(i-1) '-' num2str(i) '.tif']));
    img_mild = imread(fullfile(folder,[animal '_mild_' brainarea num2str(i-1) '-' num2str(i) '.tif']));
    img_moderate = imread(fullfile(folder,[animal '_moderate_' brainarea num2str(i-1) '-' num2str(i) '.tif']));
    imgs_3 = cat(1, img_normal, img_mild, img_moderate);
    
    if mod(i, 2) == 1
        imgs = imgs_3;

    else
        imgs = cat(2, imgs, imgs_3);
        
        imwrite(imgs,  fullfile(folder, ['combined_' animal upper(brainarea)  num2str(i/2) '.png']));
    end
    
    if i == 7
        
        imwrite(imgs,  fullfile(folder, ['combined_' animal upper(brainarea)  num2str(i/2 + 1) '.png']));
    end
    
    clear img_normal img_mild img_moderate
end
