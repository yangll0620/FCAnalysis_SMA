function m3_SKTData_spectrogram()


%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1:strfind(codefilepath, 'code') + length('code') - 1);

% add util path
addpath(genpath(fullfile(codefolder, 'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% pipelinefolder
[correspipelinefolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% input setup
folder_input = fullfile(codecorresParentfolder, 'm2_STKData_seg');

%  animal
[i,j]= regexp(folder_input, 'NHPs/[A-Za-z]*');
animal = folder_input(i + length('NHPs/'):j);
clear i j


%% save setup
savefolder = correspipelinefolder;


%% 

% pdcond = 'normal';
% files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
% plot_spectrogram_allfiles(files, savefolder, animal, pdcond)
% 
% pdcond = 'mild';
% files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
% plot_spectrogram_allfiles(files, savefolder, animal, pdcond)
% 
% pdcond = 'moderate';
% files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
% plot_spectrogram_allfiles(files, savefolder, animal, pdcond)


combined_imgs(savefolder, animal);

function plot_spectrogram_allfiles(files, savefolder, animal, pdcond)
%%

% global parameters
twin = 0.2;
toverlap = 0.15;
f_AOI = [8 50];

% load fs and chnAreas for all files
load(fullfile(files(1).folder, files(1).name),  'fs', 'chnAreas');

nwin = round(twin * fs);
noverlap = round(toverlap * fs);
nfft = noverlap;


lfpdata_allfiles = [];
for filei = 1 : length(files)
    
    % load data
    load(fullfile(files(filei).folder, files(filei).name), 'lfpdata');
    
   
    lfpdata_allfiles = cat(3, lfpdata_allfiles, lfpdata);
    
    clear lfpdata
end


% plot for each brain area
for areai = 1 : length(chnAreas)
    
    brainarea = chnAreas{areai};
    chnMask =  strcmp(chnAreas, brainarea);
    
    
    % calculate psds: nf * nt * ntrials
    psds = [];
    for triali = 1: size(lfpdata_allfiles, 3)
        
        x = lfpdata_allfiles(:, chnMask, triali);
        
        [~, f, t, ps] = spectrogram(x, nwin, noverlap,[],fs); % ps: nf * nt
        
        psds = cat(3, psds, ps);
        
        clear x ps
    end
    psd = mean(psds, 3);
    
    idx_f = (f>=f_AOI(1) &  f<=f_AOI(2));
    f_selected =  f(idx_f);
    psd_selected = psd(idx_f, :);
    t_selected = t - 1;
    
    
    psd_selected = 10 * log10(psd_selected');
    % psd_selected = imgaussfilt(psd_selected,10);
    
    figure
    set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
    imagesc(t_selected, f_selected, psd_selected);
    colorbar
    set(gca,'YDir','normal')
    
    xlabel('time/s')
    ylabel('psd')
    
    xlim([-0.94 1.5])
    
    title([animal ' ' pdcond ': ' brainarea])
    
    
    savefile = fullfile(savefolder, [animal '_' pdcond '_'  brainarea]);
    saveas(gcf, savefile, 'tiff');
   
    
    clear brainarea chnMask psds f t idx_f *_selected 
    clear savefile
end
close all

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



