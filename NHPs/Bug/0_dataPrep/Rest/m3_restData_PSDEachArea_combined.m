function m3_restData_PSDEachArea_combined()


%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code') - 1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder, 'util')));

% the corresponding pipeline folder for this code
[~, codecorresParentfolder] = code_corresfolder(codefilepath, false, false);


animal = 'Bug';

savefolder = fullfile(codecorresParentfolder, 'm3_restData_PSDEachArea_extract');

%% Code Start Here

%%% combine all figures into one %%%%

% empty figure
figure
set(gcf, 'PaperUnits', 'points',  'PaperPosition', [18, 180, 300 180]);
annotation(gcf,'textbox',...
    [0.3 0.7 0.35 0.15],...
    'String',{[animal ' Rest Data']},...
    'LineStyle','none',...
    'FontSize',15, 'FontWeight','bold',...
    'FitBoxToText','off');
saveas(gcf, fullfile(savefolder, 'text'), 'png')
img_text =  imread(fullfile(savefolder, 'text.png'));


close all

%----DBS in one figure ---%
brainarea = 'STN';
imgs_col1 = []; imgs_col2 = []; % two columns
for i = 1: 7
    img = imread(fullfile(savefolder,['psd_' brainarea '_ch' num2str(i) '.png']));
    
    if mod(i, 2) == 1
        imgs_col1 = cat(1, imgs_col1, img);
    else
        imgs_col2 = cat(1, imgs_col2, img);
    end
    
    clear img
end
imgs_col2 = cat(1, imgs_col2, img_text); % the last one in column2 is empty
imgs_STN = cat(2, imgs_col1, imgs_col2);
clear imgs_col1 imgs_col2

brainarea = 'GP';
imgs_col1 = []; imgs_col2 = []; % two columns
for i = 1: 7
    img = imread(fullfile(savefolder,['psd_' brainarea '_ch' num2str(i) '.png']));
    
    if mod(i, 2) == 1
        imgs_col1 = cat(1, imgs_col1, img);
    else
        imgs_col2 = cat(1, imgs_col2, img);
    end
    
    clear img
end
imgs_col2 = cat(1, imgs_col2, zeros(size(img_text)) + 255); % the last one in column2 is empty
imgs_GP= cat(2, imgs_col1, imgs_col2);


imgs_DBS = cat(2, imgs_STN, imgs_GP);
imwrite(imgs_DBS,  fullfile(savefolder, 'combinedDBS.png'));




%----SMA, SensoryMotor and Thalams in one figure ---%
imgs_tha = [];
thalamus = {'VA', 'VLo', 'VPLo'};
for i = 1 : length(thalamus)
    brainarea = thalamus{i};
    
    img = imread(fullfile(savefolder,['psd_' brainarea '.png']));

    imgs_tha = cat(1, imgs_tha, img);
    
    clear img brainarea
end

img_SMA = imread(fullfile(savefolder,'psd_SMA.png'));
img_SensorMotor = imread(fullfile(savefolder,'psd_SensoryMotor.png'));

imgs_SMASM = cat(1, img_SMA, img_SensorMotor, img_text);
clear img_SMA img_SensorMotor  

img_SMASMTha = cat(2, imgs_tha, imgs_SMASM);
imwrite(img_SMASMTha,  fullfile(savefolder, 'combinedM1SMATha.png'));