function m2_restDataGM_visual()
%% Visualize the preprocessed rest data recorded with grayMatter
%

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% pipelinefolder
[~, ~, pipelinefolder, ~] = exp_subfolders();
% the corresponding pipeline folder for this code
correspipelinefolder = code_corresfolder(codefilepath, true, false);

%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup
% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(pipelinefolder, ['NHP_' animal], '0_dataPrep', 'm2_restDataGM_preprocessing');

% data being cut for the first t_cut and the last t_cut
t_cut = 60;

% areas to be visulized
areas_vis = {'lSMA', 'rSMA'};

%% save setup
savefolder = correspipelinefolder;

%% starting
d = dir(fullfile(inputfolder, '*.mat'));
for filei = 1 : length(d)
    filename = d(filei).name;
        
    load(fullfile(inputfolder, filename), 'fs', 'lfpdata','T_chnsarea');
    
   
    % get all the areas
    areas = {};
    for j = 1 : height(T_chnsarea)
        area = T_chnsarea.brainarea{j};
        if sum(cellfun(@(x) strcmp(x,  area), areas)) == 0
            areas = [areas, area];
        end
        clear area
    end
        
    
    % prefix text 
    idx1 = strfind(filename, 'Pinky_preprodRestdataGM_');
    idx2 = strfind(filename, '.mat');
    tmpstr = filename(idx1 + length('Pinky_preprodRestdataGM_') : idx2 -1);
    str = strrep(tmpstr,'_','-');
    
    idx_cut = round(t_cut * fs);
    
    % visual lfpdata in each area
    for areai = 1: length(areas_vis)
        area = areas_vis{areai};
        
        % area1
        chns = T_chnsarea.chni(cellfun(@(x) strcmp(x, area),  T_chnsarea.brainarea));
        figure;
        for i = 1: length(chns)
            subplot(2,2,i);
            plot(lfpdata(idx_cut:end -idx_cut,chns(i)))
        end
        % Create textbox
        annotation(gcf,'textbox', [0.28 0.95 0.42 0.05], ...
            'String',[str '--' area],'LineStyle','none')
        
        % save
        saveas(gcf,fullfile(savefolder, [str '_' area '.png']));
        
        clear area chns i
    end
    
    close all
    clear filename
    clear fs lfpdata T_chnsarea
    clear idx1 idx2 tmpstr str
    clear idx_cut i areas
end

