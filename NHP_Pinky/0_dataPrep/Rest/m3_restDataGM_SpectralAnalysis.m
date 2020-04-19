function m2_restDataGM_SpectralAnalysis()
%% Visualize the preprocessed rest data recorded with grayMatter
%

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder = '/home/lingling/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
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
% find animal = 'Pinky'
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup
% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(pipelinefolder, ['NHP_' animal], '0_dataPrep', 'm2_restDataGM_preprocessing');

% showed frequency
freqs = [5 40];
% showed time duration
ts = [60 60*3];

% areas to be visulized
areas_vis = {'lSMA', 'rSMA'};

%% save setup
savefolder = correspipelinefolder;

%% starting
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['Spectral Analyzing']);
for filei = 1 : nfiles
    % wait bar
    waitbar(filei/nfiles,f,['Spectral Analyzing lfp data in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data
    filename = files(filei).name; 
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
    
    % spectral analysis in each area
    aidx_temporal = round(ts * fs);
    for areai = 1: length(areas_vis)
        area = areas_vis{areai};
        
        % area chns
        chns = T_chnsarea.chni(cellfun(@(x) strcmp(x, area),  T_chnsarea.brainarea));
        % mean lfp across an area in ts duration
        lfp = mean(lfpdata(idx_temporal(1):idx_temporal(2),chns),2);
        % fft
        yfft = fft(lfp);
        
        % plot the power spectrum       
        nfiles = length(lfp);
        f = (0:nfiles-1)*(fs/nfiles);
        power = abs(yfft).^2/nfiles; % power of the DFT
        % select f_show and power_show
        idx = find(f<=freqs(2) & f>=freqs(1));
        f_show = f(idx);        
        
        
        power_show = 10 * log10(power(idx));
        %power_show = power(idx);
        power_show = smooth(power_show, 1000);
        
        % plot
        figure
        plot(f_show,power_show)
        xlabel('Frequency')
        ylabel('Power')
        title([str '--' area] );
        
        % save
        saveas(gcf,fullfile(savefolder, [str '_' area '.png']));
        
        eval(['power_show_' area '= power_show;']);
        if ~exist('f_showed','var')
            f_showed = f_show;
        else
            if ~isequal(f_show, f_showed)
                disp([str '-' area ': f_show not equal']);
            end
        end
        clear area chns lfp n yfft f power idx f_show power_show
       
    end
    save(fullfile(savefolder, [str '.mat']), 'power_show*', 'f_showed')
    clear f_showed
    
    close all
end