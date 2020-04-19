function m3_restData_Spectrogram()
%% Visualize the preprocessed rest data recorded with grayMatter using spectrogram
%

%% folders generate
% the full path and the name of code file without suffix
% codefilepath = '/Users/linglingyang/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code/NHP_Pinky/0_dataPrep/Rest/m3_restData_Spectrogram'
codefilepath = mfilename('fullpath');

% find the codefolder = '/home/lingling/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code'
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1); 
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% the corresponding pipeline folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables
% find animal = 'Pinky'
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup
% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm2_restData_preprocessing');

% showed frequency
freqs = [5 40];
% showed time duration, 2min
ts = [60 60*3];

% areas to be visulized
areas_vis = {'M1', 'lSMA', 'lVA', 'lVPLo','lCd', 'lVLo', 'rSMA', 'rMC','rVLo','rVPLo','rVA'};

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'spectrogram';

%% starting
spectrogramfile = fullfile(savefolder, [savefilename_addstr '.mat']);
if(~exist(spectrogramfile, 'file'))

    % calc and store spectral analysis
    files_mild = dir(fullfile(inputfolder, '*mild*.mat'));
    [spower_showall_mild, f_showed_mild, t_showed_mild] = spectrogram_calc(files_mild, areas_vis, freqs, ts);
    
    files_normal = dir(fullfile(inputfolder, '*normal*.mat'));
    [spower_showall_normal, f_showed_normal, t_showed_normal] = spectrogram_calc(files_normal, areas_vis, freqs, ts);

    % save
    save(spectrogramfile, 'spower_showall_mild', 'f_showed_mild', 't_showed_mild', ...
        'spower_showall_normal','f_showed_normal', 't_showed_normal');
    
    clear files_mild files_normal 
    clear power_showall_* f_showed_*
else
    logtag = false;
    load(spectrogramfile, 'spower_showall_mild', 'f_showed_mild', 't_showed_mild', ...
        'spower_showall_normal','f_showed_normal', 't_showed_normal');
end

% plot spectral
savefile_prefix = fullfile(savefolder, savefilename_addstr);
%spectral_plot(f_showed_normal, power_showall_normal, f_showed_mild, power_showall_mild, freqs, savefile_prefix, logtag);
end


function [spower_show_all, f_showed, t_showed] = spectrogram_calc(files, areas_vis, freqs, ts)

nfiles = length(files);

% load lfp data from all files
for filei = 1 : nfiles
    % load data
    filename = files(filei).name;
    filefolder = files(filei).folder;
    load(fullfile(filefolder, filename), 'fs', 'lfpdata','T_chnsarea');
    
    idx_temporal = round(ts * fs);
    lfpdata = lfpdata(idx_temporal(1):idx_temporal(2)-1, :);
    if(filei == 1)
        lfpdatas = lfpdata;
        
        % get all the areas
        areas = {};
        for j = 1 : height(T_chnsarea)
            area = T_chnsarea.brainarea{j};
            if sum(cellfun(@(x) strcmp(x,  area), areas)) == 0
                areas = [areas, area];
            end
            clear area
        end
    else
        lfpdatas = cat(3, lfpdatas, lfpdata);
    end
    
    clear idx_temporal lfpdata
end


% calc and store spectrogram in each area
lfpdata = mean(lfpdatas, 3); % ntemporal * nchns
t_win = 6;
n_win = fs * t_win; 
noverlap = 0.8 * n_win;
nfft = n_win;
spower_show_all = struct();
for areai = 1: length(areas_vis)
    area = areas_vis{areai};
    
    % area chns
    chns = T_chnsarea.chni(cellfun(@(x) strcmp(x, area),  T_chnsarea.brainarea));
    
    % mean lfp across an area in ts duration
    lfp = mean(lfpdata(:,chns),2);
    
    % spectrogram
    [~,sp_f,sp_t,sp_power,~,~] = spectrogram(lfp,n_win,noverlap,nfft,fs);
    idx_f = find(sp_f>=freqs(1) & sp_f<=freqs(2));
    sp_f = sp_f(idx_f);
    sp_power = sp_power(idx_f, :);
    
    % assigen power_show to the struct power_show_all variable 
    eval(['spower_show_all.' area '= sp_power;']);
    
    if ~exist('t_showed','var')
        t_showed = sp_t;
    else
        if ~isequal(sp_t, t_showed)
            disp([str '-' area ': t_show not equal']);
        end
    end
    
    if ~exist('f_showed','var')
        f_showed = sp_f;
    else
        if ~isequal(sp_f, f_showed)
            disp([str '-' area ': f_show not equal']);
        end
    end
    
    clear area chns lfp sp_f sp_t sp_power idx_f
end
end



function spectrogram_plot(spower_showall, f_showed, t_showed, cond, freqs, savefile_prefix)
% plot
areas_vis = fieldnames(spower_showall);

for i = 1: length(areas_vis)
    area = areas_vis{i};
    eval(['spower_show = spower_showall.' area ';']);
    
       
    % plot
    figure
    
    imagesc(f_showed, t_showed, spower_show)
    
    plot(f_showed_normal,power_show_normal, 'LineWidth',2)
    xlim(freqs)
    % Create textbox to mark the frequency with highest power
    [~, idx_normal]= max(power_show_normal);
    legend("fnormal_{max} = " + num2str(f_showed_normal(idx_normal)), 'FontSize', 14);
    clear idx_*
    title([' Spectrum in ' area]);
    ylabel('Power')
    
    subplot(2,1,2)
    plot(f_showed_mild,power_show_mild, 'r', 'LineWidth',2);
    % Create textbox to mark the frequency with highest power
    [~, idx_mild]= max(power_show_mild);
    xlim(freqs)
    legend("fmild_{max} = " + num2str(f_showed_mild(idx_mild)), 'FontSize', 14);
    clear idx_*
    
    % label and title
    xlabel('Frequency/Hz')
    ylabel('Power')
    
    
    % save
    if logtag == true
        savefile_prefix = [savefile_prefix 'log'];
    end
    saveas(gcf,[savefile_prefix '_' area '.png']);
    
    clear area power_show_normal power_show_mild
end

