function m3_restData_SpectralAnalysis()
%% Visualize the preprocessed rest data recorded with grayMatter
%

%% folders generate
% the full path and the name of code file without suffix
% codefilepath = '/Users/linglingyang/yang7003@umn.edu/NMRC_umn/Projects/FCAnalysis/exp/code/NHP_Pinky/0_dataPrep/Rest/m3_restData_SpectralAnalysis'
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
% showed time duration, 5min
ts = [60 60*5];

% areas to be visulized
areas_vis = {'M1', 'lSMA', 'lVA', 'lVPLo','lCd', 'lVLo', 'rSMA', 'rMC','rVLo','rVPLo','rVA'};

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'spectral_';

%% starting
spectralfile = fullfile(savefolder, [savefilename_addstr '_' num2str(freqs(1)) '-' num2str(freqs(2)) '.mat']);
if(~exist(spectralfile, 'file'))

    % calc and store spectral analysis
    files_mild = dir(fullfile(inputfolder, ['*mild*.mat']));
    [power_showall_mild, f_showed_mild] = spectralAnalysis(files_mild, areas_vis, freqs, ts);
    
    files_normal = dir(fullfile(inputfolder, ['*normal*.mat']));
    [power_showall_normal, f_showed_normal] = spectralAnalysis(files_normal, areas_vis, freqs, ts);

    % save
    save(spectralfile, 'power_showall_mild', 'f_showed_mild', 'power_showall_normal','f_showed_normal');
    
    clear files_mild files_normal 
else
    load(spectralfile,'power_showall_mild', 'f_showed_mild', 'power_showall_normal','f_showed_normal');
end

% plot spectral
savefile_prefix = fullfile(savefolder, savefilename_addstr);
logtag = false;
spectral_plot(f_showed_normal, power_showall_normal, f_showed_mild, power_showall_mild, freqs, savefile_prefix, logtag);
end



function [power_show_all, f_showed] = spectralAnalysis(files, areas_vis, freqs, ts)
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

% calc and store spectral analysis in each area
power_show_all = struct();
lfpdata = mean(lfpdatas, 3);
for areai = 1: length(areas_vis)
    area = areas_vis{areai};
    
    % area chns
    chns = T_chnsarea.chni(cellfun(@(x) strcmp(x, area),  T_chnsarea.brainarea));
    
    % mean lfp across an area in ts duration
    lfp = mean(lfpdata(:,chns),2);
    
    % fft
    yfft = fft(lfp);
    
    % plot the power spectrum
    nfiles = length(lfp);
    f = (0:nfiles-1)*(fs/nfiles);
    power = abs(yfft).^2/nfiles; % power of the DFT
    % select f_show and power_show
    idx = find(f<=freqs(2) & f>=freqs(1));
    f_show = f(idx);
    
    
    % power_show = 10 * log10(power(idx));
    power_show = power(idx);
    power_show = smooth(power_show, 2000);
    
    % assigen power_show to the struct power_show_all variable 
    eval(['power_show_all.' area '= power_show;']);
    
    if ~exist('f_showed','var')
        f_showed = f_show;
    else
        if ~isequal(f_show, f_showed)
            disp([str '-' area ': f_show not equal']);
        end
    end
    clear area chns lfp n yfft f power idx f_show power_show
    
end
end

function spectral_plot(f_showed_normal, power_showall_normal, f_showed_mild, power_showall_mild, freqs, savefile_prefix, logtag)
% plot
areas_vis = fieldnames(power_showall_mild);

for i = 1: length(areas_vis)
    area = areas_vis{i};
    eval(['power_show_normal = power_showall_normal.' area ';']);
    eval(['power_show_mild = power_showall_mild.' area ';']);
    
    
    %low pass
    power_show_normal = lowpass(power_show_normal,5,120);
    power_show_mild = lowpass(power_show_mild,5,120);
    
    if logtag == true
        power_show_normal = 10 * log(power_show_normal);
        power_show_mild = 10 * log(power_show_mild);
    end
       
    % plot
    figure
    subplot(2,1,1)
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
close all
end