function VisualMATDTEvents()
%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder, 'toolbox', 'TDTMatlabSDK')))


% codecorresfolder, codecorresParentfolder
[codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
animal = animal_extract(codecorresfolder);

%% save setup
savefolder = codecorresfolder;


%% input: date tdtbk
dateofexp = datenum('101314','mmddyy');
tdtbk = 2;

savefig_format = 'tif'; 

% folder
folder_datadatabase = fullfile(server_NHP_extract(animal), animal, 'Recording', 'Processed', 'DataDatabase');
onedayfolder = fullfile(folder_datadatabase, [animal '_' datestr(dateofexp, 'mmddyy')]);
datebkfolder = fullfile(onedayfolder, ['Block-' num2str(tdtbk)]);

folder_rawtdt = 'Z:\root\Animals\Kitty\Recording\Raw\rawTDT';
rawtdtpath = fullfile(folder_rawtdt, ['KittyArrayDBSv2_DT1_' datestr(dateofexp, 'mmddyy')], ['Block-' num2str(tdtbk)]);


%% code start here
COTfilePattern = fullfile(datebkfolder, [animal datestr(dateofexp, 'yyyymmdd') '*_MA_COT*.mat']);
files = dir(COTfilePattern);
if length(files) ~= 1
    disp('COT file num != 1')
    return
end

% load COT data used MA System
load(fullfile(datebkfolder, files(1).name), 'COTData');

tbl_EventTimeixFromMA = table(COTData.TargetTime', COTData.ReachTimeix, COTData.TouchTimeix, COTData.ReturnTimeix, COTData.MouthTimeix, ...
    'VariableNames',{'TargetTimeix','ReachTimeix', 'TouchTimeix', 'ReturnTimeix', 'MouthTimeix'}); 

% tdt raw data
tdt = TDTbin2mat(rawtdtpath);
tdt_stpd = int16(tdt.streams.StPd.data);
tdt_stpd(tdt_stpd > 0) = 1; 
fs_tdt_stpd = tdt.streams.StPd.fs;
tdt_eventcode = int16(tdt.streams.Para.data);
fs_tdt_eventcode = tdt.streams.Para.fs;

if(fs_tdt_stpd ~= fs_tdt_eventcode)
    disp([' fs_tdt_stpd ~= fs_tdt_eventcode'])
    return;
end

fs_tdt = fs_tdt_stpd;

% 4-bit event codes
tdt_4bits = [tdt_stpd; tdt_eventcode(3, :); tdt_eventcode(2, :); tdt_eventcode(1, :)];
code_startTrial = [1;1;1;0];
code_goCue = [1; 1; 0; 1];
code_touch = [0; 1; 1; 1];

% plot ma speed
color_starttrial = 'k';
color_gocue = 'r';
color_reachonset = 'g';
color_touch = 'b';
color_mouth = 'y';
linestyle_eventma = '--';
linestyle_eventtdt = ':';
linestyle_tdtcode = '-';

ma_smoothWspeed = COTData.smoothWspeed;
fs_ma = COTData.SR;
t_bef = 2;
t_aft = 2;
for evi = 2 : height(tbl_EventTimeixFromMA)
    ma_targettimeix = tbl_EventTimeixFromMA{evi, 'TargetTimeix'};
    ma_reachtimeix = tbl_EventTimeixFromMA{evi, 'ReachTimeix'};
    ma_touchtimeix = tbl_EventTimeixFromMA{evi, 'TouchTimeix'};
    ma_mouthtimeix = tbl_EventTimeixFromMA{evi, 'MouthTimeix'};
    
    
    idx_str = ma_targettimeix - round(fs_ma * t_bef);
    idx_end = ma_mouthtimeix + round(fs_ma * t_aft);
    ma = ma_smoothWspeed(idx_str : idx_end);
    
    
    t_target = ma_targettimeix / fs_ma;
    t_reach = ma_reachtimeix / fs_ma;
    t_touch = ma_touchtimeix / fs_ma;
    t_mouth = ma_mouthtimeix / fs_ma;
    
    
    
    
    figure;
    plot([idx_str : idx_end] / fs_ma, ma, 'DisplayName','speed')
    xlabel('time/s')
    hold on
    plot([t_target t_target], ylim,  [color_starttrial linestyle_eventma], 'LineWidth',1.5, 'DisplayName', 'start trial (ma)')
    plot([t_reach t_reach], ylim,  [color_reachonset linestyle_eventma], 'LineWidth',1.5, 'DisplayName', 'reach onset (ma)')
    plot([t_touch t_touch], ylim,  [color_touch linestyle_eventma], 'LineWidth',1.5, 'DisplayName', 'touch (ma)')
    plot([t_mouth t_mouth], ylim,  [color_mouth linestyle_eventma], 'LineWidth',1.5, 'DisplayName', 'mouth (ma)')
    
    
    % extract data from tdt
    idx_str_tdt = round((t_target - t_bef) * fs_tdt);
    idx_end_tdt = round((t_mouth + t_aft) * fs_tdt);
    tdt_4bits_trial = tdt_4bits(:, idx_str_tdt: idx_end_tdt);

    mask_14StartTrial = all(tdt_4bits_trial == code_startTrial);
    mask_13goCue = all(tdt_4bits_trial == code_goCue);
    mask_11touch = all(tdt_4bits_trial == code_touch);
    idx_tdtStartTrial = find(mask_14StartTrial, 1) + idx_str_tdt;
    idx_tdtgoCue = find(mask_13goCue, 1) + idx_str_tdt;
    idx_tdttouch = find(mask_11touch, 1) + idx_str_tdt;
    
    lim = ylim();
    % plot eventcode
    plot([idx_str_tdt: idx_end_tdt] / fs_tdt, tdt_4bits_trial(1, :) * 20 + (lim(1) + lim(2)) /3 + 100, ['k' linestyle_tdtcode], 'DisplayName', 'startpad (tdt)')
    plot([idx_str_tdt: idx_end_tdt] / fs_tdt, tdt_4bits_trial(2, :) * 20 + (lim(1) + lim(2)) /3 + 60, ['r' linestyle_tdtcode], 'DisplayName', 'code3 (tdt)')
    plot([idx_str_tdt: idx_end_tdt] / fs_tdt, tdt_4bits_trial(3, :) * 20 + (lim(1) + lim(2)) /3 + 30, ['c' linestyle_tdtcode], 'DisplayName', 'code2 (tdt)')
    plot([idx_str_tdt: idx_end_tdt] / fs_tdt, tdt_4bits_trial(4, :) * 20 + (lim(1) + lim(2)) /3, ['b' linestyle_tdtcode], 'DisplayName', 'code1 (tdt)')
    
    % plot event detected from tdt
    %plot([idx_tdtStartTrial  idx_tdtStartTrial] / fs_tdt, ylim,  [color_starttrial linestyle_eventtdt], 'LineWidth',1.5, 'DisplayName', 'target onset 14(tdt)')
    plot([idx_tdtgoCue  idx_tdtgoCue] / fs_tdt, ylim,  [color_gocue linestyle_eventtdt], 'LineWidth',1.5, 'DisplayName', 'goCue 13(tdt)')
    %plot([idx_tdttouch  idx_tdttouch] / fs_tdt, ylim,  [color_touch linestyle_eventtdt], 'LineWidth',1.5, 'DisplayName', 'touch 11(tdt)')
    
    xlim([idx_str  idx_end] / fs_ma)
    legend();
    
    % save
    savefile = fullfile(savefolder, [datestr(dateofexp, 'yyyymmdd') '_trial' num2str(evi)]);
    saveas(gcf, savefile, savefig_format);
    
    close all
    clear ma_targettimeix ma_reachtimeix ma_touchtimeix ma_mouthtimeix idx_str idx_end ma
    clear t_target t_reach t_touch t_mouth 
    clear idx_str_tdt idx_end_tdt tdt_4bits_trial  
    clear mask_14StartTrial mask_13goCue  mask_11touch idx_tdtStartTrial idx_tdtgoCue idx_tdttouch
    clear lim savefile
end
