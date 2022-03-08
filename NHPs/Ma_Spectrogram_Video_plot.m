function Ma_Spectrogram_Video_plot()
savecodefolder = '';
animal = 'Kitty';

folder_video = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\Animals\Kitty\KittyMinMaxPSD';
fname_video = 'Kitty20150408_2-1394 Desktop Video Camera.avi';
file_video = fullfile(folder_video, fname_video);

folder_timexls = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\Animals\Kitty\KittyMinMaxPSD';
fname_timexls = 'KittyTrialTime.xlsx';
file_timexls = fullfile(folder_timexls, fname_timexls);

folder_malfp = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT_SegV\m2_segSKTData_SelectTrials_goodReach';


% Figure parameters
fig_left = 50;
fig_bottom = 50;
fig_width = 500;
fig_height = 800;


%% Code start here
opts = detectImportOptions(file_timexls);
opts = setvaropts(opts, 'date', 'type','string');
tb_trialstime = readtable(file_timexls, opts);
clear opts


% extract dateofexp
tmp = regexp(fname_video, '[0-9]{8}_[0-9]{1}', 'match');
dateofexp = datenum(tmp{1}(1:8), 'yyyymmdd');
bk = str2num(tmp{1}(end));

% load malfp data
pdcond = parsePDCondition(dateofexp, animal); 
fname_malfp = [animal '_TrialsWMarkers_' pdcond '_' datestr(dateofexp, 'yyyymmdd') '_bktdt' num2str(bk) '.mat'];
file_malfp = fullfile(folder_malfp, fname_malfp);
dataloaded = load(file_malfp, 'lfpdata', 'smoothWspeed_trial', 'Wrist_smooth_trial', 'T_idxevent_ma', 'T_idxevent_lfp', 'T_chnsarea', 'fs_lfp', 'fs_ma');

% select dateofexp-bk trialtime
tb_times_1daybk = tb_trialstime(tb_trialstime.date == datestr(dateofexp, 'yyyymmdd') & tb_trialstime.bk == bk, :);

% read dateofexp-bk video
obj = VideoReader(file_video);
vidframes = read(obj);
vidFrate = obj.FrameRate;
clear obj


global lfpdata T_chnsarea fs_lfp T_idxevent_lfp
global smoothWspeed_trial  Wrist_smooth_trial T_idxevent_ma  fs_ma
global tri;



lfpdata = dataloaded.lfpdata;
T_chnsarea = dataloaded.T_chnsarea; 
fs_lfp = dataloaded.fs_lfp; 
T_idxevent_lfp = dataloaded.T_idxevent_lfp;
smoothWspeed_trial = dataloaded.smoothWspeed_trial;
Wrist_smooth_trial = dataloaded.Wrist_smooth_trial; 
T_idxevent_ma = dataloaded.T_idxevent_ma; 
fs_ma = dataloaded.fs_ma;
clear dataloaded


chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);
mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
T_chnsarea = T_chnsarea(mask_chnOfI, :);

%%% plot start here  %%%
fig = figure('Name',[pdcond ', ' datestr(dateofexp, 'yyyymmdd') '-bktdt' num2str(bk)]);
set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height], 'PaperPositionMode', 'auto');


tri = 1;
plot_play_1trial(tri)



function plot_1trial_MA_Spectrogram(fig, lfp_1trial, smoothWspeed_trial, Wrist_smooth_trial, T_idxevent_ma, T_idxevent_lfp, T_chnsarea, fs_lfp, fs_ma, tri)
align2 = SKTEvent.ReachOnset;
coli_align2 = uint32(align2);
coli_mouth = uint32(SKTEvent.Mouth);


% global parameters
twin = 0.2;
toverlap = 0.15;
f_AOI = [5 100];

pos_fig = get(fig, 'Position');
fig_width = pos_fig(3);
fig_height = pos_fig(4);

% spectrogram subplot parameters
margin_left = 0.15;
margin_top = 0.05;
margin_bottom = 0.05;
margin_right = 0.05;
margin_deltaX = 0.01;
margin_deltaY = 0.02;
clim = [-40 0];
areaName_left = 0.003;


% position of Next and previous button
btn_width = 50;
btn_height = 30;
pos_BtnNext_left = fig_width - btn_width;
pos_BtnNext_bottom = (1-margin_top/2) * fig_height - btn_height/2;
pos_BtnPrev_left = (margin_left * fig_width)/2 - btn_width/2;
pos_BtnPrev_bottom = pos_BtnNext_bottom;

eventline_colors = ['c', 'r', 'g', 'y', 'k'];
sktEvents = {'TargetOnset', 'ReachOnset', 'Touch', 'ReturnOnset', 'Mouth'};

%% Start Plot
clf(fig)

% draw btn_next and btn_prev
btn_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [pos_BtnNext_left pos_BtnNext_bottom btn_width btn_height]);
btn_next.Callback = @btn_nextTrial;
btn_prev = uicontrol(fig, 'Style','pushbutton', 'String', 'Previous', 'Position', [pos_BtnPrev_left pos_BtnPrev_bottom btn_width btn_height]);
btn_prev.Callback = @btn_prevTrial;

annotation(gcf,'textbox',...
    [0.75 1-margin_top 0.21 0.03],...
    'String', ['triali = ' num2str(tri)], ...
    'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');


% calc subp_height subp_width
nchns = size(lfp_1trial, 1);
nrows = nchns + 1;
ncols = 1;
subp_height = floor((1 - margin_top - margin_bottom - margin_deltaY * (nrows -1)) / nrows * 100)/100;
subp_width = floor((1 - margin_left - margin_right - margin_deltaX * (ncols -1)) / ncols * 100)/100;
clear nrows ncols



%%% plot spectrogram %%%
for chi = 1 : nchns
    x = lfp_1trial(chi, 1:T_idxevent_lfp{tri, coli_mouth});
    
    % spectrogram
    [~, freqs, times, psd] = spectrogram(x, round(twin * fs_lfp), round(toverlap * fs_lfp),[],fs_lfp); % psd: nf * nt
    
    % select freqs_plot and corresponding psd
    idx_f = (freqs >= f_AOI(1) & freqs <= f_AOI(2));
    freqs_plot =  freqs(idx_f);
    psd_plot = psd(idx_f, :);
    times_plot = times - T_idxevent_lfp{tri, coli_align2} / fs_lfp;
    
    % convert to db and gaussfilt
    psd_plot = 10 * log10(psd_plot);
    psd_plot = imgaussfilt(psd_plot,'FilterSize',5);
    
    % calculate subplot position
    subp_left = margin_left;
    rowi_subp = chi;
    subp_bottom = margin_bottom + (rowi_subp -1) * (subp_height + margin_deltaY);
    clear rowi_subp
    
    % actual plot
    subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
    imagesc(times_plot, freqs_plot, psd_plot, 'Tag', ['spectrogram-' num2str(tri) '-' num2str(chi)]);
    hold on
    
    % set gca
    set(gca,'YDir','normal')
    ylabel('Frequency(Hz)')
    % xtick
    if chi == 1
        xts = xticks(gca);
        xtls = xticklabels;
        xtls{cellfun(@(x) strcmp(x, '0'), xtls)} = char(align2);
        xticklabels(xtls);
        clear xts xlimts xtls
    else
        xticks([]);
    end
    if isempty(clim)
        set(gca,'YDir','normal')
    else
        set(gca,'YDir','normal', 'CLim', clim)
    end
    
    colormap(jet)
    colorbar
    annotation(gcf,'textbox',...
        [areaName_left subp_bottom+subp_height/2 0.03 0.03],...
        'String', T_chnsarea.brainarea{chi}, ...
        'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
    
    % plot event lines
    for ei = 1 : width(T_idxevent_lfp)
        t_event = (T_idxevent_lfp{tri, ei} - T_idxevent_lfp{tri, coli_align2}) / fs_lfp;
        plot([t_event t_event], ylim, [eventline_colors(ei) '--'], 'LineWidth',1.5)
        clear t_event
    end
    
    clear x freqs times psd idx_f
    clear freqs_plot psd_plot times_plot psd_plot
end



%%% plot MA data %%%

% calc times_plot_ma and ma_plot(ma_WSpeed, ma_W_X/Y/Z) data
madata_1trial = smoothWspeed_trial{tri};
madata2_1trial = Wrist_smooth_trial{tri};
ma_WSpeed = madata_1trial(1:T_idxevent_ma{tri, coli_mouth});
ma_W_X =  squeeze(madata2_1trial(1:T_idxevent_ma{tri, coli_mouth}, 1));
ma_W_Y =  squeeze(madata2_1trial(1:T_idxevent_ma{tri, coli_mouth}, 2));
ma_W_Z =  squeeze(madata2_1trial(1:T_idxevent_ma{tri, coli_mouth}, 3));
times_plot_ma = ([1: length(ma_WSpeed)] - T_idxevent_ma{tri, coli_align2} )/ fs_ma;


% find one ax for spectrogram
ax_spectrogram = findobj('Tag', ['spectrogram-' num2str(tri) '-' num2str(1)]).Parent;

% calc ma subplot pos
pos_spectrogram = get(ax_spectrogram, 'Position');
subp_left_ma = pos_spectrogram(1);
rowi_subp = nchns + 1;
subp_bottom_ma = margin_bottom + (rowi_subp -1) * (subp_height + margin_deltaY);
subp_width_ma = pos_spectrogram(3);
subp_height_ma = pos_spectrogram(4);

% actual plot
ax_ma = subplot('Position', [subp_left_ma, subp_bottom_ma, subp_width_ma, subp_height_ma]);
plot(times_plot_ma, ma_WSpeed, 'b');
hold on
plot(times_plot_ma, ma_W_X, 'r', times_plot_ma, ma_W_Y, 'g', times_plot_ma, ma_W_Z, 'k');

% set gca
set(gca, 'XLim', get(ax_spectrogram, 'XLim'));
legend({'Wrist Speed','Wrist-X', 'Wrist-Y', 'Wrist-Z'},...
    'Position',[0.15 0.92 0.45 0.07], 'NumColumns',2);
xticks([]);


annotation(gcf,'textbox',...
    [areaName_left subp_bottom_ma+subp_height/2 0.03 0.03],...
    'String', 'MA', ...
    'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');

clear madata_1trial madata2_1trial
clear ax_spectrogram pos_spectrogram subp_left_ma subp_bottom_ma subp_width_ma subp_height_ma xlim_spectrogram
clear ma_WSpeed ma_W_X ma_W_Y ma_W_Z

% plot event lines
for ei = 1 : width(T_idxevent_lfp) -1
    t_event = (T_idxevent_lfp{tri, ei} - T_idxevent_lfp{tri, coli_align2}) / fs_lfp;
    plot([t_event t_event], ylim, [eventline_colors(ei) '--'], 'LineWidth',1.5, 'DisplayName',sktEvents{ei})
    clear t_event
end

end

function btn_nextTrial(src,event)
tri = tri + 1;
plot_play_1trial(tri);
end

function btn_prevTrial(src,event)
tri = tri - 1;
plot_play_1trial(tri);
end


function play_1trial_video(frames, fps)
    figs = findall(0, 'Type', 'figure', 'Name', 'Movie Player');
    if ~isempty(figs)
        close(figs)
    end
    implay(frames, fps);
end

function plot_play_1trial(tri)
    lfp_1trial = lfpdata{tri}(mask_chnOfI, :);
    tinVid_str = tb_times_1daybk.tStrInVideo(tb_times_1daybk.triali == tri);
    tinVid_end = tb_times_1daybk.tEndInVideo(tb_times_1daybk.triali == tri);
    vidframe_1trial = vidframes(:, :, :, tinVid_str*vidFrate: tinVid_end*vidFrate);
    
    plot_1trial_MA_Spectrogram(fig, lfp_1trial, smoothWspeed_trial, Wrist_smooth_trial, T_idxevent_ma, T_idxevent_lfp, T_chnsarea, fs_lfp, fs_ma, tri);
    play_1trial_video(vidframe_1trial, vidFrate);  
end

end