inputfolder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Kitty\0_dataPrep\SKT\m1_SKTData_avgArea';


t_bef = -0.5;
t_end = 1;

%%% normal case
files = dir(fullfile(inputfolder, '*.mat'));
for fi = 1 : length(files)
    load(fullfile(files(fi).folder, files(fi).name),  'fs', 'lfpdata', 'T_chnsarea', 'T_idxevent');
    
    idx0 = T_idxevent.ReachTimeix;
    mask_M1 = strcmp('M1', T_chnsarea.brainarea);
    
    % extract lfp_phase (ntemp * ntrials) from current file
    ntrials = length(idx0);
    for tri = 1 : ntrials
        idx_str = round(idx0(tri) + t_bef * fs);
        idx_end = round(idx0(tri) + t_end * fs);
        lfp_phase_1trial = squeeze(lfpdata(mask_M1, idx_str:idx_end, tri));
        
        if tri == 1
            lfp_phase = zeros(length(lfp_phase_1trial), ntrials);
        end
        lfp_phase(:, tri) = lfp_phase_1trial;
        
        
        clear idx_str idx_end lfp_phase_1trial
    end
    
    
    tmp = strsplit(files(fi).name, '.');
    filename = strrep(tmp{1}, '_', '-');
    clear tmp
    
    % plot mean lfp_phase
    times = [1: size(lfp_phase, 1)]/fs + t_bef;
    plot(times, mean(lfp_phase, 2))
    title(filename)
    saveas(gcf, fullfile(inputfolder, filename), 'png');
    close all
    
    clear fs lfpdata T_chnsarea T_idxevent
    clear idx0 mask_M1 filename
    clear times lfp_phase
end