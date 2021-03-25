function m5_imCohUsingFFT()
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


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables

% animal
if ismac
    % Code to run on Mac platform
elseif isunix
    % Code to run on Linux platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '/', '[A-Za-z]*']);
elseif ispc
    % Code to run on Windows platform
    
    [fi, j] = regexp(codecorresfolder, ['NHPs', '\\', '[A-Za-z]*']);
else
    disp('Platform not supported')
end
animal = codecorresfolder(fi + length('NHPs') + 1:j);

%% save setup
savefolder = codecorresfolder;


%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');

% extract lfptrials align2 Event
align2 = SKTEvent.ReachOnset;


twin = 0.2;
toverlap = 0.15;
f_AOI = [8 50];

cond_cell = cond_cell_extract(animal);
if strcmpi(animal, 'bug')
    t_minmax_reach_normal = [0.5, 1];
    t_minmax_return_normal = [0.5, 1];
    t_minmax_reach_mild = [0.5, 1];
    t_minmax_return_mild = [0.5, 1];
    
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
end
if strcmpi(animal, 'jo')
    t_minmax_reach_normal = [0.5, 0.8];
    t_minmax_return_normal = [0.4, 0.8];
    t_minmax_reach_mild = [0.6 1];
    t_minmax_return_mild = [0.8 1.3];
    t_minmax_reach_moderate = [0.6 1];
    t_minmax_return_moderate = [0.8 1.4];
    
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_mild = [-0.5 0.5];
    tdur_trial_moderate = [-0.5 0.5];
    
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    t_minmax_reach_normal = [0.5, 10];
    t_minmax_return_normal = [0.5, 10];
    t_minmax_reach_moderate = [0.5, 10];
    t_minmax_return_moderate = [0.5, 10];
    
    tdur_trial_normal = [-0.5 0.5];
    tdur_trial_moderate = [-0.5 0.5];
end

if strcmpi(animal, 'pinky')
    t_minmax_reach_normal = [0.5, 1];
    t_minmax_return_normal = [0.5, 1];
    t_minmax_reach_mild = [0.5 1];
    t_minmax_return_mild = [0.5 1];
    t_minmax_reach_moderate = [0.7 1.2];
    t_minmax_return_moderate = [0.8 1.2];
    
    tdur_trial_normal = [-0.6 1];
    tdur_trial_mild = [-0.6 1];
    tdur_trial_moderate = [-.6 1];
end

for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    if strcmp(pdcond, 'normal')
        t_minmax_reach = t_minmax_reach_normal;
        t_minmax_return = t_minmax_return_normal;
        tdur_trial = tdur_trial_normal;
    else
        if strcmp(pdcond, 'mild')
            t_minmax_reach = t_minmax_reach_mild;
            t_minmax_return = t_minmax_return_mild;
            tdur_trial = tdur_trial_mild;
        else
            if strcmp(pdcond, 'moderate')
                t_minmax_reach = t_minmax_reach_moderate;
                t_minmax_return = t_minmax_return_moderate;
                tdur_trial = tdur_trial_moderate;
            end
        end
        
    end
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    [lfptrials, fs, T_chnsarea] = lfp_goodTrials_align2(files, align2, tdur_trial, t_minmax_reach, t_minmax_return);
    
    [nchns, ~, ntrials] = size(lfptrials);
    
    % calculate imaginary of Coherency
    for chni = 1 : nchns-1
        lfptrialsi = squeeze(lfptrials(chni, :, :));
        for chnj = chni : nchns
            lfptrialsj = squeeze(lfptrials(chnj, :, :));
            
            cross_density_sum = 0;
            for triali = 1: ntrials
                x = lfptrialsi(: , triali);
                y = lfptrialsj(: , triali);
                
                [Sx, fx, tx, ~] = spectrogram(x, round(twin * fs), round(toverlap * fs),[],fs); % Sx: nf * nt
                [Sy, fy, ty, ~] = spectrogram(y, round(twin * fs), round(toverlap * fs),[],fs); % Sx: nf * nt
                
                if chni == 1 && chnj == chni && triali == 1
                    freqs = fx;
                    times = tx;
                    idx_f = (freqs>=f_AOI(1) &  freqs<=f_AOI(2));
                    f_selected =  freqs(idx_f);
                end
                
                phix = angle(Sx(idx_f, :));
                phiy = angle(Sy(idx_f, :));
                cross_density_sum = cross_density_sum + exp(1i * (phix - phiy));
                
                clear x y Sx fx tx Sy fy ty
                clear phix phiy
                
            end
            
            
            
            if chni == 1 && chnj == chni + 1
                [nf, nt] = size(cross_density_sum);
                cross_densitys = zeros(nchns, nchns, nf, nt);
                clear nf nt
            end
            cross_densitys(chni, chnj, :, :) = cross_density_sum / ntrials; % cross_densitys: nchns * nchns * nf  * nt
            cross_densitys(chnj, chni, :, :) = cross_densitys(chni, chnj, :, :);
            
            
            clear cross_density_sum lfptrialsj
        end
        
        clear lfptrialsi
    end
    iCoh = imag(cross_densitys);
    
    lfp1 = squeeze(lfptrials(1, :, :));
    lfp2 = squeeze(lfptrials(10, :, :));
    [mus, stds] = psedo_SKTLFP_Test(lfp1, lfp2, 100, fs, twin, toverlap, f_AOI);
    clear lfp1 lfp2
    
    % pvalues using permutation test
    [nchns, ~, nf, nt] = size(iCoh);
    pvals = zeros(size(iCoh));
    for fi = 1 : nf
        for ti = 1 : nt
            mu = mus(fi, ti);
            std = stds(fi, ti);
            pd = makedist('Normal','mu',mu,'sigma',std);
            
            for chni = 1: nchns -1
                for chnj = chni : nchns
                    x = iCoh(chni, chnj, fi, ti);
                    pvals(chni, chnj, fi, ti) = (1 - cdf(pd,abs(x)));
                    pvals(chnj, chni, fi, ti) = pvals(chni, chnj, fi, ti);
                    clear x
                end
            end
            
            clear mu std pd
        end
    end
    % Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR) 
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);
    
    % set values not significant as 0
    iCoh(h==0) = 0;
    
    % show and save iCoh images and a converted video
    tmpfolder = fullfile(savefolder, pdcond);
    if exist(tmpfolder, 'dir')
        rmdir(tmpfolder, 's');
    end
    mkdir(tmpfolder)
    video = VideoWriter(fullfile(savefolder, [pdcond '.avi']));
    video.FrameRate = 10;
    open(video);
    for ti = 1 : size(iCoh, 4)
        % generate chnPairNames, such as M1-stn0-1
        nf = size(iCoh, 3);
        chnPairNames = {};
        iCoh_1time = zeros(nchns * (nchns -1)/2, nf);
        ci = 0;
        for chni = 1 : nchns -1
            for chnj = chni + 1  : nchns
                chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
                
                ci = ci + 1;
                iCoh_1time(ci, :) = iCoh(chni, chnj, :, ti);
            end
        end
        
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
        
        usedChnPairsMask = M1DBS_mask | STN2GP_mask;
        showData = iCoh_1time(usedChnPairsMask, :);
        
        figure('WindowState','maximized');
        imagesc(showData)
        [npairs, nf] = size(showData);
        xticks([1:nf])
        xticklabels(f_selected)
        yticks([1:npairs]);
        yticklabels(chnPairNames(usedChnPairsMask));
        xlabel('freqs')
        title([animal '-'  pdcond ' timei = ' num2str(times(ti) + tdur_trial(1))  's' ', ntrials = ' num2str(ntrials)])
        set(gca,'CLim', [0 1])
        colorbar
        
        savefile = fullfile(tmpfolder, [animal '_' pdcond '_t' num2str(ti) '.png']);
        saveas(gcf, savefile, 'png');
        close all  
        
        img = imread(savefile);
        writeVideo(video,img); %write the image to file
    end
    close(video); %close the file

    
    clear pdcond files lfptrials fs T_chnsarea
    clear nchns ntrials iCoh tmpfolder video
end
end