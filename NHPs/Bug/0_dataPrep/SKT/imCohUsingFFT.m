clear
inputfolder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Bug\0_dataPrep\SKT\m1_SKTData_avgArea';
cond = 'normal';
files = dir(fullfile(inputfolder, ['*_' cond '_*.mat']));

% extract lfptrials align2 Event
t_minmax_reach = [0.6, 0.8];
t_minmax_return = [0.5, 0.9];
align2 = SKTEvent.ReachOnset;
tdur_trial = [-0.6 1];

[lfptrials, fs, T_chnsarea] = lfp_align2(files, align2, tdur_trial, t_minmax_reach, t_minmax_return);


% 
twin = 0.2;
toverlap = 0.15;
f_AOI = [8 50];

nwin = round(twin * fs);
noverlap = round(toverlap * fs);

[nchns, ~, ntrials] = size(lfptrials);
for chni = 1 : nchns-1
    lfptrialsi = squeeze(lfptrials(chni, :, :));
    for chnj = chni : nchns
        lfptrialsj = squeeze(lfptrials(chnj, :, :));
        
        cross_density_sum = 0;
        for triali = 1: ntrials
            x = lfptrialsi(: , triali);
            y = lfptrialsj(: , triali);
            
            [Sx, fx, tx, ~] = spectrogram(x, nwin, noverlap,[],fs); % Sx: nf * nt
            [Sy, fy, ty, ~] = spectrogram(y, nwin, noverlap,[],fs); % Sx: nf * nt
            
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
        cross_densitys(chni, chnj, :, :) = cross_density_sum / ntrials; % cross_density: nf * nt
        cross_densitys(chnj, chni, :, :) = cross_densitys(chni, chnj, :, :);
        
        
        clear cross_density_sum lfptrialsj
    end

    clear lfptrialsi
end

% imaginary of Coherency
iCoh = abs(cross_densitys);


% plot
imagesc(iCoh(:,:, 1,1));
set(gca, 'YDir', 'normal');
xticks([1:nchns]);
xticklabels(T_chnsarea.brainarea)
yticks([1:nchns]);
yticklabels(T_chnsarea.brainarea)


savefolder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Bug\0_dataPrep\SKT\test';
for chni = 1 : nchns -1 
        for chnj = chni + 1  : nchns
            chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];

        end
    end
for ni = 1 : size(iCoh, 4)
    nf = size(iCoh, 3);
    chnPairNames = {};
    iCoh_1time = zeros(nchns * (nchns -1)/2, nf);
    ci = 0;
    for chni = 1 : nchns -1 
        for chnj = chni + 1  : nchns
            chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];

            ci = ci + 1;
            iCoh_1time(ci, :) = iCoh(chni, chnj, :, ni);
        end
    end

    M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);

    chnPairsMask = M1DBS_mask;
    showData = iCoh_1time(chnPairsMask, :);
    
    figure('WindowState','maximized');
    imagesc(showData)
    [npairs, nf] = size(showData);
    xticks([1:nf])
    xticklabels(f_selected)
    yticks([1:npairs]);
    yticklabels(chnPairNames(chnPairsMask));
    xlabel('freqs')
    title(['timei = ' num2str(times(ni) + tdur_trial(1))  's'])
    set(gca,'CLim', [0 1])
    colorbar
    
    savefile = fullfile(savefolder, [num2str(ni) '.png']);
    saveas(gcf, savefile, 'png');
    close all
    
end




