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
inputfolder = fullfile(codecorresParentfolder, 'm3_restData_rmChns_avgArea');

fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;


twin = 0.2;
toverlap = 0.15;
f_AOI = [8 50];

cond_cell = cond_cell_extract(animal);


for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
    
    % extract all segments
    lfpAllSegs = [];
    for fi = 1: length(files)
        load(fullfile(files(fi).folder, files(fi).name), 'data_segments');
        
        lfpAllSegs = [lfpAllSegs data_segments];
        clear data_segments
    end
    
    
    
    
    % calculate imaginary of Coherency
    load(fullfile(files(1).folder, files(1).name), 'fs', 'T_chnsarea');
    nsegs = length(lfpAllSegs);
    for segi = 1 : nsegs
        lfpdata = lfpAllSegs(segi).lfp;
        nchns = size(lfpdata, 1);

        for chni = 1 : nchns-1
            lfpi = squeeze(lfpdata(chni, :));
            [Si, fi, ~, ~] = spectrogram(lfpi, round(twin * fs), round(toverlap * fs),[],fs); % Si: nf * nt
            
            % extract idx_f, f_selected and iCoh_eachSeg at the first time
            if segi == 1 && chni == 1
                idx_f = (fi>=f_AOI(1) &  fi<=f_AOI(2));
                f_selected =  fi(idx_f);
                
                nf = length(f_selected);
                crossDensity_eachSeg = zeros(nchns, nchns, nf, nsegs);
                
                clear nf
            end
            
            phii = angle(Si(idx_f, :));
            
            for chnj = chni + 1 : nchns
                lfpj = squeeze(lfpdata(chnj, :));
                
                [Sj, ~, ~, ~] = spectrogram(lfpj, round(twin * fs), round(toverlap * fs),[],fs); % Sj: nf * nt
                phij = angle(Sj(idx_f, :));
                
                crossDensity_eachSeg(chni, chnj, :, segi) = mean(exp(-1i * (phii - phij)), 2);
                crossDensity_eachSeg(chnj, chni, :, segi) = crossDensity_eachSeg(chni, chnj, :, segi);
                clear lfpj Sj phij
            end
            
            clear lfpi Si fi  phii
        end
        
        clear lfpdata nchns
    end
    iCoh = imag(mean(crossDensity_eachSeg, 4)); % iCoh: nchns * nchns * nf
    
    
    lfp1 = lfpAllSegs(1).lfp(1,:);
    lfp2 = lfpAllSegs(1).lfp(10,:);
    [mus, stds] = psedo_restLFP_Test(lfp1, lfp2, 100, fs, twin, toverlap, f_AOI);
    clear lfp1 lfp2
    
    % pvalues using permutation test
    [nchns, ~, nf] = size(iCoh);
    pvals = zeros(size(iCoh));
    for fi = 1 : nf
        mu = mus(fi, 1);
        std = stds(fi, 1);
        pd = makedist('Normal','mu',mu,'sigma',std);
        
        for chni = 1: nchns -1
            for chnj = chni : nchns
                x = iCoh(chni, chnj, fi);
                pvals(chni, chnj, fi) = (1 - cdf(pd,abs(x)));
                pvals(chnj, chni, fi) = pvals(chni, chnj, fi);
                clear x
            end
        end
        
        clear mu std pd
        
    end
    % Benjamini & Hochberg (1995) procedure for controlling the false discovery rate (FDR)
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);
    
    % set values not significant as 0
    iCoh(h==0) = 0;
    
    % show and save iCoh images 
    % generate chnPairNames, such as M1-stn0-1
    nf = size(iCoh, 3);
    chnPairNames = {};
    iCoh_pair = zeros(nchns * (nchns -1)/2, nf);
    ci = 0;
    for chni = 1 : nchns -1
        for chnj = chni + 1  : nchns
            chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
            
            ci = ci + 1;
            iCoh_pair(ci, :) = iCoh(chni, chnj, :);
        end
    end
    
    M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
    STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
    
    usedChnPairsMask = M1DBS_mask | STN2GP_mask;
    
    showData = iCoh_pair(usedChnPairsMask, :);
    
    figure;
    set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
    imagesc(showData)
    set(gca, 'Position', [0.07 0.05 0.92 0.9])
    [npairs, nf] = size(showData);
    xticks([1:nf])
    xticklabels(round(f_selected,3))
    yticks([1:npairs]);
    yticklabels(chnPairNames(usedChnPairsMask));
    xlabel('freqs')
    title([ animal ' connectivity -'  pdcond ' in Rest'])
    set(gca,'CLim', [0 1])
    colorbar
    
    savefile = fullfile(savefolder, [animal ' connectivity_' pdcond '.png']);
    saveas(gcf, savefile, 'png');
    close all

    clear image
    clear pdcond files lfptrials fs T_chnsarea
    clear nchns ntrials iCoh tmpfolder video
end

% combined
images_combined = [];
conds = {'normal', 'mild', 'moderate'};
for ci = 1 : length(conds)
    pdcond = conds{ci};
    try
        image = imread(fullfile(savefolder, [animal ' connectivity_' pdcond '.png']));
        [m, n, c] =  size(image);
    catch
        if exist('m', 'var')
            image = ones(m,n,c);
        end
    end
    
    images_combined = cat(2, images_combined, image);
    imwrite(images_combined,fullfile(savefolder, [animal ' connectivity_combined.png']))
end


end