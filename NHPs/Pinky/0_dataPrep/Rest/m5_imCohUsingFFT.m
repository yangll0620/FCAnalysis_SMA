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
image_type = 'bmp';

for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    savefile_prefix = fullfile(savefolder, [animal ' connectivity_' pdcond]);
    
    if ~exist([savefile_prefix, '.mat'])
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
                    f_selected =  round(fi(idx_f), 3);
                    
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

        % save data
    	save([savefile_prefix '.mat'], 'iCoh_pair', 'f_selected', 'chnPairNames')
    	
    else
        load([savefile_prefix '.mat'], 'iCoh_pair', 'f_selected',  'chnPairNames')
    end
    
    
    M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
    STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
    usedChnPairsMask = M1DBS_mask | STN2GP_mask;
    clear M1DBS_mask STN2GP_mask
    
    showData = abs(iCoh_pair(usedChnPairsMask, :));
    chnPairNames_show = chnPairNames(usedChnPairsMask);
    
    % plot connectivity
    figure;
    set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
    imagesc(showData)
    colormap(jet)
    set(gca, 'Position', [0.09 0.05 0.9 0.9])
    [npairs, nf] = size(showData);
    xticks([1:nf])
    xticklabels(round(f_selected,2) )
    yticks([1:npairs]);
    set(gca,'YTickLabel',chnPairNames_show,'fontsize',12,'FontWeight','bold')
    xlabel('freqs')
    title([ animal ' connectivity -'  pdcond ' in Rest'])
    set(gca,'CLim', [0 1])
    colorbar
    
    chnPair_prev = '';
    for ci = 1: length(chnPairNames_show)
        chnPair = chnPairNames_show{ci};
        
        % replace M1-stn0-1 to M1-STN
        s_stn = regexp(chnPair, 'stn[0-9]*-[0-9]*', 'match');
        if ~isempty(s_stn)
            for si = 1 : length(s_stn)
                chnPair = strrep(chnPair, s_stn{si}, 'STN');
            end
        end
        % replace M1-stn0-1 to M1-STN
        s_gp = regexp(chnPair, 'gp[0-9]*-[0-9]*', 'match');
        if ~isempty(s_gp)
            for si = 1 : length(s_gp)
                chnPair = strrep(chnPair, s_gp{si}, 'GP');
            end
        end
        
        if ~strcmp(chnPair_prev, '') && ~strcmp(chnPair_prev, chnPair) % a new site pairs
            hold on; plot(gca, xlim, [(ci + ci -1)/2 (ci + ci -1)/2], 'w--')
            % Create line
        end
        chnPair_prev = chnPair;
        
        clear s_stn s_gp chnPair
    end
    
    
    
    % save png
    saveas(gcf, [savefile_prefix '.' image_type], image_type);
    close all
    
    
    
    clear image
    clear pdcond files lfptrials fs T_chnsarea
    clear nchns ntrials iCoh tmpfolder video
end

% combined
images_combined = [];
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    try
        image = imread(fullfile(savefolder, [animal ' connectivity_' pdcond '.' image_type]));
        [m, n, c] =  size(image);
    catch
        if exist('m', 'var')
            image = zeros(m,n,c) + 255;
        end
    end
    
    images_combined = cat(2, images_combined, image);
    imwrite(images_combined,fullfile(savefolder, [animal ' connectivity_combined.' image_type]))
end



%%%-------  plot Peak ---------%
patterns = {'ko', 'g+', 'r*'};
cond_cell = cond_cell_extract(animal);
figure
set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    load(fullfile(savefolder, [animal ' connectivity_' pdcond '.mat']), 'chnPairNames', 'iCoh_pair');
    if ~exist('f_selected_show', 'var')
        load(fullfile(savefolder, [animal ' connectivity_' pdcond '.mat']), 'f_selected');
        f_selected_show = f_selected;
    else
        load(fullfile(savefolder, [animal ' connectivity_' pdcond '.mat']), 'f_selected');
        if any(~(f_selected_show == f_selected))
            disp([pdcond ' f_selected_show not equal f_selected'])
            continue;
        end
        clear f_selected
    end
    

    
    % select used chns
    M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
    STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
    usedChnPairsMask = M1DBS_mask | STN2GP_mask;
    showData = abs(iCoh_pair(usedChnPairsMask, :));
    if ~exist('chnPairNames_show', 'var')
        chnPairNames_show = chnPairNames(usedChnPairsMask);
    else
        chnPairNames_show_new = chnPairNames(usedChnPairsMask);
        cellcmp = strcmp(chnPairNames_show, chnPairNames_show_new);
        if any(~cellcmp)
            disp([pdcond ' chnPairNames_show_new not equal chnPairNames_show'])
            continue;
        end
        clear chnPairNames_show_new
    end
   
    if ~exist('npairs', 'var')
        [npairs, nf] = size(showData);
    else
        [m, n] = size(showData);
        if m ~= npairs || n ~= nf 
            disp([pdcond ' m ~= npairs or n ~= nf'])
            continue;
        end
        clear m n
    end
    
    [iCoh_Peak, idxs] = max(abs(showData), [], 2);
    idxs = idxs(find(iCoh_Peak ~= 0));
    pairs = find(iCoh_Peak ~= 0);
    
    % plot
    plot(f_selected_show(idxs), pairs, patterns{ci}, 'DisplayName',pdcond)
    ylim([1, npairs])
    xlim([f_selected_show(1), f_selected_show(end)])
    hold on
    
    
    
    if ci == length(cond_cell)
        chnPair_prev = '';
        for ci = 1: length(chnPairNames_show)
            chnPair = chnPairNames_show{ci};
            
            % replace M1-stn0-1 to M1-STN
            s_stn = regexp(chnPair, 'stn[0-9]*-[0-9]*', 'match');
            if ~isempty(s_stn)
                for si = 1 : length(s_stn)
                    chnPair = strrep(chnPair, s_stn{si}, 'STN');
                end
            end
            % replace M1-stn0-1 to M1-STN
            s_gp = regexp(chnPair, 'gp[0-9]*-[0-9]*', 'match');
            if ~isempty(s_gp)
                for si = 1 : length(s_gp)
                    chnPair = strrep(chnPair, s_gp{si}, 'GP');
                end
            end
            
            if ~strcmp(chnPair_prev, '') && ~strcmp(chnPair_prev, chnPair) % a new site pairs
                hold on; 
                plot(gca, xlim, [(ci + ci -1)/2 (ci + ci -1)/2], 'k--')
                % Create line
            end
            chnPair_prev = chnPair;
            
            clear s_stn s_gp chnPair
        end
    end
    
    
    clear M1DBS_mask STN2GP_mask usedChnPairsMask
    clear chnPairNames f_selected iCoh_pair
end
yticks([1:npairs])
set(gca,'YTickLabel',chnPairNames_show,'fontsize',12,'FontWeight','bold','Position', [0.09 0.05 0.9 0.9])
set(gca, 'YDir','reverse')
xticks(f_selected_show)
xticklabels(round(f_selected_show,2) )
set(gca, 'Color', 'w');
hleg = legend('show', 'Color', 'w');
hleg.String(end-1:end) = [];
title([ animal 'Peak connectivity in Rest'])

% save png
saveas(gcf, fullfile(savefolder,[animal '_peakFC' '.' image_type]), image_type);
close all


end