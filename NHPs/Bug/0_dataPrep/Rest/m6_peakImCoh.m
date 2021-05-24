input_folder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Bug\0_dataPrep\Rest\m5_imCohUsingFFT';

animal = 'Bug';

patterns = {'ro', 'g+', 'b*'};
cond_cell = cond_cell_extract(animal);
figure
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    load(fullfile(input_folder, [animal ' connectivity_' pdcond '.mat']), 'chnPairNames', 'iCoh_pair');
    if ~exist('f_selected_show', 'var')
        load(fullfile(input_folder, [animal ' connectivity_' pdcond '.mat']), 'f_selected');
        f_selected_show = f_selected;
    else
        load(fullfile(input_folder, [animal ' connectivity_' pdcond '.mat']), 'f_selected');
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
    showData = iCoh_pair(usedChnPairsMask, :);
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
    plot(idxs, pairs, patterns{ci})
    ylim([1, npairs])
    hold on
    
    clear M1DBS_mask STN2GP_mask usedChnPairsMask
    clear chnPairNames f_selected iCoh_pair
end
yticks([1:npairs])
yticklabels(chnPairNames_show)
set(gca, 'YDir','reverse')
xticks([1: nf])
xticklabels(f_selected_show )
legend(cond_cell)


