function m5_imCohUsingFFT_EventPhase()
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
animal = animal_extract(codecorresfolder);

%% save setup
savefolder = codecorresfolder;


%%  input setup
inputfolder_rest = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\Jo\0_dataPrep\Rest\m5_imCohUsingFFT';
inputfolder_move = fullfile(codecorresParentfolder, 'm5_imCohUsingFFT_EventPhase');

image_type = 'tif';


unwanted_DBS = unwanted_DBS_extract(animal);
noisy_chns = noisy_chns_extract(animal);
removed_chns = [unwanted_DBS noisy_chns];
clear unwanted_DBS noisy_chns


EventPhases = {'preMove'; 'Anticipated';'earlyReach'; 'Return';'lateReach'};


fig_left = 50;
fig_bottom = 50;
fig_width = 1200;
fig_height = 600;

%% Code start here
cond_cell = cond_cell_extract(animal);

for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    % load iCOH_base
    iCOHfile_rest = fullfile(inputfolder_rest, ['Jo connectivity_' pdcond '.mat']);
    load(iCOHfile_rest, 'iCoh_pair', 'f_selected',  'chnPairNames');

    removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
    M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
    STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
    usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
    clear M1DBS_mask STN2GP_mask removedChns_mask

    iCoh_rest = abs(iCoh_pair(usedChnPairsMask, :));
    chnPairNames_rest = chnPairNames(usedChnPairsMask);
    f_selected_rest = round(f_selected,2);
    clear usedChnPairsMask usedChnPairsMask
    clear('iCoh_pair', 'f_selected',  'chnPairNames');
    
    for ei = 1: length(EventPhases)
        event = EventPhases{ei};
        
        
        if strcmpi(event, 'preMove')
            align2 = SKTEvent.TargetOnset;
            t_AOI = [-1 -0.8];
        end
        if strcmpi(event, 'Anticipated')
            align2 = SKTEvent.TargetOnset;
            t_AOI = [-0.2 0];
        end
        if strcmpi(event, 'earlyReach')
            align2 = SKTEvent.ReachOnset;
            t_AOI = [0 0.2];
        end
        if strcmpi(event, 'Return')
            align2 = SKTEvent.ReturnOnset;
            t_AOI = [0 0.2];
        end
        
        if strcmpi(event, 'lateReach')
            align2 = SKTEvent.Reach;
            t_AOI = [-0.2 0];
        end
        align2name = char(align2);
        if strcmpi(align2, 'Reach')
            align2name = 'Touch';
        end
        
        iCOHfile_move = fullfile(inputfolder_move, [animal '_' pdcond '_' event '_align2' align2name '.mat']);
        load(iCOHfile_move, 'iCoh_1time', 'f_selected',  'chnPairNames', 'ntrials');
        
        removedChns_mask = cellfun(@(x) contains(x, removed_chns), chnPairNames);
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
        usedChnPairsMask = (M1DBS_mask | STN2GP_mask) & ~removedChns_mask;
        clear M1DBS_mask STN2GP_mask removedChns_mask
        
        iCoh_move = abs(iCoh_1time(usedChnPairsMask, :));
        chnPairNames_move = chnPairNames(usedChnPairsMask);
        f_selected_move = round(f_selected,2);
        clear usedChnPairsMask usedChnPairsMask
               
        
        % calculate changes relative to Rest base
        if ~all(strcmp(chnPairNames_rest, chnPairNames_move))
            disp('chnPairNames_rest not equal chnPairNames_move' )
        else
            chnPairNames_show = chnPairNames_rest;
        end
        
        if ~all(f_selected_rest == f_selected_move)
            disp('f_selected_rest ~= f_selected_move')
        else
            f_selected = f_selected_rest;
        end
        digiCOH_move = iCoh_move;
        digiCOH_move(digiCOH_move>0) = 1;
        digiCOH_rest = iCoh_rest;
        digiCOH_rest(digiCOH_rest>0) = 1;
        
        diffiCOH = digiCOH_move - digiCOH_rest;
        
        
        % plot diffiCOH
        showData = diffiCOH;
        figure;
        set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
        imagesc(showData)
        pinky = [255,192,203 ]/ 255;
        red = [255 0 0] / 255;
        blue = [0 0 255]/ 255;
        cy = [0 255 255] / 255;
        gray = [211 211 211] / 255;
        cincrease = red;
        cnochange = gray;
        cdecrease = blue;
        cmap = [cdecrease; cnochange; cincrease]; colormap(cmap); colorbar
        set(gca,'CLim', [-1 1])
        
        set(gca, 'Position', [0.09 0.05 0.9 0.88])
        [npairs, nf] = size(showData);
        xticks([1:nf])
        xticklabels(round(f_selected,2))
        yticks([1:npairs]);
        set(gca,'YTickLabel',chnPairNames_show,'fontsize',12,'FontWeight','bold')
        xlabel('freqs')
        title([animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' 'align2 = ' align2name ', ntrials = ' num2str(ntrials)], ...
            'FontSize', 15, 'FontWeight', 'normal')
        set(gca,'CLim', [-1 1])
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
        
        % save image
        saveas(gcf, fullfile(savefolder, [animal '_diffciCOH_' event '_' pdcond '_align2' char(align2) '.' image_type]), image_type);
        close all
        
        
        clear('iCoh_1time', 'f_selected',  'chnPairNames', 'ntrials');
    end
    
end







