function m6_FC_show()
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

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm5_imCohUsingFFT_EventPhase');


f_AOIs = [[23 28];[10 16]];


event = 'Reach';
align2 = 'ReachOnset';


%% code Start Here

cond_cell = cond_cell_extract(animal);
for ci = 1 : length(cond_cell)
    pdcond = cond_cell{ci};
    
    if strcmpi(event, 'preMove')
        if strcmpi(animal, 'Kitty') && strcmpi(pdcond, 'normal')
            align2 = 'ReachOnset';
        else
            align2 = 'TargetOnset';
        end
    end
    
    file = [animal '_' pdcond '_' event '_align2' align2 '.mat'];
    load(fullfile(inputfolder, file));
    
    for fi = 1 : size(f_AOIs, 1)
        f_AOI = f_AOIs(fi, :);
        
        idx_f = (f_selected >= f_AOI(1) & f_selected <= f_AOI(2));
        f_used = f_selected(idx_f);
        iCoh_used = iCoh_1time(:, idx_f);
        
        
        iCohDig_used = zeros(size(iCoh_used));
        iCohDig_used(iCoh_used >0 ) = 1;
        
        % majority selected
        avg_iCohDig = mean(iCohDig_used, 2);
        fConn = zeros(size(avg_iCohDig));
        fConn(avg_iCohDig > 0.5) =1;
        
        
        unwanted_DBS = unwanted_DBS_extract(animal);
        noisy_chns = noisy_chns_extract(animal);
        removed_chns = [unwanted_DBS; noisy_chns];
        clear unwanted_DBS noisy_chns
        
        % remove removed_chns and extract inter-cluster
        mask_removedPairs = cellfun(@(x) contains(x, removed_chns), chnPairNames);
        M1DBS_mask = cellfun(@(x) contains(x, 'M1-stn') || contains(x, 'M1-gp'), chnPairNames);
        STN2GP_mask = cellfun(@(x) contains(x, 'stn') && contains(x, 'gp'), chnPairNames);
        usedChnPairsMask = M1DBS_mask | STN2GP_mask;
        
        chnPairNames_used = chnPairNames(~mask_removedPairs & usedChnPairsMask);
        fConn_used = fConn(~mask_removedPairs & usedChnPairsMask);
        clear mask_removedPairs M1DBS_mask STN2GP_mask usedChnPairsMask
        
        % extract chnNames
        chnNames = {};
        for cpi = 1: length(chnPairNames_used)
            chnPairN = chnPairNames_used{cpi};
            
            % M1-stn0-1 style
            M_STNMatch = regexp(chnPairN, 'M1-stn[0-9]{1}-[0-9]{1}', 'match');
            if ~isempty(M_STNMatch)
                if isempty(chnNames) || ~any(contains(chnNames, 'M1'))
                    chnNames = [chnNames; {'M1'}];
                end
                
                tmp = M_STNMatch{1};
                stnName = tmp(4:end);
                if isempty(chnNames) || ~any(contains(chnNames, stnName))
                    chnNames = [chnNames; {stnName}];
                end
                clear tmp stnName
            end
            clear M_STNMatch
            
            
            % M1-gp0-1 style
            M_GPMatch = regexp(chnPairN, 'M1-gp[0-9]{1}-[0-9]{1}', 'match');
            if ~isempty(M_GPMatch)
                if isempty(chnNames) || ~any(contains(chnNames, 'M1'))
                    chnNames = [chnNames; {'M1'}];
                end
                
                tmp = M_GPMatch{1};
                gpName = tmp(4:end);
                if isempty(chnNames) || ~any(contains(chnNames, gpName))
                    chnNames = [chnNames; {gpName}];
                end
                clear tmp stnName
            end
            clear M_STNMatch
            
            
            clear chnPairN
        end
        
        
        % extract paired symmetric fConn_Matrix
        nchns = length(chnNames);
        fConn_Matrix = zeros(nchns, nchns);
        for i = 1 : nchns -1
            chnNamei = chnNames{i};
            for j = i +1 : nchns
                chnNamej = chnNames{j};
                
                % skip same cluster
                if (contains(chnNamei, 'stn') && contains(chnNamej, 'stn')) || (contains(chnNamei, 'gp') && contains(chnNamej, 'gp'))
                    continue;
                end
                
                idx = find(contains(chnPairNames_used, [chnNamei '-' chnNamej]));
                if ~isempty(idx)
                    fConn_Matrix(i, j) = fConn_used(idx);
                    fConn_Matrix(j, i) = fConn_Matrix(i, j);
                else
                    disp([chnNamei '-' chnNamej])
                end
                
                clear chnNamej
            end
            clear chnNamei
        end
        
        
        %%% plot FC
        fig = figure;
        set(fig, 'PaperUnits', 'points',  'Position', [700 600 400 350]);
        axes(fig, 'Position', [0.05 0.05 0.85 0.92])
        % plot M1
        pos_M1 = [-30, 0];
        
        chnName = 'M1';
        pos = pos_M1;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        hold on
        
        % stn0-1, 1-2 and 2-3
        pos_stn01 = [0, 50];
        pos_stn12 = [0, 45];
        pos_stn23 = [0, 40];
        
        chnName = 'stn0-1';
        pos = pos_stn01;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        chnName = 'stn1-2';
        pos = pos_stn12;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        chnName = 'stn2-3';
        pos = pos_stn23;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        
        % gp1-2, 2-3, 3-4, 4-5, 5-6, 6-7
        pos_gp12 = [30, 40];
        pos_gp23 = [30, 35];
        pos_gp34 = [30, 30];
        pos_gp45 = [30, 25];
        pos_gp56 = [30, 20];
        pos_gp67 = [30, 15];
        
        
        chnName = 'gp1-2';
        pos = pos_gp12;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        chnName = 'gp2-3';
        pos = pos_gp23;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        chnName = 'gp3-4';
        pos = pos_gp34;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        chnName = 'gp4-5';
        pos = pos_gp45;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        
        chnName = 'gp5-6';
        pos = pos_gp56;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        
        chnName = 'gp6-7';
        pos = pos_gp67;
        plot(pos(1), pos(2), 'ko', 'MarkerFaceColor', 'r', 'Tag', chnName);
        text(pos(1) + 1,pos(2), chnName)
        
        set(gca,'visible','off')
        
        % plot connection
        for i = 1 : length(chnPairNames_used)
            chnPairN = chnPairNames_used{i};
            fcValue = fConn_used(i);
            
            % M1-stn0-1 style
            M_STNMatch = regexp(chnPairN, 'M1-stn[0-9]{1}-[0-9]{1}', 'match');
            if ~isempty(M_STNMatch)
                tmp = M_STNMatch{1};
                stnName = tmp(4:end);
                
                oM1 = findobj(gcf, 'Tag', 'M1');
                oSTN = findobj(gcf, 'Tag', stnName);
                
                if fcValue > 0
                    plot([oM1.XData oSTN.XData], [oM1.YData oSTN.YData], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.2)
                end
                
                clear tmp stnName oM1 oSTN
                
            end
            clear M_STNMatch
            
            
            % M1-gp0-1 style
            M_GPMatch = regexp(chnPairN, 'M1-gp[0-9]{1}-[0-9]{1}', 'match');
            if ~isempty(M_GPMatch)
                tmp = M_GPMatch{1};
                gpName = tmp(4:end);
                
                oM1 = findobj(gcf, 'Tag', 'M1');
                oGP = findobj(gcf, 'Tag', gpName);
                
                if fcValue > 0
                    plot([oM1.XData oGP.XData], [oM1.YData oGP.YData], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.2)
                end
                
                clear tmp gpName oM1 oGP
                
            end
            clear M_GPMatch
            
            
            
            % stn0-1-gp0-1 style
            STN_GPMatch = regexp(chnPairN, 'stn[0-9]{1}-[0-9]{1}-gp[0-9]{1}-[0-9]{1}', 'match');
            if ~isempty(STN_GPMatch)
                tmp = STN_GPMatch{1};
                
                stnName = tmp(1:6);
                gpName = tmp(8:end);
                
                oSTN = findobj(gcf, 'Tag', stnName);
                oGP = findobj(gcf, 'Tag', gpName);
                
                if fcValue > 0
                    plot([oSTN.XData oGP.XData], [oSTN.YData oGP.YData], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.2)
                end
                
                clear tmp stnName gpName oSTN oGP
                
            end
            clear STN_GPMatch
        end
        
        text(0, 0, [animal ': ' event ' ' pdcond ' [' num2str(f_AOI(1)) ' ' num2str(f_AOI(2)) '] Hz'])
        
        savefile = [animal '-fc-' event ' ' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) '-'  pdcond '.png'];
        saveas(gcf, fullfile(savefolder, savefile));
        close gcf
    end
    
    
end








