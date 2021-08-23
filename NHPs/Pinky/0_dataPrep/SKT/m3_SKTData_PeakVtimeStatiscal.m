function m3_SKTData_PeakVtimeStatiscal()


%% extract the corresponding pipeline folder for this code
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');
% code folder
codefolder = codefilepath(1:strfind(codefilepath, 'code') + length('code') - 1);

% add util path
addpath(genpath(fullfile(codefolder, 'util')));

% add NexMatablFiles path
addpath(genpath(fullfile(codefolder, 'toolbox', 'NexMatlabFiles')))

% pipelinefolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% input setup: For Kitty normal and moderate in two folders
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');

%  animal 
animal = animal_extract(codecorresfolder);


tasks = {'reachonset2peakVs', 'peakV2reachs'};

%% save setup
savefolder = codecorresfolder;
savefilename = [animal '_peakV_timeStaComp'];




%% code Starting Here
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);


%%% extract t_reachonset2peakVs_normal/moderate,t_peakV2reachs_moderate %%%
cond_cell = cond_cell_extract(animal);
for ci = 1: length(cond_cell)
    pdcond = cond_cell{ci};
    files = dir(fullfile(inputfolder, ['*' pdcond '*.mat']));
    
    t_reachonset2peakVs = [];
    t_peakV2reachs = [];
    for fi = 1 : length(files)
        file = fullfile(files(fi).folder, files(fi).name);
        load(file, 'lfpdata', 'goodTrials', 'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma');

        [~, ~, ntrials] = size(lfpdata);
        for tri = 1: ntrials
            
            % ignore trials marked with 0
            if ~goodTrials(tri)
                continue
            end
            
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            
            
            [~, idx] = max(smoothWspeed_trial(idx_reachonset_ma: idx_reach_ma, tri));
            idx_peakV = idx + idx_reachonset_ma -1;
            clear idx
            
            t_reachonset2peakV = (idx_peakV - idx_reachonset_ma)/ fs_ma;
            t_peakV2reach = (idx_reach_ma - idx_peakV)/ fs_ma;
            
            t_reachonset2peakVs = [t_reachonset2peakVs; t_reachonset2peakV];
            t_peakV2reachs = [t_peakV2reachs; t_peakV2reach];
            
            clear idx_reachonset_ma idx_reach_ma idx_peakV
            clear t_reachonset2peakV t_peakV2reach
        end
        
        
        clear file 
        clear('lfpdata', 'goodTrials', 'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma');
        clear nchns ntrials tri
    end
    
    
    eval(['t_reachonset2peakVs_' pdcond ' = t_reachonset2peakVs;'])
    eval(['t_peakV2reachs_' pdcond ' = t_peakV2reachs;'])
    
    clear pdcond files t_reachonset2peakVs t_peakV2reachs
end



%%% plot pdcond comparison %%%
figure('Position',[675 549 570* 2 413]);

% position of each subplot
x0 = 0.05; y = 0.11;
w =  0.4; h = 0.815;  gap = 0.1;


alpha = 0.05;

% extract group vector gs
gs =[];
for ci = 1 : length(cond_cell)
    pdcond  = cond_cell{ci};
    
    eval(['n = length(t_peakV2reachs_' pdcond ');'])
    g = repmat({pdcond}, n,1);
    gs = [gs; g];
    clear pdcond
end


%----- time comparison ----------%
for ti = 1 : length(tasks)
    task = tasks{ti};
    position = [x0 + (w + gap) * (ti-1) y w h];
    
    
    % extract t_normal, t_mild and t_moderate
    ts = [];
    for ci = 1 : length(cond_cell)
        pdcond  = cond_cell{ci};
        eval(['t_' pdcond ' = t_' task '_' pdcond ';']);
        eval(['ts = [ts; t_' pdcond '];'])
        clear pdcond
    end
    
    % wilcoxon rank sum test
    if exist('t_normal', 'var') && exist('t_mild', 'var')
        p1 = ranksum(t_normal, t_mild);
    end
    if exist('t_mild', 'var') && exist('t_moderate', 'var')
        p2 = ranksum(t_mild, t_moderate);
    end
    if exist('t_normal', 'var') && exist('t_moderate', 'var')
        p3 = ranksum(t_normal, t_moderate);
    end
    

    % actual plot
    subplot('Position', position);
    boxplot(ts, gs); hold on
    title([animal ':' task ' time comparison'])
    
    
    % significant part
    xt = get(gca, 'XTick');
    yt = get(gca, 'YTick');
    axis([xlim    0  max(yt) * 1.1])
    
    if exist('p1', 'var') && exist('p2', 'var') && exist('p3', 'var')
        if p1 < alpha
            plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
        end
        if p2 < alpha
            plot(xt([2 3]), [1 1]*max(yt)*0.97, '-k',  mean(xt([2 3])), max(yt) * 0.99, '*k')
        end
        if p3 < alpha
            plot(xt([1 3]), [1 1]*max(yt)*1.02, '-k',  mean(xt([1 2])), max(yt) * 1.04, '*k')
        end
    else
        if exist('p1', 'var')
            if p1 < alpha
                plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
            end
        end
        
        if exist('p2', 'var')
            if p2 < alpha
                plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
            end
        end
        
        if exist('p3', 'var')
            if p3 < alpha
                plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
            end
        end
    end
     
    
    
    
    % Create textbox
    if ti == 1
        pos1 = [0.23 0.1 0.06 0.07];
        pos2 = [0.25 0 0.2 0.2];
    else
        if ti == 2
            pos1 = [0.73 0.1 0.06 0.07];
            pos2 = [0.75 0 0.2 0.2];
        end
    end
	
    % annotation 'Median'
    annotation(gcf,'textbox',...
        pos1,...
        'String','Median',...
        'LineStyle','none',...
        'HorizontalAlignment','right',...
        'FitBoxToText','on');
    
    % annotation Median Values
    medianValueString = [''];
    for ci = 1 : length(cond_cell)
        pdcond  = cond_cell{ci};
        if ci > 1
            medianValueString = [medianValueString, newline];
        end
        eval(['t_pdcond = t_' pdcond ';'])
        medianValueString = [medianValueString [pdcond ': ' num2str(median(t_pdcond * 1000)) ' ms, ' num2str(length(t_pdcond)) ' trials']];
        
        clear pdcond t_pdcond
    end
    annotation(gcf,'textbox',...
        pos2,...
        'String',medianValueString,...
        'LineStyle','none',...
        'HorizontalAlignment','right',...
        'FitBoxToText','off');
    clear medianValueString
    
    
    
    % --- text for 25% and 75% --- %
    bp = findobj(gca, 'Tag', 'boxplot');
    boxes = findobj(bp, 'Tag', 'Box');
    
    for ci = 1 : length(cond_cell)
        pdcond  = cond_cell{ci};
        eval(['t_pdcond = t_' pdcond ';'])
        bi = length(cond_cell) - ci + 1;
        text(boxes(bi).XData(3)* 1.01, boxes(bi).YData(3), [num2str(quantile(t_pdcond, 0.75) * 1000) ' ms'])
        text(boxes(bi).XData(4)* 1.01, boxes(bi).YData(4), [num2str(quantile(t_pdcond, 0.25) * 1000) ' ms'])
        
        clear t_pdcond idx pdcond
    end
    
    set(gca, 'XLim', get(gca, 'XLim')+ 0.2)
    clear bp boxes
    
    
    clear task position ts
    clear t_normal t_mild t_moderate
    clear ps p1 p2 p3 xt yt 
end

%---- save figure ---%
savefile = fullfile(savefolder, savefilename);
% save figure
saveas(gcf, savefile, 'png');









