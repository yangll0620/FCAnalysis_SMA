function m2_SKTData_timeStatiscal()


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
folder_input = fullfile(codecorresParentfolder, 'm1_SKTData_avgArea');

%  animal 
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


tasks = {'reach', 'return'};

%% save setup
savefolder = codecorresfolder;
savefilename = [animal '_timeStaComp'];




%% code Starting Here
coli_targetonset = 1;
coli_reachonset = 2;
coli_touch = 3;
coli_returnonset = 4;
coli_mouth = 5;


%%% extract tbl_normal, tbl_mild and tbl_moderate %%%
cond_cell = cond_cell_extract(animal);
for ci = 1 : length(cond_cell)
    
    pdcond  = cond_cell{ci};
    files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
    t_event = [];
    for fi = 1: length(files)
        load(fullfile(folder_input, files(fi).name), 'T_idxevent_lfp', 'fs_lfp');
        
        t_event = cat(1, t_event, T_idxevent_lfp{:, :}/ fs_lfp);
        
        clear T_idxevent_lfp fs_lfp
    end
    t_reaction = t_event(:, coli_reachonset) - t_event(:, coli_targetonset);
    t_reach = t_event(:, coli_touch) - t_event(:, coli_reachonset);
    t_return = t_event(:, coli_mouth) - t_event(:, coli_returnonset);
    eval(['tbl_' pdcond ' = table(t_reaction, t_reach, t_return);'])
    
    clear t_reaction t_reach t_return pdcond files t_event
end



%%% plot pdcond comparison %%%
figure('Position',[675 549 570* 2 413]);

% position of each subplot
x0 = 0.05; y = 0.11;
w =  0.4; h = 0.815;  gap = 0.1;


alpha = 0.05;
gs =[];
for ci = 1 : length(cond_cell)
    pdcond  = cond_cell{ci};
    
    eval(['n = height(tbl_' pdcond ');'])
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
        eval(['t_' pdcond ' = tbl_' pdcond '.t_' task ';']);
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









