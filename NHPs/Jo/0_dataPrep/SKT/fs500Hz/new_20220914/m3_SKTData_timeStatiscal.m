function m3_SKTData_timeStatiscal()

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


[~, codefilename]= fileparts(codefilepath);

% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% global variables

% animal
animal = animal_extract(codecorresfolder);

tasks = {'reaction','reach', 'return'};

%% input

folder_input = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

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
    files = dir(fullfile(folder_input, ['*' pdcond '_*.mat']));
    t_event = [];
    for fi = 1: length(files)
        file = fullfile(folder_input, files(fi).name);

        load(file, 'tbl_goodTrialsMarks', 'T_idxevent_lfp', 'fs_lfp');
        goodTrials_allGs = tbl_goodTrialsMarks{:, :};
        goodTrials = ones(size(goodTrials_allGs, 1), 1);
        for gi = 1 : size(goodTrials_allGs, 2)
            goodTrials = goodTrials & goodTrials_allGs(:, gi);
        end
        
        % use only the good Trials
        T_idxevent_lfp = T_idxevent_lfp(goodTrials, :);
        
        t_event = cat(1, t_event, T_idxevent_lfp{:, :}/ fs_lfp);
        
        clear T_idxevent_lfp fs_lfp listOfVariables goodTrials
    end
    t_reaction = t_event(:, coli_reachonset) - t_event(:, coli_targetonset);
    t_reach = t_event(:, coli_touch) - t_event(:, coli_reachonset);
    t_return = t_event(:, coli_mouth) - t_event(:, coli_returnonset);
    eval(['tbl_' pdcond ' = table(t_reaction, t_reach, t_return);'])
    
    clear t_reaction t_reach t_return pdcond files t_event
end



%%% plot pdcond comparison %%%

sigalpha = 0.05;
gs =[];
for ci = 1 : length(cond_cell)
    pdcond  = cond_cell{ci};
    
    eval(['n = height(tbl_' pdcond ');'])
    g = repmat({pdcond}, n,1);
    gs = [gs; g];
    clear pdcond
end


%----- time comparison ----------%
figure('Position',[675 549 570*3 400]);

% position of each subplot
ncols = length(tasks);
nrows = 1;
marg_left = 0.05; marg_right = 0.05;
marg_top = 0.1; marg_bottom = 0.05; 
deltax = 0.1; deltay = 0.05;

x0 = marg_left; y0 = marg_bottom;
w_subplot =  floor(((1-marg_left-marg_right-(ncols-1)*deltax)/ncols)*100)/100;
h_subplot =  floor(((1-marg_top-marg_bottom-(nrows-1)*deltay)/nrows)*100)/100;
for ti = 1 : length(tasks)
    task = tasks{ti};
    position = [x0 + (w_subplot + deltax) * (ti-1) y0 w_subplot h_subplot];
    
    
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
    ylabel('time/s')
    
    
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
    
    % --- set ylim --- %
    bp = findobj(gca, 'Tag', 'boxplot');
    uppers = findobj(bp, 'Tag', 'Upper Whisker');
    maxYs = [];
    for upi = 1 : length(uppers)
        maxYs = [maxYs; uppers(upi).YData(2)];
    end
    ylimit = ylim;
    if ylimit(2) > 10
        ylim([0 max(maxYs)*1.1])
    end
    clear bp uppers maxYs ylimit
    
    
    % --- significant part --- %
    xt = get(gca, 'XTick');
    yt = get(gca, 'YTick');
    %axis([xlim    0  max(yt) * 1.1])
    
    if exist('p1', 'var') && exist('p2', 'var') && exist('p3', 'var')
        if p1 < sigalpha
            plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
        end
        if p2 < sigalpha
            plot(xt([2 3]), [1 1]*max(yt)*0.97, '-k',  mean(xt([2 3])), max(yt) * 0.99, '*k')
        end
        if p3 < sigalpha
            plot(xt([1 3]), [1 1]*max(yt)*1.02, '-k',  mean(xt([1 2])), max(yt) * 1.04, '*k')
        end
    else
        if exist('p1', 'var')
            if p1 < sigalpha
                plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
            end
        end
        
        if exist('p2', 'var')
            if p2 < sigalpha
                plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
            end
        end
        
        if exist('p3', 'var')
            if p3 < sigalpha
                plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
            end
        end
    end
     
    
    
    
    % Create textbox
    medianValueString = ['Median'];
    for ci = 1 : length(cond_cell)
        pdcond  = cond_cell{ci};
        eval(['t_pdcond = t_' pdcond ';'])
        
        medianValueString = [medianValueString, newline];
        medianValueString = [medianValueString [pdcond ': ' num2str(median(t_pdcond * 1000)) ' ms, ' num2str(length(t_pdcond)) ' trials']];
        
        clear pdcond t_pdcond
    end
    xlimit = xlim; ylimit = ylim;
    x = (xlimit(1) + xlimit(2))/2;
    y = (ylimit(1) + ylimit(2))/3*2;
    text(x, y, medianValueString, 'LineStyle', 'none');
    clear medianValueString x y xlimit ylimit
    

    
    clear task position ts
    clear t_normal t_mild t_moderate
    clear ps p1 p2 p3 xt yt 
end

%---- save figure ---%
savefile = fullfile(savefolder, savefilename);
% save figure
saveas(gcf, savefile, 'png');









