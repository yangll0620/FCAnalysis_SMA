function fig_timeStatiscal()


codefilepath = mfilename('fullpath');


% find the codefolder
tmp = regexp(codefilepath, '.*\code', 'match');
if length(tmp) ~= 1
    disp('can not find code path correctly.')
    return;
end
codefolder = tmp{1};
clear tmp

% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));

%% Input & save
[~, ~, pipelinefolder, outputfolder] = exp_subfolders();
input_folder_J = fullfile(pipelinefolder, 'NHPs', 'Jo', '0_dataPrep', 'SKT', 'fs500Hz', 'm2_SKTData_SelectTrials');


savefolder = fullfile(outputfolder, 'results', 'figures');
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


[~, funcname, ~]= fileparts(codefilepath);
savefilename = funcname;

%% parameters
coli_reachonset = 2;
coli_touch = 3;

%% Code Start here

animal = 'Jo';

% extract tbl_normal, tbl_mild and tbl_moderate 
cond_cell_J = cond_cell_extract('Jo');
nconds = length(cond_cell_J);
for ci = 1 : nconds
    
    pdcond  = cond_cell_J{ci};
    files = dir(fullfile(input_folder_J, ['*_' pdcond '_*.mat']));
    t_event = [];
    for fi = 1: length(files)
        load(fullfile(input_folder_J, files(fi).name), 'T_idxevent_lfp', 'fs_lfp', 'goodTrials'); 
        t_event = cat(1, t_event, T_idxevent_lfp{goodTrials, :}/ fs_lfp);
        clear T_idxevent_lfp fs_lfp
    end
    t_reaction.(pdcond) = t_event(:, coli_touch) - t_event(:, coli_reachonset);
    
    clear pdcond files t_event
    clear fi 
end
clear ci



%%% plot pdcond comparison %%%
figure('Position',[675 549 570 413]);

% grouping variables and ts for boxplot
ts = [];
gs =[];
for ci = 1 : length(cond_cell_J)
    pdcond  = cond_cell_J{ci};
    
    n = length(t_reaction.(pdcond));
    pdconds = repmat({[pdcond]}, n,1);
    gs = [gs; pdconds];
    
    ts = [ts; t_reaction.(pdcond)];
    
    clear pdcond n pdconds
end


% actual box plot
boxplot(ts, gs); hold on
title(['Animal ' animal(1)], 'FontSize', 11, 'FontWeight', 'bold')
ylabel('Reach Time/s', 'FontSize', 10, 'FontWeight', 'bold')

ax = gca;
xtlabels = ax.XTickLabel;
ax.XTickLabel = '';
for xti = 1 : length(xtlabels)
    pdcond = xtlabels{xti};
    text2 = ['n=' num2str(length(t_reaction.(pdcond)))];
    newlabel = {pdcond; text2};
    text(xti, ax.YLim(1), sprintf('%s\n%s\n%s', newlabel{:,:}), 'horizontalalignment', 'center', 'verticalalignment', 'top', 'FontSize', 10, 'FontWeight', 'bold'); 
end


% significant part: wilcoxon rank sum test and plot
alpha = 0.5;
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
hs = [0.92 0.97 1.02];
for ci = 1 : nconds-1
    pdcondi  = cond_cell_J{ci};
    ti = t_reaction.(pdcondi);
    for cj = ci+1 : nconds
        pdcondj  = cond_cell_J{cj};
        tj = t_reaction.(pdcondj);
        p = ranksum(ti, tj);
        
        if p < alpha % plot sig * 
            if ci == 1 && cj == 2
                h = hs(1);
            elseif ci == 2 && cj == 3
                h = hs(2);
            elseif ci == 1 && cj == 3
                h = hs(3);
            end
            plot(xt([ci cj]), [1 1]*max(yt)*h, '-k',  mean(xt([ci cj])), max(yt) * (h+0.02), '*k')
            clear h
        end
        
        clear pdcondj tj 
    end
    clear cj
    clear pdcondi ti
end
clear ci


%%% save
print(fullfile(savefolder, savefilename), '-dpng', '-r1000')
close all

