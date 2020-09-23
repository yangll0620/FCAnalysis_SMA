function m1_SKTData_timeStatiscal()


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
[correspipelinefolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% input setup
folder_input = fullfile(codecorresParentfolder, 'm0_SKTData_extract');

%  animal 
[i,j]= regexp(folder_input, 'NHP_\w*');
animal = folder_input(i + length('NHP_'):j);


%% save setup
savefolder = correspipelinefolder;
savefilename = [animal '_timeStaComp'];




%% code Starting Here

%%% extract tbl_normal, tbl_mild and tbl_moderate %%%
pdcond  = 'normal';
files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
t_event = [];
for fi = 1: length(files)
    load(fullfile(folder_input, files(fi).name), 'T_idxevent', 'fs');
    
    t_event = cat(1, t_event, T_idxevent{:, :}/ fs);
    
    clear T_idxevent
end
t_reach = t_event(:, 3) - t_event(:, 2);
t_return = t_event(:, 5) - t_event(:, 3);
tbl_normal = table(t_reach, t_return);
clear t_reach t_return pdcond files t_event


pdcond  = 'mild';
files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
t_event = [];
for fi = 1: length(files)
    load(fullfile(folder_input, files(fi).name), 'T_idxevent', 'fs');
    
    t_event = cat(1, t_event, T_idxevent{:, :}/ fs);
    
    clear T_idxevent
end
t_reach = t_event(:, 3) - t_event(:, 2);
t_return = t_event(:, 5) - t_event(:, 3);
tbl_mild = table(t_reach, t_return);
clear t_reach t_return pdcond files t_event



pdcond  = 'moderate';
files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
t_event = [];
for fi = 1: length(files)
    load(fullfile(folder_input, files(fi).name), 'T_idxevent', 'fs');
    
    t_event = cat(1, t_event, T_idxevent{:, :}/ fs);
    
    clear T_idxevent
end
t_reach = t_event(:, 3) - t_event(:, 2);
t_return = t_event(:, 5) - t_event(:, 3);
tbl_moderate = table(t_reach, t_return);
clear t_reach t_return pdcond files t_event





%%% plot normal, mild and moderate comparison %%%
figure('Position',[675 549 570* 2 413]);

% position of each subplot
x0 = 0.05; y = 0.11;
w =  0.4; h = 0.815;  gap = 0.1;


alpha = 0.05;

n_normal = height(tbl_normal);
n_mild = height(tbl_mild);
n_moderate = height(tbl_moderate);

g1 = repmat({'normal'}, n_normal,1);
g2 = repmat({'mild'}, n_mild,1);
g3 = repmat({'moderate'}, n_moderate,1);
g = [g1; g2; g3];


%----- reach time comparison ----------%
task = 'reach'; i = 1;
position = [x0 + (w + gap) * (i-1) y w h];

% extract t_normal, t_mild and t_moderate
eval(['t_normal = tbl_normal.t_' task ';']);
eval(['t_mild = tbl_mild.t_' task ';']);
eval(['t_moderate = tbl_moderate.t_' task ';']);


ts = [t_normal; t_mild; t_moderate];


% wilcoxon rank sum test
p1 = ranksum(t_normal, t_mild);
p2 = ranksum(t_mild, t_moderate);
p3 = ranksum(t_normal, t_moderate);

% actual plot 
subplot('Position', position); 
boxplot(ts, g); hold on
title([animal ':' task ' time comparison'])

% significant part
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
axis([xlim    0  max(yt) * 1.1])

if p1 < alpha
    plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
end
if p2 < alpha
    plot(xt([2 3]), [1 1]*max(yt)*0.97, '-k',  mean(xt([2 3])), max(yt) * 0.99, '*k')
end
if p3 < alpha
    plot(xt([1 3]), [1 1]*max(yt)*1.02, '-k',  mean(xt([1 2])), max(yt) * 1.04, '*k')
end


% Create textbox
annotation(gcf,'textbox',...
    [0.25 0.05 0.2 0.2],...
    'String',[['normal: ' num2str(median(t_normal * 1000)) ' ms, ' num2str(n_normal) ' trials'], newline, ...
              ['mild: ' num2str(median(t_mild * 1000)) ' ms, ' num2str(n_mild) ' trials'], newline, ...
              ['moderate: ' num2str(median(t_moderate * 1000)) ' ms, ' num2str(n_moderate) ' trials']],...
    'LineStyle','none',...
    'HorizontalAlignment','right',...
    'FitBoxToText','off');




%----- return time comparison ----------%
task = 'return';  i = 2;
position = [x0 + (w + gap) * (i-1) y w h];

% extract t_normal, t_mild and t_moderate
eval(['t_normal = tbl_normal.t_' task ';']);
eval(['t_mild = tbl_mild.t_' task ';']);
eval(['t_moderate = tbl_moderate.t_' task ';']);


ts = [t_normal; t_mild; t_moderate];


% wilcoxon rank sum test
p1 = ranksum(t_normal, t_mild);
p2 = ranksum(t_mild, t_moderate);
p3 = ranksum(t_normal, t_moderate);


% actual plot 
subplot('Position', position);
boxplot(ts, g); hold on
title([animal ':' task ' time comparison'])

% significant part
xt = get(gca, 'XTick');
yt = get(gca, 'YTick');
axis([xlim    0  max(yt) * 1.1])

if p1 < alpha
    plot(xt([1 2]), [1 1]*max(yt)*0.92, '-k',  mean(xt([1 2])), max(yt) * 0.94, '*k')
end
if p2 < alpha
    plot(xt([2 3]), [1 1]*max(yt)*0.97, '-k',  mean(xt([2 3])), max(yt) * 0.99, '*k')
end
if p3 < alpha
    plot(xt([1 3]), [1 1]*max(yt)*1.02, '-k',  mean(xt([1 2])), max(yt) * 1.04, '*k')
end

% Create textbox
% Create textbox
annotation(gcf,'textbox',...
    [0.75 0.05 0.2 0.2],...
    'String',[['normal: ' num2str(median(t_normal * 1000)) ' ms, ' num2str(n_normal) ' trials'], newline, ...
              ['mild: ' num2str(median(t_mild * 1000)) ' ms, ' num2str(n_mild) ' trials'], newline, ...
              ['moderate: ' num2str(median(t_moderate * 1000)) ' ms, ' num2str(n_moderate) ' trials']],...
    'LineStyle','none',...
    'HorizontalAlignment','right',...
    'FitBoxToText','off');




%---- save figure ---%
savefile = fullfile(savefolder, savefilename);
saveas(gcf, savefile, 'png');







