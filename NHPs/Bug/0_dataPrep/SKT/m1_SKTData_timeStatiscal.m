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
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% input setup: For Kitty normal and mild in two folders
folder_input = fullfile(codecorresParentfolder, 'm0_SKTData_extract');

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
savefilename = [animal '_timeStaComp'];




%% code Starting Here

coli_reachonset = 2;
coli_touch = 3;
coli_returnonset = 4;
coli_mouth = 5;


%%% extract tbl_normal, tbl_mild and tbl_moderate %%%
pdcond  = 'normal';
files = dir(fullfile(folder_input, ['*_' pdcond '_*.mat']));
t_event = [];
for fi = 1: length(files)
    load(fullfile(folder_input, files(fi).name), 'T_idxevent', 'fs');
    
    t_event = cat(1, t_event, T_idxevent{:, :}/ fs);
    
    clear T_idxevent
end
t_reach = t_event(:, coli_touch) - t_event(:, coli_reachonset);
t_return = t_event(:, coli_mouth) - t_event(:, coli_returnonset);
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
t_reach = t_event(:, coli_touch) - t_event(:, coli_reachonset);
t_return = t_event(:, coli_mouth) - t_event(:, coli_returnonset);
tbl_mild = table(t_reach, t_return);
clear t_reach t_return pdcond files t_event


%%% remove the trials with super longer reach and return
t_maxreach = 2;
t_maxreturn = 2;
tbl_normal = tbl_normal(tbl_normal.t_reach<=t_maxreach & tbl_normal.t_return<=t_maxreturn,:);
tbl_mild = tbl_mild(tbl_mild.t_reach<=t_maxreach & tbl_mild.t_return<=t_maxreturn,:);


%%% plot normal, mild and mild comparison %%%
figure('Position',[675 549 570* 2 413]);

% position of each subplot
x0 = 0.05; y = 0.11;
w =  0.4; h = 0.815;  gap = 0.1;


alpha = 0.05;

n_normal = height(tbl_normal);
n_mild = height(tbl_mild);

g1 = repmat({'normal'}, n_normal,1);
g2 = repmat({'mild'}, n_mild,1);
g = [g1; g2];


%----- reach time comparison ----------%
task = 'reach'; i = 1;
position = [x0 + (w + gap) * (i-1) y w h];

% extract t_normal and t_moderate
eval(['t_normal = tbl_normal.t_' task ';']);
eval(['t_mild = tbl_mild.t_' task ';']);


ts = [t_normal; t_mild];


% wilcoxon rank sum test
p1 = ranksum(t_normal, t_mild);

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


% Create textbox
annotation(gcf,'textbox',...
    [0.32 0.23 0.06 0.07],...
    'String','Median',...
    'LineStyle','none',...
    'HorizontalAlignment','right',...
    'FitBoxToText','on');
annotation(gcf,'textbox',...
    [0.25 0.05 0.2 0.2],...
    'String',[['normal: ' num2str(median(t_normal * 1000)) ' ms, ' num2str(n_normal) ' trials'], newline, ...
              ['mild: ' num2str(median(t_mild * 1000)) ' ms, ' num2str(n_mild) ' trials']],...
    'LineStyle','none',...
    'HorizontalAlignment','right',...
    'FitBoxToText','off');


% --- text for 25% and 75% --- %
bp = findobj(gca, 'Tag', 'boxplot');
boxes = findobj(bp, 'Tag', 'Box');

% mild
text(boxes(1).XData(3)* 1.01, boxes(1).YData(3), [num2str(quantile(t_mild, 0.75) * 1000) ' ms'])
text(boxes(1).XData(4)* 1.01, boxes(1).YData(4), [num2str(quantile(t_mild, 0.25) * 1000) ' ms'])
% normal
text(boxes(2).XData(3)* 1.01, boxes(2).YData(3), [num2str(quantile(t_normal, 0.75) * 1000) ' ms'])
text(boxes(2).XData(4)* 1.01, boxes(2).YData(4), [num2str(quantile(t_normal, 0.25) * 1000) ' ms'])
set(gca, 'XLim', get(gca, 'XLim')+ 0.2)
clear bp boxes



%----- return time comparison ----------%
task = 'return';  i = 2;
position = [x0 + (w + gap) * (i-1) y w h];

% extract t_normal, t_mild and t_moderate
eval(['t_normal = tbl_normal.t_' task ';']);
eval(['t_mild = tbl_mild.t_' task ';']);


ts = [t_normal; t_mild];


% wilcoxon rank sum test
p1 = ranksum(t_normal, t_mild);


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


% Create textbox
annotation(gcf,'textbox',...
    [0.82 0.23 0.06 0.07],...
    'String','Median',...
    'LineStyle','none',...
    'HorizontalAlignment','right',...
    'FitBoxToText','on');
annotation(gcf,'textbox',...
    [0.75 0.05 0.2 0.2],...
    'String',[['normal: ' num2str(median(t_normal * 1000)) ' ms, ' num2str(n_normal) ' trials'], newline, ...
              ['mild: ' num2str(median(t_mild * 1000)) ' ms, ' num2str(n_mild) ' trials']],...
    'LineStyle','none',...
    'HorizontalAlignment','right',...
    'FitBoxToText','off');

% --- text for 25% and 75% --- %
bp = findobj(gca, 'Tag', 'boxplot');
boxes = findobj(bp, 'Tag', 'Box');

% mild
text(boxes(1).XData(3)* 1.01, boxes(1).YData(3), [num2str(quantile(t_mild, 0.75) * 1000) ' ms'])
text(boxes(1).XData(4)* 1.01, boxes(1).YData(4), [num2str(quantile(t_mild, 0.25) * 1000) ' ms'])
% normal
text(boxes(2).XData(3)* 1.01, boxes(2).YData(3), [num2str(quantile(t_normal, 0.75) * 1000) ' ms'])
text(boxes(2).XData(4)* 1.01, boxes(2).YData(4), [num2str(quantile(t_normal, 0.25) * 1000) ' ms'])
set(gca, 'XLim', get(gca, 'XLim')+ 0.2)
clear bp boxes



%---- save-----%
savefile = fullfile(savefolder, savefilename);

% save figure
saveas(gcf, savefile, 'png');

% save task time 
rowNames = {'normalCOT', 'moderateSKT'};

coli = 1;
t_median = [median(tbl_normal{:,coli});  median(tbl_mild{:,coli})];
t_10 = [quantile(tbl_normal{:,coli}, 0.1);  quantile(tbl_mild{:,coli}, 0.1)];
t_25 = [quantile(tbl_normal{:,coli}, 0.25);  quantile(tbl_mild{:,coli}, 0.25)];
t_75 = [quantile(tbl_normal{:,coli}, 0.75); quantile(tbl_mild{:,coli}, 0.75)];
t_90 = [quantile(tbl_normal{:,coli}, 0.9);  quantile(tbl_mild{:,coli}, 0.9)];
tbl_staReachT = table(t_median, t_10, t_25, t_75, t_90,'rowNames', rowNames);

coli = 2;
t_median = [median(tbl_normal{:,coli});  median(tbl_mild{:,coli})];
t_10 = [quantile(tbl_normal{:,coli}, 0.1);  quantile(tbl_mild{:,coli}, 0.1)];
t_25 = [quantile(tbl_normal{:,coli}, 0.25); quantile(tbl_mild{:,coli}, 0.25)];
t_75 = [quantile(tbl_normal{:,coli}, 0.75);  quantile(tbl_mild{:,coli}, 0.75)];
t_90 = [quantile(tbl_normal{:,coli}, 0.9);  quantile(tbl_mild{:,coli}, 0.9)];
tbl_staReturnT = table(t_median, t_10, t_25, t_75, t_90, 'rowNames', rowNames);

save(fullfile(savefolder, [animal '_SKTTime.mat']), 'tbl_normal', 'tbl_mild', 'tbl_staReachT', 'tbl_staReturnT')