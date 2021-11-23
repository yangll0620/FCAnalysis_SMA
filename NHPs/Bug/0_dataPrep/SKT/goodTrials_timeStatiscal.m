clear
% add util path
codefolder = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\code';
addpath(genpath(fullfile(codefolder,'util')));


mafolder_input = 'C:\Users\nmrc3\Desktop\Bug Data\';
animal = 'Bug';

folder_input2 = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\Bug\Recording\Processed\DataDatabase';

%% code Start Here


coli_reachonset = 2;
coli_touch = 3;
coli_returnonset = 4;
coli_mouth = 5;



%%% extract tbl_normal, tbl_mild and tbl_moderate %%%
folders1 = dir(fullfile(mafolder_input, [animal '_' '*']));
folders2 = dir(fullfile(folder_input2, [animal '_' '*']));
folders = [folders1; folders2];
t_event_normal = [];
t_event_mild = [];
for fi = 1: length(folders)
    folder = fullfile(folders(fi).folder, folders(fi).name);
    
    mafolderstruct = dir(fullfile(folder, 'Block-*'));
    if isempty(mafolderstruct)
        continue;
    end
    
    mafolder = fullfile(mafolderstruct.folder, mafolderstruct.name);
    mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));
    
    % ma file does not exist
    if isempty(mafilestruct)
        continue;
    end
    
    % get the pd conditioon for the date of experiment
    tmps = regexp(mafolder, [animal '_[0-9]*'],'match');
    dateofexp = datenum(tmps{1}(length([animal '_']) + 1:end), 'mmddyy');
    pdcond = parsePDCondition(dateofexp, animal);
    tmps = regexp(mafolder, 'Block-[0-9]*','match');
    bktdt = str2num(tmps{1}(length('Block-') + 1:end));
    
    
    % load SingleTargetKluverMAData
    load(fullfile(mafolder, mafilestruct.name), 'SingleTargetKluverMAData');
    
    
    % the tag of good reach trials
    tag_goodreach = SingleTargetKluverMAData.goodix_reach;
    % the tag of good return trials
    tag_goodreturn = SingleTargetKluverMAData.goodix_return;
    % extract indices of good trials (both have good reach and return)
    idx_goodtrials = find(tag_goodreach .* tag_goodreturn == 1);
  
    disp([pdcond ', ' datestr(dateofexp, 'yyyymmdd') '-Block' num2str(bktdt) ', number of good trials = ' num2str(length(idx_goodtrials))])
      
    
    % ma sample rate
    fs_ma = SingleTargetKluverMAData.SR;
    
    % time indices for target onset, reach onset, touch screen, return and mouth
    TargetTime = SingleTargetKluverMAData.TargetTime;
    ReachTimeix = SingleTargetKluverMAData.ReachTimeix;
    TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
    ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix
    MouthTimeix = SingleTargetKluverMAData.MouthTimeix;
    
    timeixtbl_ma = [table(TargetTime) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];
    timeixtbl_ma = timeixtbl_ma(idx_goodtrials,:);
    
    t_event_today = timeixtbl_ma{:, :}/ fs_ma;
    t_reach_today = t_event_today(:, coli_touch) - t_event_today(:, coli_reachonset);
    t_return_today = t_event_today(:, coli_mouth) - t_event_today(:, coli_returnonset);
    disp(['reach time range = [' num2str(min(t_reach_today)) ' ' num2str(max(t_reach_today)) ']s' ', ' ...
          'return time range = [' num2str(min(t_return_today)) ' ' num2str(max(t_return_today)) ']s'])
    disp(['avg reach time= ' num2str(mean(t_reach_today)) 's' ', ' ...
          'avg return time = ' num2str(mean(t_return_today)) 's'])
    disp('\n')
    
    
    if strcmp(pdcond, 'normal')
        t_event_normal = cat(1, t_event_normal, t_event_today);
    end
    if strcmp(pdcond, 'mild')
        t_event_mild = cat(1, t_event_mild, t_event_today);
    end
    
    clear folder mafolderstruct mafolder mafilestruct
    clear tmps dateofexp pdcond bktdt
    clear SingleTargetKluverMAData  tag_goodreach  tag_goodreturn idx_goodtrials
    clear fs_ma TargetTime ReachTimeix TouchTimeix ReturnTimeix MouthTimeix timeixtbl_ma
    clear t_event_today t_reach_today t_return_today
end

if ~isempty(t_event_normal)
    t_reach = t_event_normal(:, coli_touch) - t_event_normal(:, coli_reachonset);
    t_return = t_event_normal(:, coli_mouth) - t_event_normal(:, coli_returnonset);
    tbl_normal = table(t_reach, t_return);
    clear t_reach t_return
end
if ~isempty(t_event_mild)
    t_reach = t_event_mild(:, coli_touch) - t_event_mild(:, coli_reachonset);
    t_return = t_event_mild(:, coli_mouth) - t_event_mild(:, coli_returnonset);
    tbl_mild = table(t_reach, t_return);
    clear t_reach t_return
end



%%% plot normal, mild and moderate comparison %%%
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


