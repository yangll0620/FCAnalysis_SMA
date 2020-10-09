clear
close all
[reachTimes_normal, returnTimes_normal, reachTimes_mild, returnTimes_mild]= reachReturnAnalysis();

%% output
fileID = fopen('reachReturnStatics.txt','w');            %here I create the text file

%% welch ttest for mean values if normal
if kstest(zscore(reachTimes_normal)) || kstest(zscore(reachTimes_mild)) 
    reach_welch = 'reach Time of normal and mild are not both normal';
else
    
    % welch t-test
    [~, p_reach] = ttest2(returnTimes_mild, returnTimes_normal);
    
    reach_welch = ['mean reach:'  'normal = ' num2str(mean(reachTimes_normal)) ', mild = ' num2str(mean(reachTimes_mild)), ...
        ', welch ttest p = ' num2str(p_reach)];
end

if kstest(zscore(returnTimes_normal)) || kstest(zscore(returnTimes_mild))
    return_welch = 'return Time of normal and mild are not both normal';
else
    
    [~, p_return] = ttest2(returnTimes_mild, returnTimes_normal);
    
    return_welch = ['mean return:'  'normal = ' num2str(mean(returnTimes_normal)) ', mild = ' num2str(mean(returnTimes_mild)), ...
        ', welch ttest p = ' num2str(p_return)];
end

%% Wilcoxon rank sum test
p_reach = ranksum(reachTimes_normal,reachTimes_mild);
reach_wilcoxon = ['median reach:'  'normal = ' num2str(median(reachTimes_normal)) ', mild = ' num2str(median(reachTimes_mild)), ...
        ', welch ttest p = ' num2str(p_reach)];
    
p_return = ranksum(returnTimes_normal,returnTimes_mild); 
return_wilcoxon = ['median return:'  'normal = ' num2str(median(returnTimes_normal)) ', mild = ' num2str(median(returnTimes_mild)), ...
        ', welch ttest p = ' num2str(p_return)];

    
%% display and write to text
disp(reach_welch)
disp(return_welch)
disp(reach_wilcoxon)
disp(return_wilcoxon)

fprintf(fileID,  reach_welch);
fprintf(fileID,  '\n');
fprintf(fileID,  return_welch);
fprintf(fileID,  '\n\n');
fprintf(fileID,  reach_wilcoxon);
fprintf(fileID,  '\n');
fprintf(fileID,  return_wilcoxon);

fclose(fileID);


%% plot
reachTime = cat(1, reachTimes_normal, reachTimes_mild);
returnTime = cat(1, returnTimes_normal, returnTimes_mild);
ntrials_normal = length(reachTimes_normal);
ntrials_mild = length(reachTimes_mild);
g_normal = repmat({'normal'}, ntrials_normal,1);
g_mild = repmat({'mild'}, ntrials_mild,1);
g = cat(1, g_normal, g_mild);


% reach
h1 = figure;
boxplot(reachTime, g)
title('Reach Time Comparison')
string = {"$t_{normal}=$" + num2str(median(reachTimes_normal)), ...
    "$t_{mild}=$" + num2str(median(reachTimes_mild)),...
    "",...
    "Wilcoxon rank sum test",...
    "$p=$" + num2str(p_reach)};
annotation(h1,'textbox',[0.3875 0.77 0.05 0.05],...
    'String', string, 'interpreter','latex', 'FitBoxToText','on', 'LineStyle', 'none');

saveas(h1, 'reach_Time_Compare', 'png');


% return
h2 = figure;
boxplot(returnTime, g)
title('Return Time Comparison')
string = {"$t_{normal}=$" + num2str(median(returnTimes_normal)), ...
    "$t_{mild}=$" + num2str(median(returnTimes_mild)),...
    "", ...
    "Wilcoxon rank sum test",...
    "$p=$" + num2str(p_return)};
annotation(h2,'textbox',[0.3875 0.77 0.05 0.05],...
    'String', string, 'interpreter','latex', 'FitBoxToText','on', 'LineStyle', 'none');

% save
saveas(h2, 'return_Time_Compare', 'png');