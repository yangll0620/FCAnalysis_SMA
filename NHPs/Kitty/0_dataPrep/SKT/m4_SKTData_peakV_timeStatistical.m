function m4_SKTData_peakV_timeStatistical()
%  extract lfp data respect to reachonset
% 
%  return:
%        lfptrials: nchns * ntemp * ntrials


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
animal = animal_extract(codecorresfolder);

%%  input setup

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_SKTData_SelectTrials');


saveFig_format = 'tif';

%% save setup
savefolder = codecorresfolder;
savefilename = [animal 'peakV_timeStaComp'];

%% code start here
coli_targetonset = uint32(SKTEvent.TargetOnset);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);
coli_returnonset = uint32(SKTEvent.ReturnOnset);
coli_mouth = uint32(SKTEvent.Mouth);


cond_cell = cond_cell_extract(animal);
for ci = 1: length(cond_cell)
    pdcond = cond_cell{ci};
    files = dir(fullfile(inputfolder, ['*' pdcond '*.mat']));
    
    t_reachonset2peakVs = [];
    t_peakV2reachs = [];
    for fi = 1 : length(files)
        file = fullfile(files(fi).folder, files(fi).name);
        load(file, 'fs_lfp', 'lfpdata', 'goodTrials', 'T_chnsarea', 'T_idxevent_lfp',...
            'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma');

        [nchns, ~, ntrials] = size(lfpdata);
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
        clear('fs_lfp', 'lfpdata', 'goodTrials', 'T_chnsarea', 'T_idxevent_lfp',...
              'smoothWspeed_trial', 'T_idxevent_ma', 'fs_ma');
        clear nchns ntrials tri
    end
    
    
    eval(['t_reachonset2peakVs_' pdcond ' = t_reachonset2peakVs;'])
    eval(['t_peakV2reachs_' pdcond ' = t_peakV2reachs;'])
    
    clear pdcond files t_reachonset2peakVs t_peakV2reachs
end


%%%  plot t_reachonset2peakV
ts = [];
gs = [];
for ci = 1: length(cond_cell)
    pdcond = cond_cell{ci};
    
    eval(['t_tmp = t_reachonset2peakVs_' pdcond ';'])
    ts = [ts; t_tmp];
        
    g_tmp = repmat({pdcond}, length(t_tmp),1);
    gs = [gs; g_tmp];
    
    eval(['t_' pdcond ' = t_tmp;'])
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

% 
alpha = 0.05; 
figure;
boxplot(ts, gs);
hold on;

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


title('t-duration reachonset to peakV')
saveas(gcf, fullfile(savefolder, [savefilename '_reachonset2peakV']), saveFig_format);



% plot t_peakV2reach
ts = [];
gs = [];
for ci = 1: length(cond_cell)
    pdcond = cond_cell{ci};
    
    eval(['t_tmp = t_peakV2reachs_' pdcond ';'])
    ts = [ts; t_tmp];
    
    g_tmp = repmat({pdcond}, length(t_tmp),1);
    gs = [gs; g_tmp];
end

figure
boxplot(ts, gs);
title('t-duration peakV to reach')
saveas(gcf, fullfile(savefolder, [savefilename '_peakV2reach']), saveFig_format);









