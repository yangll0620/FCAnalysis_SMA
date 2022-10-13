function fig3_uNHP_SKT_PSD(varargin)
%
%%

% parse params
p = inputParser;
addParameter(p, 'pos_ifig', [50 50 400 150], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));
addParameter(p, 'newsavefolder', false, @(x) assert(islogical(x) && isscalar(x)));


parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;
newsavefolder = p.Results.newsavefolder;

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
[~, funcname, ~]= fileparts(codefilepath);

% variables for plotting
plotF_AOI = [8 40];




%% save setup
savefolder = fullfile(outputfolder, 'results', 'figures', 'current', funcname);
savecodefolder = fullfile(savefolder, 'code');
if(exist(savefolder, 'dir') && newsavefolder)
    rmdir(savefolder, 's')
end
if exist(savecodefolder, "dir")
    rmdir(savecodefolder, 's');
end
if ~exist(savefolder, "dir")
    mkdir(savefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);



%% Code Start here
animals = {'Jo'; 'Kitty'};

for ai = 1 : length(animals)
    animal = animals{ai};
    inputfolder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', 'm3_SKTData_PSD');
    load(fullfile(inputfolder, [animal '_psd.mat']), 'pxxs', 'f_selected', 'brainareas');
    EventPhases = fieldnames(pxxs);
    
    for ei = 1: length(EventPhases)
        eventname = EventPhases{ei};
        for bi = 1:length(brainareas)
            brainarea = brainareas{bi};
            
            % load normal and PD data
            psd_normal = pxxs.(eventname).normal.(brainarea);
            psd_PD = pxxs.(eventname).PD.(brainarea);
            
            %
            psd_normal = exp(psd_normal/10);
            psd_PD = exp(psd_PD/10);
            
            %  psd_normal, psd_mild: nfs * nsegs
            plotPSD_comp_1chn(pos_ifig, psd_normal, psd_PD, f_selected, plotF_AOI, savefolder, brainarea, animal, eventname)
            
            clear brainarea psd_normal psd_PD
        end
        clear eventname
    end
    
    clear animal inputfolder
    clear pxxs f_selected brainareas EventPhases
    close all
end



function plotPSD_comp_1chn(pos_ifig, psd_normal, psd_PD, F_all, plotF_AOI, savefolder, brainarea, animal, eventname)
%%  plot the psd comparison of normal and mild

% find the idx for F_AOI
idx_AOI = find(F_all >= plotF_AOI(1) & F_all <= plotF_AOI(2));

% colors setup
color_normal_range = [220, 220, 220] / 255;
color_normal_mean = [128, 128, 128] / 255;
color_PD_range = [255, 192, 203] / 255;
color_PD_mean = [255, 0, 0] / 255;

% plot setup
linewidth = 1.5;

% F_AOI
F_AOI = F_all(idx_AOI);
% reshape F_AOI
[m, n] = size(F_AOI); F_AOI = reshape(F_AOI, 1, m * n); clear m n

% extract the  psd of AOI
psd_normal_FAOI = psd_normal(idx_AOI, :);
psd_moderate_FAOI = psd_PD(idx_AOI, :);

psd_normal_high = max(psd_normal_FAOI, [], 2);
psd_normal_low = min(psd_normal_FAOI, [], 2);
psd_normal_mean = mean(psd_normal_FAOI, 2);

psd_PD_high = max(psd_moderate_FAOI, [], 2);
psd_PD_low = min(psd_moderate_FAOI, [], 2);
psd_PD_mean = mean(psd_moderate_FAOI, 2);

% reshape, psd_*_high/low/mean into 1 * nfs

[m, n] = size(psd_normal_high); psd_normal_high = reshape(psd_normal_high, 1, m * n); clear m n
[m, n] = size(psd_normal_low); psd_normal_low = reshape(psd_normal_low, 1, m * n); clear m n
[m, n] = size(psd_normal_mean); psd_normal_mean = reshape(psd_normal_mean, 1, m * n); clear m n


[m, n] = size(psd_PD_high); psd_PD_high = reshape(psd_PD_high, 1, m * n); clear m n
[m, n] = size(psd_PD_low); psd_PD_low = reshape(psd_PD_low, 1, m * n); clear m n
[m, n] = size(psd_PD_mean); psd_PD_mean = reshape(psd_PD_mean, 1, m * n); clear m n

% plot range
ifig = figure('Position', pos_ifig);
set(ifig, 'PaperUnits', 'points');
fill([F_AOI flip(F_AOI)], [psd_normal_high flip(psd_normal_low)], color_normal_range, 'LineStyle', 'none')
hold all
fill([F_AOI flip(F_AOI)], [psd_PD_high flip(psd_PD_low)], color_PD_range, 'LineStyle', 'none')

% plot mean
h1 = plot(F_AOI, psd_normal_mean, 'Color', color_normal_mean, 'LineWidth', linewidth);
h3 = plot(F_AOI, psd_PD_mean, 'Color', color_PD_mean, 'LineWidth', linewidth);


xlim([min(F_AOI) max(F_AOI)])
ylim([0 0.3])
set(gca,'XTick',[10 15 20 25 30 35 40 45 50],'YTick',[0 0.1 0.2]);


xlabel('Frequency (Hz)', 'FontWeight','bold')
ylabel('Power', 'FontWeight','bold')

% legend
legend([h1, h3], {'Normal',  'PD'})


% save figure
savename = fullfile(savefolder, [animal '_SKTPsd_' brainarea '_' eventname]);
print(gcf, savename, '-painters', '-depsc')
print(gcf, savename, '-dpng', '-r1000')



clear psd_allsegs_normal  psd_allsegs_PD
clear psd_normal_FAOI  psd_PD_FAOI
clear psd_normal_high psd_normal_low psd_normal_mean psd_PD_high psd_PD_low psd_PD_mean
clear h1 h2 h3 maxPSD F_maxPSD idx_max
clear savename

