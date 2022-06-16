function m4_fs500Hz_freezeSKTData_Time_Histogram()
% Objective:
%       extract freezing episode
%       data structure
%           episodes{}


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

% animal
animal = animal_extract(codecorresfolder);



%%  input setup

inputfolder = fullfile(codecorresParentfolder, 'm3_fs500Hz_freezeSKTData_EpisodeExtract');
pdcond = 'moderate';

img_format = 'jpg';

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);


%% Code Start Here
optFreezeTypes = optFreezeTypes_extract('codesavefolder', savecodefolder);

% extract tdur_freezes: nepisodes * length(optFreezeTypes)
tdur_freezes = zeros(1, length(optFreezeTypes));
datebkstrs = {};
files = dir(fullfile(inputfolder, ['*' pdcond '*.mat']));
for fi = 1 : length(files)
    filename = files(fi).name;
    tmp = regexp(filename, '\d{8}_bktdt\d{1}', 'match');
    datebkstr = tmp{1};
    datebkstr = strrep(datebkstr, '_bktdt', '-bk');
    clear tmp
    
    load(fullfile(inputfolder, filename), 'freezStruct', 'selectedTrials');
    freezEpisodes = freezStruct.freezEpisodes;
    for frzi = 1 : length(freezEpisodes)
        tri = freezEpisodes{frzi}.triali;
        if ~selectedTrials(tri)
            clear tri
            continue;
        end
        
        freezeType = freezEpisodes{frzi}.freezeType;
        t_freeze = freezEpisodes{frzi}.freezeTPhaseS(:, 2) - freezEpisodes{frzi}.freezeTPhaseS(:, 1);
        mask_freeze = strcmp(freezeType, optFreezeTypes);
        tmp = zeros(1, length(optFreezeTypes));
        tmp(mask_freeze) = t_freeze;
        tdur_freezes = [tdur_freezes; tmp];
        datebkstrs{end+1,1} = datebkstr;
        clear freezeType t_freeze mask_freeze tmp
    end
    
    clear filename freezStruct datebkstr
end 
tdur_freezes(1,:)= [];



%%% plot tdur_freezes

tbl_DateFreezeT = [table(datebkstrs), array2table(tdur_freezes, 'VariableNames', optFreezeTypes)];
tdursum_freeze_datestrs = [];
uniq_datebkstrs = unique(tbl_DateFreezeT.datebkstrs);
for ui = 1 : length(uniq_datebkstrs)
    datebkstr = uniq_datebkstrs{ui};
    masks = cellfun(@(x) strcmp(x, datebkstr), tbl_DateFreezeT.datebkstrs);
    
    tdursum_freeze_datestr = sum(tbl_DateFreezeT{masks, 2:end}, 1);
    tdursum_freeze_datestrs = [tdursum_freeze_datestrs; tdursum_freeze_datestr];
    
    clear datebkstr tdursum_freeze_datestr
end

% add up React-Reach and Reach
freezeTypes_plot = {'initFreeze', 'reachFreeze', 'maniFreeze'};
tdur_sumDates_freeze= [tdursum_freeze_datestrs(:, 1) tdursum_freeze_datestrs(:, 2) + tdursum_freeze_datestrs(:, 3) tdursum_freeze_datestrs(:, 4)];
tdur_sum_freeze = [tdur_sumDates_freeze; sum(tdur_sumDates_freeze,1)];
tbl_sumDateFreezeT = array2table(tdur_sum_freeze, 'VariableNames', freezeTypes_plot, 'RowNames', [uniq_datebkstrs; {'Sum'}]);
clear tdur_sumDates_freeze tdur_sum_freeze

%%% plot bar chart for each day and sum
datPlot = tbl_sumDateFreezeT{:, :};
bar(datPlot, 'stacked')
xticklabels(tbl_sumDateFreezeT.Properties.RowNames)
xtickangle(45)
ylabel('freeze duration/s')
legend(freezeTypes_plot, 'Location','northwest')
title([animal ' Freezing Time Accumulation'])

for di = 1 : size(datPlot, 1)
    Ys = round(datPlot(di, :));
    str = {};
    for ci = length(Ys): -1 : 1
        str = [str, num2str(Ys(ci))];
    end
    if di < size(datPlot, 1)
        text(di-0.25, ceil(sum(Ys)) + 100, str);
    else
        for ci = 1: length(Ys)
            text(di+0.55, sum(Ys(1:ci-1))+(Ys(ci)/2), num2str(Ys(ci)));
        end
    end
end

saveas(gcf, fullfile(savefolder, [animal '_FreezeTimeAcc.' img_format]));



ts_freezInit = tbl_DateFreezeT{:, 2}(tbl_DateFreezeT{:, 2}~= 0); 
tmp1 = tbl_DateFreezeT{:, 3}(tbl_DateFreezeT{:, 3}~= 0);
tmp2 = tbl_DateFreezeT{:, 4}(tbl_DateFreezeT{:, 4}~= 0);
ts_freezReach = [tmp1;tmp2]; 
ts_freezMani = tbl_DateFreezeT{:, 5}(tbl_DateFreezeT{:,5}~= 0);
n_init = length(ts_freezInit);
n_reach = length(ts_freezReach);
n_mani = length(ts_freezMani);
clear tmp1 tmp2

% histogram plot freeze time
figure
subplot(3,1,1);
histogram(ts_freezInit, 30); 
title([animal ' init freeze histogram, init# =' num2str(n_init)])
subplot(3,1,2);
histogram(ts_freezReach, 30)
title([animal ' reach freeze histogram, init# =' num2str(n_reach)])
subplot(3,1,3);
histogram(ts_freezMani, 30)
title([animal ' manipulate freeze histogram, init# =' num2str(n_mani)])
saveas(gcf, fullfile(savefolder, [animal '_FreezeTimeHist.' img_format]));

% box plot freeze time
x = [ts_freezInit; ts_freezReach; ts_freezMani];
freeTypes = cell(size(x));
freeTypes(1 : n_init) = {'freezInit'};
freeTypes(n_init+1 : n_init+n_reach) = {'freezReach'};
freeTypes(n_init+n_reach+1 : n_init+n_reach+n_mani) = {'freezMani'};
figure;
boxplot(x, freeTypes)
title([animal ' freeze time, init# =' num2str(n_init) ', reach# =' num2str(n_reach) ', mani# =' num2str(n_mani)])
ylabel('freeze duration time/s')
saveas(gcf, fullfile(savefolder, [animal '_FreezeTimePlot.' img_format]));

close all
