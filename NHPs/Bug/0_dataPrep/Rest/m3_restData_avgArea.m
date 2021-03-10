function m3_restData_avgArea()
%
%   averaged the lfp in one area (except the DBS)
%       only chns of Layer 5  are used for  PMC and M1
%
%
%   Input:
%
%       \m2_restData_selectSeg_Power
%
%       'data_segments', 'fs', 'T_chnarea'
%
%   Output variables:
%
%
%        lfpsegs: lfp segments in all areas
%
%        fs: samping rate (default 500Hz)
%
%        chnAreas: chnAreas cell for used in python


%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code') - 1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder, 'util')));


[datafolder, ~, ~, ~] = exp_subfolders();
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% global variables
% animal
[i, j] = regexp(codecorresfolder, 'NHPs/[A-Za-z]*');
animal = codecorresfolder(i + length('NHPs/'):j);


%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_Power');


% depth file folder
folder_depthFile = fullfile(datafolder, animal);
file_chnDepth_normal = fullfile(folder_depthFile, [animal '_channelDepth_normal.csv']);
file_chnDepth_mild = fullfile(folder_depthFile, [animal '_channelDepth_mild.csv']);
file_chnDepth_moderate = fullfile(folder_depthFile, [animal '_channelDepth_moderate.csv']);

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = 'avgArea';

strformat_date_save = 'yyyymmdd';

%% starting: narrow filter the lfp data of all the files
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);

for filei = 1:nfiles
    filename = files(filei).name;
    
    
    load(fullfile(files(filei).folder, filename), 'T_chnsarea', 'data_segments', 'fs');
    
    if contains(filename, 'normal')
        pdcondition = 'normal';
    end
    if contains(filename, 'mild')
        pdcondition = 'mild';
    end
    if contains(filename, 'moderate')
        pdcondition = 'moderate';
    end
    
    %%%  add depth variable to T_chnsarea %%%
    T_GMchnsarea = T_chnsarea(strcmp(T_chnsarea.electype, 'Gray Matter'), :);
    T_DBSchnsarea = T_chnsarea(strcmp(T_chnsarea.electype, 'DBS'), :);
    
    % find dateofexp
    idx = strfind(filename, '_tdt');
    dateofexp = datenum(filename(idx-8:idx-1), 'yyyymmdd');
    
    switch pdcondition
        case 'normal'
            file_chnDepth = file_chnDepth_normal;
        case 'mild'
            file_chnDepth = file_chnDepth_mild;
        case 'moderate'
            file_chnDepth = file_chnDepth_moderate;
    end
    
    T_GMchnsarea_Depth = add_daily_GMDepth(T_GMchnsarea, file_chnDepth, dateofexp);
    
    T_chnsarea = combine_GMDBSChns(T_GMchnsarea_Depth, T_DBSchnsarea);
    clear idx T_GMchnsarea T_DBSchnsarea T_GMchnsarea_Depth dateofexp
    
    
    
    %%%  avgGMArea %%%
    [data_segments_new, T_chnsarea_new]= avgGMArea(data_segments, T_chnsarea);
    
    data_segments = data_segments_new;
    T_chnsarea = T_chnsarea_new;
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx + tmpn - 1) savefilename_addstr ...
        upper(filename(idx + tmpn)) filename(idx + tmpn + 1:end)];
    
    save(fullfile(savefolder, savefilename), 'data_segments',  'T_chnsarea', 'fs');
end
end




function [data_segments_new, T_chnsarea_new]= avgGMArea(data_segments, T_chnsarea)
%  average across each GM area and keep DBS contact channels 
%

depth_M1Layer5 = [10 14];
depth_PMCLayer5 = [10 14];

depthUsed_M1 = depth_M1Layer5;
depthUsed_PMC = depth_PMCLayer5;



%%% T_GMchnsarea_new %%%

T_GMchnsarea = T_chnsarea(strcmp(T_chnsarea.electype, 'Gray Matter'), :);

% extract uniqGMAreas without "", "GPe" and "GPi" adn DBS contacts
mask_removearea = cellfun(@(x) isempty(x)||strcmp(x, 'GPe')||strcmp(x, 'GPi'), T_GMchnsarea.brainarea);
uniqGMAreas = unique(T_GMchnsarea.brainarea(~mask_removearea));

% T_GMchnsarea_new: note variable used for mask_area
T_GMchnsarea_new = T_GMchnsarea([], :);
for areai = 1:length(uniqGMAreas)
    brainarea = uniqGMAreas{areai};
    
    
    if strcmp(brainarea, 'M1') || strcmp(brainarea, 'PMC')
        
        if strcmp(brainarea, 'M1')
            depthUsed = depthUsed_M1 ;
        else
            depthUsed = depthUsed_PMC ;
        end
        
        % mask_area in T_chnsarea
        mask_area = strcmp(T_chnsarea.brainarea, brainarea) & T_chnsarea.depth >= depthUsed(1) & T_chnsarea.depth <= depthUsed(2);
        
        if ~any(mask_area) % no chnnels in depthUsed
            continue;
        end
        
    else
        % extract the idx of brainarea
        mask_area = strcmp(T_chnsarea.brainarea, brainarea);
    end
    
    
    % T_chnsarea
    T_area = T_chnsarea(mask_area, :);    
    T_area = T_area(1,:);
    T_area.recordingchn = NaN;
    T_area.depth = NaN;
    T_area.notes{1} = mask_area;

    T_GMchnsarea_new = [T_GMchnsarea_new; T_area];
end
T_GMchnsarea_new.chni = [1: height(T_GMchnsarea_new)]';




%%% extract data_segments based T_GMchnsarea_new and  T_DBSchnsarea %%%

T_DBSchnsarea = T_chnsarea(strcmp(T_chnsarea.electype, 'DBS'), :);
for segi = 1:length(data_segments)
    
    lfp_1seg = data_segments(segi).lfp;
    
    % average across each GM area
    avgGMlfp = [];
    for areai = 1:height(T_GMchnsarea_new)
        mask_area = T_GMchnsarea_new.notes{areai};
        
        % extract the avg lfp data of segi in brainarea, avglfp_area: ntemp * 1
        avglfp_area = mean(lfp_1seg(:, mask_area), 2); 
        
        avgGMlfp = cat(2, avgGMlfp, avglfp_area);  
    end
    
    % combine channels from DBS
    lfp_1seg = cat(2, avgGMlfp, lfp_1seg(:, T_DBSchnsarea.chni)); 
    
    data_segments_new(segi).lfp  = lfp_1seg;
    
    
   clear lfp_1seg avgGMlfp
end
T_chnsarea_new = [T_GMchnsarea_new; T_DBSchnsarea];
T_chnsarea_new.chni = [1: height(T_chnsarea_new)]';
end  
    

function T_chnsarea = combine_GMDBSChns(T_GMchnsarea_Depth, T_DBSchnsarea)
% 	combine T_chnsarea_GM_Depth and T_chnsarea_DBS (have depth or no depth variable)
%
%

if ~any(strcmp('depth', T_DBSchnsarea.Properties.VariableNames))
    % add depth variable if T_chnsarea_DBS doesn't have
    
    T_DBSchnsarea = addvars(T_DBSchnsarea, NaN(height(T_DBSchnsarea), 1), 'NewVariableNames', 'depth', 'After', 'recordingchn');
end

T_chnsarea = vertcat(T_GMchnsarea_Depth,T_DBSchnsarea);

% adjust chi
T_chnsarea.chni = [1:height(T_chnsarea)]';

end