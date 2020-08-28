function m6_restData_averagedAcrossOneArea()
%% averaged the lfp in one area
%

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% the corresponding pipeline and the parent folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);


%%  input setup

% input folder: extracted raw rest data with grayMatter
inputfolder = fullfile(codecorresParentfolder, 'm5_restData_segNarrowedDownsample');

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['avg_'];




%% Code Starts Here

files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['Seging Filtered Rest Data....']);


for filei = 1 : nfiles
    disp(files(filei).name)
    
    % wait bar
    waitbar(filei/nfiles,f,['Averaging Rest Data of each area in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data
    filename = files(filei).name;
    load(fullfile(inputfolder, filename), 'lfpsegs_m1', 'lfpsegs_GM', 'lfpsegs_stn', 'lfpsegs_gp','fs', 'GMChnAreas');
    
    
    
    % averaged lfp across M1
    avglfpsegs_m1 = mean(lfpsegs_m1,2);
    
    % extract the unique GM Areas
    uniqGMAreas = [];
    for i = 1: length(GMChnAreas)
        GMArea = string(GMChnAreas{i});
        
        if isempty(uniqGMAreas) || (~isempty(GMArea)&&~contains(GMArea, uniqGMAreas))
            uniqGMAreas = [uniqGMAreas, GMArea];
        end
        
        clear GMArea
    end
    
    % average lfp across each GM area in uniqGMAreas
    avglfpsegs_GM = [];
    for i = 1: length(uniqGMAreas)
        GMArea = uniqGMAreas(i);
        
        % find the channel numbers of GMArea in GMChnAreas(i.e. in lfpsegs_GM)
        indices = cellfun(@(x) strcmp(x, GMArea), GMChnAreas);
        
        % average across GMArea
        avglfpsegs = mean(lfpsegs_GM(:,indices, :),2);
        
        % concatenate to avglfpsegs_GM
        avglfpsegs_GM = cat(2, avglfpsegs_GM, avglfpsegs);
        
        clear GMArea indices avglfpsegs
    end
    
    % concatenate the lfpsegs and chnAreas of all the areas
    lfpsegs = cat(2, avglfpsegs_m1, avglfpsegs_GM, lfpsegs_stn, lfpsegs_gp);
    chnAreas = cellstr(["M1", uniqGMAreas, ...
        "stn0-1", "stn1-2", "stn2-3", "stn3-4", "stn4-5", "stn5-6", "stn6-7", ...
        "gp0-1", "gp1-2", "gp2-3", "gp3-4", "gp4-5", "gp5-6", "gp6-7"]);
    
    
    
     
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpsegs','chnAreas', 'fs');
    
    
        
    clear lfpsegs_stn lfpsegs_gp fs GMChnAreas
    clear avglfpsegs_m1 avglfpsegs_GM uniqGMAreas
    clear lfpsegs chnAreas
end











