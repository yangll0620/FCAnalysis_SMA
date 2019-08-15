%% 
clear

%% parse datadrive and googledrive base on operation system
if ispc
    datadrive = 'Y:/root2';
    googledrive = '';
end
if isunix
    datadrive = '/home/lingling/root2';
    googledrive = '/home/lingling/yang7003@umn.edu/NMRC_umn';
end

%% Read skb information in /Projects/FCAnalysis/metainf/pinky_skbinf.csv
inffolder = fullfile(googledrive, 'Projects', 'FCAnalysis', 'metainf', 'Pinky');
inffilename = 'pinky_skbinf.csv';
inffile = fullfile(inffolder, inffilename);

% open the text file
fid = fopen(inffile, 'r');

% extract the column name
varNames = split(fgetl(fid), ',');

% Format for each line of text:
%   column1: test (%s)
%	column2: int8 (%d8)
%   column3: int8 (%d8)
%	column4: categorical (%C)
%   column5: categorical (%C)
%	column6: categorical (%C)
formatSpec = '%s%d8%d8%s%s%s%[^\n\r]';
delimiter = ',';
dataArray = textscan(fid, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false, 'EndofLine', '\r\n');

fclose(fid);

% Create output variable
pinkyskbinf = table(dataArray{1:end-1}, 'VariableNames', varNames);

% Clear temporary variables
clearvars inffolder inffilename inffile delimiter startRow formatSpec fid dataArray varNames;


%% load
eventdir = fullfile(datadrive, 'Animals2', 'Pinky', 'Recording', 'Processed', 'DataDatabase');

% extract the dateofexp, bkma, and bktdt of used skbs which are marked
% 'Yes' in the column 'YingUsed'
validskbs = pinkyskbinf{strcmp(pinkyskbinf.YingUsed, 'Yes'), {'dateofexp','bkma', 'bktdt'}};

% load all the reach and return time
for i = 1 : length(validskbs)
    dateofexp = datenum(validskbs(i, 1), 'yymmdd') ;
    bkma = char(validskbs(i, 2));
    bktdt = char(validskbs(i, 3));
    
    % event folder
    eventfolder = fullfile(eventdir, ['Pinky_' datestr(dateofexp, 'mmddyy')], ['Block-' bktdt]);
    
    % event file
    eventfile = ['pinky_' datestr(dateofexp, 'yyyymmdd') '_' bkma '_cleaned_MA_SingleTargetKluver_Analyze2.mat'];
    
    if exist(fullfile(eventfolder, eventfile)) ~= 2
        disp([eventfile 'does not exist'])
        continue;
    end
    
    
    % load event .mat file
    load(fullfile(eventfolder, eventfile))
    
    reach_time = SingleTargetKluverMAData.reach_time(find(SingleTargetKluverMAData.goodix_reach == 1));
    return_time = SingleTargetKluverMAData.return_time(find(SingleTargetKluverMAData.goodix_return == 1));
    
    % append reach time 
    tblreach = table(repmat(eventfile, [length(reach_time),1]), reach_time, ...
        'VariableNames', {'filenames', 'reachtime'});
    if exist('tbl_reachtime','var') ~= 1
        tbl_reachtime = tblreach;
    else
        tbl_reachtime = [tbl_reachtime; tblreach];
    end
    
    % append return time
    tblreturn = table(repmat(eventfile, [length(return_time),1]), return_time, ...
        'VariableNames', {'filenames', 'returntime'});
    if exist('tbl_returntime','var') ~= 1
        tbl_returntime = tblreturn;
    else
        tbl_returntime = [tbl_returntime; tblreturn];
    end

    
    clear dateofexp bkma bktdt eventfolder eventfile reach_time return_time tblreach tblreturn;
    clear SingleTargetKluverMAData
    
end

