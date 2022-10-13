clear 

%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\pipeline\NHPs\JoKittyGoodBadTrialDate.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 19-Sep-2022 14:45:45

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 5);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:E49";

% Specify column names and types
opts.VariableNames = ["Animal", "State", "Date", "bktdt", "VarName5"];
opts.VariableTypes = ["categorical", "categorical", "double", "categorical", "categorical"];

% Specify variable properties
opts = setvaropts(opts, ["Animal", "State", "bktdt", "VarName5"], "EmptyFieldRule", "auto");

% Import the data
xlsfile = "H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\JoKittyGoodBadTrialDate.xlsx";
JoKittyGoodBadTrialDate = readtable(xlsfile, opts, "UseExcel", false);


%% Clear temporary variables
clear opts

%% code Start here
folder_server_Jo = fullfile('Z:', 'root2', 'Animals', 'Jo', 'Recording', 'Processed', 'DataDatabase');
folder_server_Kitty = fullfile('Z:', 'root2', 'Animals', 'Kitty', 'Recording', 'Processed', 'DataDatabase');
mafolders = {}; 
mafiles = {};
for ti = 1 : height(JoKittyGoodBadTrialDate) 
    animal = char(JoKittyGoodBadTrialDate.Animal(ti));
    dateofexp = datevec(num2str(JoKittyGoodBadTrialDate.Date(ti)), 'yyyymmdd');
    tmp = char(JoKittyGoodBadTrialDate.bktdt(ti));
    bktdt = str2num(tmp(6:end));
    pdcond = char(JoKittyGoodBadTrialDate.State(ti));
    
    if strcmpi(animal, 'Jo')
        processedfolder_inroot2 = folder_server_Jo;
    end
    
    if strcmpi(animal, 'Kitty')
        processedfolder_inroot2 = folder_server_Kitty;
    end
    
    
     
    onedaypath = fullfile(processedfolder_inroot2, [animal '_' datestr(dateofexp, 'mmddyy')]);
    mafolder = fullfile(onedaypath, ['Block-' num2str(bktdt)]);
    
    if strcmpi(animal, 'Kitty') && strcmpi(pdcond, 'normal')
        mafilestruct = dir(fullfile(mafolder, '*COT_Analyze2.mat'));
    else
        mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));
    end
    
    if length(mafilestruct) > 1
        disp([mafolder 'has more than 1 ' num2str(length(mafilestruct)) ' files!'])
        continue;
    end
    
    if isempty(mafilestruct)
        disp([mafolder ' has ' num2str(length(mafilestruct)) ' files!'])

        mafolders = [mafolders; mafolder];
        mafiles = [mafiles; 'empty'];
    else
              
        mafolders = [mafolders; mafolder];
        mafiles = [mafiles; fullfile(mafolder, mafilestruct(1).name)];
    end
    
    
    
    
    clear animal dateofexp tmp bktdt 
    clear processedfolder_inroot2 onedaypath mafolder mafilestruct
end

JoKittyGoodBadTrialDate = [JoKittyGoodBadTrialDate cell2table(mafolders) cell2table(mafiles)];

% save
writetable(JoKittyGoodBadTrialDate,xlsfile,'Sheet',1)

