function skbinf_extract()
% extract the date, bkma and tdtma information for all the single kluver
% board (SKB) task from the google master sheet.
% 
% Output:
%       save the date, bkma and tdtma information for all SKB tasks into
%       F:\yang7003@umn\NMRC_umn\Projects\FCAnalysis\metainf\Pinky\pinky_skbdatebk.csv file

%
addpath(fullfile('..','..','util'))
animal = 'Pinky';
googlesheetlink_Pinky = '1Mn_HcvWt4FVc2kcvRMbDj5jQffyGNrwppNOK6T5W-zI';

%% get the google sheet
googlesheetdata = GetGoogleSpreadsheet(googlesheetlink_Pinky);

%% variable name for each column in google sheet
% Brief Description column record task discription
varname_taskdisp = 'Brief Description'; 
varname_tdtbk = 'TDT Block';
varname_mabk = 'MA File';
varname_date = 'OutputFolderName';
% different discriptions of Kluver board task
taskdisps = {'Kluver', 'Single Target Kluver', 'SKB', 'Single KB'  ...
    'Single Kluver', 'Single', 'Single ' ,'Single Target', ...
    'Single Target Kluver', 'Single-Target Kluver', 'Single-target Kluver',...
    'single', 'single Kluver', 'single-target Kluver'};

%%  parse the exp date, bkma, bktdt for single kluverboard task
% varname_googlesheet: column name of google sheet
varname_googlesheet = googlesheetdata(1,:);

% extract the column number for 'Brief Description', 'tdtbk' and 'date'
col_match = @(varname) find(cell2mat(cellfun(@(x) contains(x, varname),varname_googlesheet,'UniformOutput', 0)));
coli_taskdisp = col_match(varname_taskdisp);
coli_tdtbk = col_match(varname_tdtbk);
coli_date =  col_match(varname_date);
coli_mabk = col_match(varname_mabk);

% extract the task discription
strs_taskdispinsheet = googlesheetdata(1:end,coli_taskdisp);
% define abstract matching function
row_match = @(strpattern) find(cell2mat(cellfun(@(x) ~isempty(find(strcmp(strpattern, x))), strs_taskdispinsheet, 'UniformOutput',0)));
% extract all the row numbers for single kluver board task: idxrow
idxrow = row_match(taskdisps); 
clear strs_taskdispinsheet

% extract the date, MAblock and TDTblock for each single kluver board task
nrows = length(idxrow);
skbtable = table('Size',[nrows, 3],'VariableTypes',{'string', 'uint8', 'uint8'}, ...
    'VariableNames',{'dateexp','bkma', 'bktdt'});
for i = 1:nrows
    rowi = idxrow(i);
    
    strs = split(googlesheetdata{rowi, coli_date},'_');
    skbtable.dateexp(i) = datestr(datevec(strs{2},'mmddyy'),'yymmdd');
    skbtable.bkma(i) = str2num(googlesheetdata{rowi,coli_mabk});
    skbtable.bktdt(i) = str2num(googlesheetdata{rowi,coli_tdtbk});
    
    clear rowi strs
end

%% save
savefolder = ['F:\yang7003@umn\NMRC_umn\Projects\FCAnalysis\metainf\' animal];
savename_skbtbl = [lower(animal) '_skbinf.csv'];
writetable(skbtable, fullfile(savefolder, savename_skbtbl));
