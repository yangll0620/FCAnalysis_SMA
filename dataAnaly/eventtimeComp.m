%%
% this script compares the reach and return time between normal and mild
% states using two-sample t-test (ttest2)
%
% Output:
%     figures showing the statistical analysis and autosaved to
%     /Projects/FCAnalysis/results/Pinky/behaviors

%% 
clear

%% add path
addpath(fullfile('..', 'util'))

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
pinkyskbinf = table(dataArray{1:end-1}, 'VariableNames',varNames(1:6));

% Clear temporary variables
clearvars inffolder inffilename inffile delimiter startRow formatSpec fid dataArray varNames;


%% load
eventdir = fullfile(datadrive, 'Animals2', 'Pinky', 'Recording', 'Processed', 'DataDatabase');

% extract the dateofexp, bkma, and bktdt of used skbs which are marked
% 'Yes' in the column 'YingUsed'
validskbs = pinkyskbinf{strcmp(pinkyskbinf.YingUsed, 'Yes'), {'dateofexp','bkma', 'bktdt'}};

for i = 1 :  length(validskbs)
    dateofexp = datenum(validskbs(i, 1), 'yymmdd') ;
    bkma = char(validskbs(i, 2));
    bktdt = char(validskbs(i, 3));
    
    % event folder
    eventfolder = fullfile(eventdir, ['Pinky_' datestr(dateofexp, 'mmddyy')], ['Block-' bktdt]);
    
    % event file
    eventfile = ['pinky_' datestr(dateofexp, 'yyyymmdd') '_' bkma '_cleaned_MA_SingleTargetKluver_Analyze2.mat'];
    
    if exist(fullfile(eventfolder, eventfile)) ~= 2
        disp([eventfile 'does not exist'])
        %continue;
    end
    
    
    % load event .mat file
    load(fullfile(eventfolder, eventfile))
    
    % extract the reach_time and return_time of good trials
    reach_time = SingleTargetKluverMAData.reach_time(find(SingleTargetKluverMAData.goodix_reach == 1));
    return_time = SingleTargetKluverMAData.return_time(find(SingleTargetKluverMAData.goodix_return == 1));
    
    
    % convert to table structure
    tblreach = table(repmat(eventfile, [length(reach_time),1]), reach_time, ...
        'VariableNames', {'filenames', 'reachtime'});
    tblreturn = table(repmat(eventfile, [length(return_time),1]), return_time, ...
        'VariableNames', {'filenames', 'returntime'});
    
    % append return/return time to normal or mild condition
    if strcmp(parsePDCondition_Pinky(dateofexp), 'normal')
        % ---- normal condition
        
        % reach time to normal state
        if exist('tbl_normalreachtime','var') ~= 1
            tbl_normalreachtime = tblreach;
        else
            tbl_normalreachtime = [tbl_normalreachtime; tblreach];
        end
        
        % return time to normal state
        if exist('tbl_normalreturntime','var') ~= 1
            tbl_normalreturntime = tblreturn;
        else
            tbl_normalreturntime = [tbl_normalreturntime; tblreturn];
        end
        
    else
        if strcmp(parsePDCondition_Pinky(dateofexp), 'mild')
            % --- mild condition
            
            % append reach time to mild state
            if exist('tbl_mildreachtime','var') ~= 1
                tbl_mildreachtime = tblreach;
            else
                tbl_mildreachtime = [tbl_mildreachtime; tblreach];
            end
            
            % append return time to mild state
            if exist('tbl_mildreturntime','var') ~= 1
                tbl_mildreturntime = tblreturn;
            else
                tbl_mildreturntime = [tbl_mildreturntime; tblreturn];
            end
            
        end
    end
    
    clear dateofexp bkma bktdt eventfolder eventfile reach_time return_time tblreach tblreturn;
    clear SingleTargetKluverMAData
    
end


%% statistical analysis
stdi = 2; 

% +-2std, 95% normalreachtime
normalreachtime = tbl_normalreachtime{:,2};

valmin = mean(normalreachtime)-stdi* std(normalreachtime);
valmax = mean(normalreachtime)+stdi* std(normalreachtime);
idx = find(normalreachtime>=valmin & normalreachtime <=valmax);

tbl_normalreachtime = tbl_normalreachtime(idx,:);
clear valmin valmax idx normalreachtime;


% +-2std, 95% mildreachtime
mildreachtime = tbl_mildreachtime{:,2};

valmin = mean(mildreachtime)-stdi* std(mildreachtime);
valmax = mean(mildreachtime)+stdi* std(mildreachtime);
idx = find(mildreachtime>=valmin & mildreachtime <=valmax);

tbl_mildreachtime = tbl_mildreachtime(idx,:);
clear valmin valmax idx mildreachtime;

% reach time statistical analysis
[~,p_reachtime, ~, ~] = ttest2(tbl_normalreachtime{:,2}, tbl_mildreachtime{:,2});


% +-2std, 95% normalreturntime
normalreturntime = tbl_normalreturntime{:,2};

valmin = mean(normalreturntime)-stdi* std(normalreturntime);
valmax = mean(normalreturntime)+stdi* std(normalreturntime);
idx = find(normalreturntime>=valmin & normalreturntime <=valmax);

tbl_normalreturntime = tbl_normalreturntime(idx,:);
clear valmin valmax idx normalreturntime;


% +-2std, 95% mildreturntime
mildreturntime = tbl_mildreturntime{:,2};

valmin = mean(mildreturntime)-stdi* std(mildreturntime);
valmax = mean(mildreturntime)+stdi* std(mildreturntime);
idx = find(mildreturntime>=valmin & mildreturntime <=valmax);

tbl_mildreturntime = tbl_mildreturntime(idx,:);
clear valmin valmax idx mildreturntime;

% return time statistical analysis
[~,p_returntime, ~, ~] = ttest2(tbl_normalreturntime{:,2}, tbl_mildreturntime{:,2});


%% plot t-test statistical analysis

%%%%%%   reach time figure   %%%%%%%%%
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box(axes1,'on');

% plot normal and mild reach time
plot(tbl_normalreachtime{:,2},'b.');
plot(tbl_mildreachtime{:,2},'r.');

% legend and title of reach time
legend('normal', 'mild')
title('reach time')

% Create textbox for mean reach time in normal and mild states
annotation(figure1,'textbox', [0.5 0.2 0.3 0.06], ...
    'String',{['normal = '  num2str(mean(tbl_normalreachtime{:,2}))]}, ...
    'LineStyle', 'none');

annotation(figure1,'textbox', [0.5 0.15 0.3 0.06], ...
    'String',{['mild = '  num2str(mean(tbl_mildreachtime{:,2}))]}, ...
    'LineStyle', 'none');


% Create textbox for p-value of reach time
annotation(figure1,'textbox', [0.5 0.8 0.2 0.06], ...
    'String',{['p = '  num2str(p_reachtime)]}, ...
    'LineStyle', 'none');


%%%%%%   return time figure   %%%%%%%%%
figure2 = figure;

% Create axes for return time figure
axes1 = axes('Parent',figure2);
hold(axes1,'on');
box(axes1,'on');

% plot normal and mild return time
plot(tbl_normalreturntime{:,2},'b.');
plot(tbl_mildreturntime{:,2},'r.');

% legend and title of return time figure
legend('normal', 'mild')
title('return time')

% Create textbox for mean return time in normal and mild states
annotation(figure2,'textbox', [0.5 0.2 0.3 0.06], ...
    'String',{['normal = '  num2str(mean(tbl_normalreturntime{:,2}))]}, ...
    'LineStyle', 'none');

annotation(figure2,'textbox', [0.5 0.15 0.3 0.06], ...
    'String',{['mild = '  num2str(mean(tbl_mildreturntime{:,2}))]}, ...
    'LineStyle', 'none');

% Create textbox for p-value of return time
annotation(figure2,'textbox', [0.5 0.8 0.2 0.06], ...
    'String',{['p = '  num2str(p_returntime)]}, ...
    'LineStyle', 'none');

%% save figures

% save folder
savefolder = fullfile('..', '..', '..', ... % /NMRC_umn/
    'Projects', 'FCAnalysis', 'results', 'Pinky', 'behavior');

% save figure1: reach time
print(figure1, fullfile(savefolder, 'reachtime_YYData'),'-dpng')

% save figure2: return time
print(figure2, fullfile(savefolder, 'returntime_YYData'),'-dpng')

% disp information
disp(['Figures are saved to ' fullfile(savefolder, 'reachtime_YYData') ' and returntime_YYData'])