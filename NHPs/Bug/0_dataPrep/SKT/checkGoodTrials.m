clear
mafolder = 'C:\Users\nmrc3\Desktop\Bug_062818\Block-3';

mafilestruct = dir(fullfile(mafolder, '*SingleTargetKluver_Analyze2.mat'));

% ma file does not exist
if isempty(mafilestruct)
    disp([mafolder 'has no Analyze2.mat'])
   
end


% load SingleTargetKluverMAData
load(fullfile(mafolder, mafilestruct.name), 'SingleTargetKluverMAData'); 


% the tag of good reach trials
tag_goodreach = SingleTargetKluverMAData.goodix_reach;
% the tag of good return trials
tag_goodreturn = SingleTargetKluverMAData.goodix_return;

% extract indices of good trials (both have good reach and return)
idx_goodtrials = find(tag_goodreach .* tag_goodreturn == 1);

disp(['number of good trials = ' num2str(length(idx_goodtrials))])


% get the pd conditioon for the date of experiment
tmps = split(mafolder, '_');
dateofexp = datenum(tmps{2}(1:6), 'mmddyy');
pdcondition = parsePDCondition(dateofexp, 'Bug');

% ma sample rate
fs_ma = SingleTargetKluverMAData.SR;

% time indices for target onset, reach onset, touch screen, return and mouth
TargetTime = SingleTargetKluverMAData.TargetTime;
ReachTimeix = SingleTargetKluverMAData.ReachTimeix;
TouchTimeix = SingleTargetKluverMAData.TouchTimeix;
ReturnTimeix = round(SingleTargetKluverMAData.ReturnTimeix); % not integer in .ReturnTimeix
MouthTimeix = SingleTargetKluverMAData.MouthTimeix;

timeixtbl_ma = [table(TargetTime) table(ReachTimeix) table(TouchTimeix) table(ReturnTimeix) table(MouthTimeix)];
timeixtbl_ma = timeixtbl_ma(idx_goodtrials,:);


onedayReachTimes = (timeixtbl_ma.TouchTimeix - timeixtbl_ma.ReachTimeix)/ fs_ma;
onedayReturnTimes = (timeixtbl_ma.MouthTimeix - timeixtbl_ma.ReturnTimeix)/ fs_ma;
disp([datestr(dateofexp, 'yyyymmdd') pdcondition ' avg reach time = ' num2str(mean(onedayReachTimes)) ...
    ', return time = ' num2str(mean(onedayReturnTimes))])

