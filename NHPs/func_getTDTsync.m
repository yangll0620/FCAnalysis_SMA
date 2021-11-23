function [ SyncTDT ] = func_getTDTsync( filename, plotON )
%func_getTDTsync Accesses TDT block sync (startpad) signal


%% Author Luke Johnson created on 5/19/2016
% change log: 

%%
% INPUT
% filename: full path of TDT block you want to open Y:\Luke\Jo\Recording\Raw\rawTDT\Jo_CR1_DT1_020216\Block-1
% plotON: 1 or 0, 1 if you want to plot the sync (startpad) data

% OUTPUT

% SyncTDT = 
% 
%        sync: [8726922x1 single] sync signal of 1's and 0's
%      sync_t: [8726922x1 double]  time in seconds
%     sync_SR: 2.4414e+04 5000 sampling rate

%     header.MonkeyName     =  'Jo'; % Define which animal
%     header.server_dir = 'Y:\'; %enter appropriate directory on your PC for the server
%     header.TDTSessionName       =   txt(k+1, ismember(heads,'TDT SessionName')==1); % e.g. '\Jo_CR1_DT1_020116'
%     header.TDTBlockName         =   cell2mat(raw(k+1, ismember(heads,'TDT Block')==1)); % e.g. 5
%     header.rawTDTFolder         =   txt(k+1, ismember(heads,'rawTDTfolder')==1); % e.g. 'Y:\Luke\Jo\Recording\Raw\rawTDT'
    
    
%% Access Data Tank
%     filename = [header.server_dir header.rawTDTFolder{files_ix} header.TDTSessionName{files_ix} '\Block-' num2str(header.TDTBlockName(files_ix))]
[ TTX, FileInfoTD ] = func_AccessTDTDataTanks( filename )
 

 EventName = 'StPd'; % Define array variable event name
ix = find(strcmp(FileInfoTD.EventName,EventName));              
SR = FileInfoTD.SampleFrequency(ix);
TTX.SetGlobalV('Channel',1);
StPd = TTX.ReadWavesV(EventName);

threshold = 0.3; 
StPd(StPd<threshold) = 0;
StPd(StPd>=threshold) = 1;
    
tStPd = [0:length(StPd)-1]'./SR;
if plotON
    figure;plot(tStPd,StPd)
end


SyncTDT.sync = StPd;
SyncTDT.sync_t = tStPd;
SyncTDT.sync_SR = SR;


TTX.CloseTank
TTX.ReleaseServer
close all
    
end