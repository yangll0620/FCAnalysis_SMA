% Created by Luke Johnson 05/20/2016
% update log:

%06/01/2016 LAJ If passive movement times output file already exists, confirm that you want to overwrite
%06/01/2016 LAJ Add jiggler single marker functionality
% Adapted from run_GetPassiveMovementTimes for kluver board reaching
%01/19/2017 Ying Add value of wristix in line 229&237.

% run_GetSingleTrialKluverMovementTimes.m  is used to load marker data from defined MA
% files and identify the active movement times target, reach, return etc
% The sync signal is loaded from MA and TDT and the movement times
% are adjusted to match the TDT physiology recordings

% Resultant movement times are saved in a .mat file and .nex file


%% add path YLL
addpath(genpath(fullfile('H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\code', 'toolbox', 'NexMatlabFiles')))

%% clear variables
close all; clear all; 
clc;
clear all

%% START USER INPUT
% Initialize parameters

header.MonkeyName     =  'Kitty'; % Define which animal
header.server_dir = 'H:\My Drive\NMRC_umn\Projects\FCAnalysis\exp\data\'; %enter appropriate directory on your PC for the server
% SessionsBlocksToAnalyze{1}.Session = 'Jo_011116' % Corresponds to OutputFolderName in spreadsheet
SessionsBlocksToAnalyze{1}.Session = 'Kitty_041515' % Corresponds to OutputFolderName in spreadsheet
% SessionsBlocksToAnalyze{1}.Session = 'Pinky_030217' % Corresponds to OutputFolderName in spreadsheet

% SessionsBlocksToAnalyze{1}.Blocks = [2 5 8 10 12 14 18] % the blocks you wish to process
SessionsBlocksToAnalyze{1}.Blocks = [1] % %NOTE THAT BLOCKS REFERS TO THE MA FILES, NOT NECESSARILY THE TDT BLOCKS -- USUALLY THE SAME BUT NOT ALWAYS!!!
% SessionsBlocksToAnalyze{1}.Blocks = [11] % %NOTE THAT BLOCKS REFERS TO THE MA FILES, NOT NECESSARILY THE TDT BLOCKS -- USUALLY THE SAME BUT NOT ALWAYS!!!

% SessionsBlocksToAnalyze{2}.Session = 'Jo_020416'
% SessionsBlocksToAnalyze{2}.Blocks = [1 3]



%     header.MonkeyName     =  'Kitty'; % Define which animal
%     header.server_dir = 'X:\'; %enter appropriate directory on your PC for the server
% SessionsBlocksToAnalyze{1}.Session = 'Kitty_011515' % Corresponds to OutputFolderName in spreadsheet
% % SessionsBlocksToAnalyze{1}.Session = 'Kitty_121714' % Corresponds to OutputFolderName in spreadsheet
% % SessionsBlocksToAnalyze{1}.Session = 'Kitty_010815' % Corresponds to OutputFolderName in spreadsheet
% 
% % SessionsBlocksToAnalyze{1}.Blocks = [2 5 8 10 12 14 18] % the blocks you wish to process
% SessionsBlocksToAnalyze{1}.Blocks = [20 28 30] % %NOTE THAT BLOCKS REFERS TO THE MA FILES, NOT NECESSARILY THE TDT BLOCKS -- USUALLY THE SAME BUT NOT ALWAYS!!!
% % SessionsBlocksToAnalyze{2}.Session = 'Jo_020416'
% % SessionsBlocksToAnalyze{2}.Blocks = [1 3]
plotON = 0;
%% END USER INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Access Speadsheet Database
if strcmp(header.MonkeyName, 'Jo')
    
    % Select spreadsheet and worksheet to read from your google drive
    spreadsheetName='JoMasterDatabase';
    sheetName = 'MasterList' %choose appropriate excel sheet
    web_link ='1DtQo5qMbkvTDDQlHdPUcFMIOWUP09pZMiZen0t6yEXM';  %% it would be different for other user
elseif strcmp(header.MonkeyName, 'Pinky')
    web_link ='1Mn_HcvWt4FVc2kcvRMbDj5jQffyGNrwppNOK6T5W-zI';  %% it would be different for other user
    spreadsheetName='PinkyMasterDatabase';
    sheetName = 'MasterList' %choose appropriate excel sheet
elseif strcmp(header.MonkeyName, 'Bug')
    web_link ='1mvTwDnWtN2-UhF4o_NSC9NjF1681lIss9zI1vn48GlU';  %% it would be different for other user
    spreadsheetName='BugMasterDatabase';
    sheetName = 'MasterList' %choose appropriate excel sheet
elseif strcmp(header.MonkeyName, 'Kitty')
    
    % Select spreadsheet and worksheet to read from your google drive
    spreadsheetName='KittyMasterDatabase';
    sheetName = 'MasterList' %choose appropriate excel sheet
    web_link ='1B4lnBGlYiI-O_tZWEz9Thvj3hVkczh7iWNS7mQ3PNsc';  %% it would be different for other user
end
    
    % Access the spreadsheet i.e. google drive worksheet
    result = GetGoogleSpreadsheet(web_link);
    [raw, txt] = rawToString(result );
    heads = txt(1,:);


    %% Find files in spreadsheet to analyze based on Session and Block information provided
    k = [];
    SessionList    = (raw(:,ismember(heads,'OutputFolderName')==1)); %find excel rows that have been marked for processing
    BlockList    = (raw(:,ismember(heads,'MA File')==1)); %find excel rows that have been marked for processing
    BlockList(1) = [];
    BlockList = cell2mat(BlockList); %stupid workaround
 
    for i = 1:length(SessionsBlocksToAnalyze)
        for j = 1:length(SessionsBlocksToAnalyze{i}.Blocks)
            
            isSession = strcmp(SessionsBlocksToAnalyze{i}.Session, SessionList); isSession(1) = [];
            isBlock = (SessionsBlocksToAnalyze{i}.Blocks(j)== BlockList)
            k = cat(1,k,find(isSession&isBlock))
            
        end
    end
%%

% define information necessary for loading TDT data (Sesssion, Block etc
% based on spreadsheet)
header.TDTSessionName       =   txt(k+1, ismember(heads,'TDT SessionName')==1); % e.g. '\Jo_CR1_DT1_020116'
header.TDTBlockName         =   cell2mat(raw(k+1, ismember(heads,'TDT Block')==1)); % e.g. 5
header.rawTDTFolder         =   txt(k+1, ismember(heads,'rawTDTfolder')==1); % e.g. 'Y:\Luke\Jo\Recording\Raw\rawTDT'
% output folder name
header.DataOutputDir      =   txt(k+1, ismember(heads,'OutputFolderDirectory')==1); % e.g. '\Jo_CR1_DT1_020116'
header.DataOutputFolder      =   txt(k+1, ismember(heads,'OutputFolderName')==1);
    
    
% define information necessary for loading MA data (Session, File etc based on spreadsheet)
header.MASessionFolder       =   txt(k+1, ismember(heads,'MA SessionFolder')==1); % e.g. '\MA20160202'
header.MASessionName       =   txt(k+1, ismember(heads,'MA SessionName')==1); % e.g. 'Jo_20160202'
header.MAFileName         =   cell2mat(raw(k+1, ismember(heads,'MA File')==1)); % e.g. 5
header.MACleanedFileName         =   txt(k+1, ismember(heads,'MA CleanedFileName')==1); % '_cleaned'
header.rawMAFolder         =   txt(k+1, ismember(heads,'rawMAfolder')==1); % e.g. 'Y:\Luke\Jo\Recording\Raw\rawMA'

header.nMovements       = cell2mat(raw(k+1, ismember(heads,'nMovements')==1)); %
header.passiveOrder        = raw(k+1, ismember(heads,'passiveOrder')==1); %
try
    header.passiveCondition        = raw(k+1, ismember(heads,'passiveCondition')==1); %
catch
    header.passiveCondition = {};
end
%%
for files_ix = 1:length(k)
        
    % make output folder if not already there (e.g. may have been created
    % during spike processing or other processing)
    savedir_session =  [header.server_dir header.DataOutputDir{files_ix} header.DataOutputFolder{files_ix}];
    savedir_block = [savedir_session '\Block-' num2str(header.TDTBlockName(files_ix))];
    mkdir(savedir_block);
    
    
    %% If kluver movement times output file already exists, confirm that you want to overwrite
    %.mat file
    savefilename = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_SingleTargetKluver'];
    savefilename1 = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_Analyze'];
    
    if exist([savefilename '.mat']) == 2 ||(exist([savefilename1 '.mat']))
        fileinfo = [header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_SingleTargetKluver.mat']
        output = overwrite_file(fileinfo)
        if output==1
            disp(['overwriting file:  ' fileinfo])
        elseif output ==0
            disp('go to next iteration in for loop')
            continue % go to next iteration in for loop
        end
    else %do nothing, 
    end
    
    
    %%
    
%%     get MA Marker Data
%      load .mat file if already saved, otherwise get from .trc file

    % presumed name of .mat MA_MarkerData filename if already saved
%     assuming saved marker data ends in _cleaned_MA_MarkerData_wPadding.mat meaning
%     that marker data has already been shifted according to TDT&MA sync
%     signals and padded.
    savefilename = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_MarkerData_wPadding']; 
    
    try        
        load(savefilename) % variable loaded should be MarkerDataMA  
        
    catch         
        Padding = 1;    
        % Load MA Sync Data
        % load .mat file if already saved, otherwise get from .anc file 
        % presumed name of .mat _MA_StPd filename if already saved
        StPdfilename = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) '_MA_StPd']; 
        if exist([StPdfilename '.mat'])==2        
            load(StPdfilename) % variable loaded should be MarkerDataMA
            SyncMA.sync = Sync; SyncMA.sync_t = [0:length(Sync)-1]'./SyncSR; SyncMA.sync_SR = SyncSR;
        else        
            filename = [header.server_dir header.rawMAFolder{files_ix} header.MASessionFolder{files_ix} '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) '.anc']
            SyncMA = func_getMAsync(filename, 0);

            Sync = SyncMA.sync; SyncSR = SyncMA.sync_SR;
            save([StPdfilename],'Sync','SyncSR')
        end

    % Load TDT Sync Data
         filename = [header.server_dir header.rawTDTFolder{files_ix} header.TDTSessionName{files_ix} '\Block-' num2str(header.TDTBlockName(files_ix))]
        [SyncTDT] = func_getTDTsync( filename, 0 )

    % upsample MA and identify MA correction
        SyncMA.sync_interp = interp1(SyncMA.sync_t,SyncMA.sync,[0:1/24414.0625:SyncMA.sync_t(end)],'linear');
        D = finddelay(SyncTDT.sync,SyncMA.sync_interp);
        MA_Correction = -D/24414.0625 % almost always should be positive since MA is typically started after TDT system)

        if plotON
            figure
            subplot(2,1,1)
            plot(SyncMA.sync_t +MA_Correction,SyncMA.sync,'r.-'); hold on;
            plot(SyncTDT.sync_t,SyncTDT.sync,'b--')
            findix = find(SyncMA.sync>.1);
            set(gca,'Xlim',[SyncMA.sync_t(findix(1)-20) SyncMA.sync_t(findix(1)+20)]+MA_Correction)
            title('beginning')
            subplot(2,1,2)
            plot(SyncMA.sync_t  +MA_Correction,SyncMA.sync,'r.-'); hold on;
            plot(SyncTDT.sync_t,SyncTDT.sync,'b--')
            findix = find(SyncMA.sync>.1);
            set(gca,'Xlim',[SyncMA.sync_t(findix(end)-20) SyncMA.sync_t(findix(end)+20)]+MA_Correction)
            xlabel('time (s)')
            title('end')
        end    

        PaddingInfo.TDTtime = SyncTDT.sync_t; PaddingInfo.TDTSR = SyncTDT.sync_SR;
        PaddingInfo.downsamplefactor = 240; % new sampling rate of marker data output from func_getMAMarkerData will be PaddingInfo.TDTSR/PaddingInfo.downsamplefactor
        PaddingInfo.MA_Correction = MA_Correction;
        
        % get MA Marker Data
        filename = [header.server_dir header.rawMAFolder{files_ix} header.MASessionFolder{files_ix} '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '.trc']
        MarkerDataMA = func_getMAMarkerData_yy(filename, plotON, PaddingInfo ); 
        
    % save MA Marker Data
    
        MAMarkerDatasavefilename = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_MarkerData_wPadding']; 
        %output of func_getMAMarkerData contains redundant information, will
        %not save all of that data.  Dont save time vector because that can be
        %created from the length of data and sampling rate
        MarkerDataMA = rmfield(MarkerDataMA,{'ByMarker','MarkerData_t','MarkerData_pad_t','ByMarker_pad'})
        save(MAMarkerDatasavefilename,'MarkerDataMA')
    end
        
%% identify number of markers used
    
    t = [0:length(MarkerDataMA.MarkerData_d(:,1))-1]./MarkerDataMA.MarkerData_SR;
    if length(MarkerDataMA.MarkerData_headers) == 3  %figure out how many markers were used
        header.nMarkers(files_ix)= 1;  
         wristix = 1:3;
    elseif length(MarkerDataMA.MarkerData_headers) == 6  %figure out how many markers were used
        header.nMarkers(files_ix)= 2;
        wristix = 1:3;
        kluverix = 4:6;

    elseif length(MarkerDataMA.MarkerData_headers) == 9  %figure out how many markers were used
        header.nMarkers(files_ix)= 3; 
         wristix = 4:6;
         kluverix = 7:9;
    elseif length(MarkerDataMA.MarkerData_headers) == 12  %figure out how many markers were used
        header.nMarkers(files_ix)= 4; 
        wristix = 7:9;
        kluverix = 10:12;
    end
    
    %% identify order of passive movements assessed
%     if header.nMovements(files_ix) == 1 
%         passiveOrder{1} = header.passiveOrder{files_ix}(1:5)
%     elseif header.nMovements(files_ix) == 2 
%         passiveOrder{1} = header.passiveOrder{files_ix}(1:5);
%         passiveOrder{2} = header.passiveOrder{files_ix}(8:12)
%     elseif header.nMovements(files_ix) == 3 
%         passiveOrder{1} = header.passiveOrder{files_ix}(1:5);
%         passiveOrder{2} = header.passiveOrder{files_ix}(8:12);
%         passiveOrder{3} = header.passiveOrder{files_ix}(15:19)
%     elseif header.nMovements(files_ix) == 4 % 
%         passiveOrder{1} = header.passiveOrder{files_ix}(1:5);
%         passiveOrder{2} = header.passiveOrder{files_ix}(8:12);
%         passiveOrder{3} = header.passiveOrder{files_ix}(15:19);
%         passiveOrder{4} = header.passiveOrder{files_ix}(22:26)
%     elseif header.nMovements(files_ix) == 5 % 
%         passiveOrder{1} = header.passiveOrder{files_ix}(1:5);
%         passiveOrder{2} = header.passiveOrder{files_ix}(8:12);
%         passiveOrder{3} = header.passiveOrder{files_ix}(15:19);
%         passiveOrder{4} = header.passiveOrder{files_ix}(22:26);
%         passiveOrder{5} = header.passiveOrder{files_ix}(29:33)
%     elseif header.nMovements(files_ix) == 6 % 
%         passiveOrder{1} = header.passiveOrder{files_ix}(1:5);
%         passiveOrder{2} = header.passiveOrder{files_ix}(8:12);
%         passiveOrder{3} = header.passiveOrder{files_ix}(15:19);
%         passiveOrder{4} = header.passiveOrder{files_ix}(22:26);
%         passiveOrder{5} = header.passiveOrder{files_ix}(29:33);
%         passiveOrder{6} = header.passiveOrder{files_ix}(36:40)
%     end
%     
%     passiveOrder = {};
%     passiveOrder_temp = header.passiveOrder{files_ix};
%     findcomma = strfind(passiveOrder_temp,',')
%     indices = cat(2,0,findcomma,length( passiveOrder_temp)+1)
%     for xx = 1:length(findcomma)+1
%         passiveOrder{xx} = strtrim(passiveOrder_temp(indices(xx)+1:indices(xx+1)-1));
%     end
    
    
    %% begin creating .nex data structure for passive movement times
        dnex_SingleKluver = [];
        dnex_SingleKluver.comment = 'Single Target Kluver movement data aligned to TDT';
        dnex_SingleKluver.version = 101;  
        A = exist('SyncTDT', 'var')
        if ~A % Load TDT Sync so to get SR and timepoints for file
            filename = [header.server_dir header.rawTDTFolder{files_ix} header.TDTSessionName{files_ix} '\Block-' num2str(header.TDTBlockName(files_ix))]
            [SyncTDT] = func_getTDTsync( filename, 0 )
        end
        dnex_SingleKluver.freq = SyncTDT.sync_SR; %NOT MA SR, RATHER MAKE SURE THIS MATCHES THE SU DATA SAMPLING RATE, which is important when merging files later on.
        dnex_SingleKluver.tbeg = 0;
        dnex_SingleKluver.tend = SyncTDT.sync_t(end);  
        dnex_SingleKluver = nexAddInterval(dnex_SingleKluver,dnex_SingleKluver.tbeg,dnex_SingleKluver.tend,'All');
        
%     %% passive movement conditions defined in spreadsheet (i.e. Pre, DBS, Post)
%     if isempty(header.passiveCondition);
%         passiveCondition = 'not specified';
%     else
%         passiveCondition_temp = header.passiveCondition{files_ix};
%         findcomma = strfind(passiveCondition_temp,',')
%         indices = cat(2,0,findcomma,length( passiveCondition_temp)+1)
%         for xx = 1:length(findcomma)+1
%             passiveCondition{xx} = strtrim(passiveCondition_temp(indices(xx)+1:indices(xx+1)-1));
%         end
%     end
    
%     %% extract passive movement timepoints (e.g. max flexion, max extension)
%     passiveTimestamps = {};
%     for j = 1:length(passiveOrder)
%         j
%         passiveType = passiveOrder{j}; nMarkers = header.nMarkers(i);        
%         output = func_GetPassiveMovements(nMarkers, passiveType, passiveOrder, MarkerDataMA.MarkerData_pad_d, MarkerDataMA.MarkerData_pad_SR )
%         close all;
%         
% %         % hand exceptions
% %         fnamecompare = [header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_MarkerData']; 
% %         if strcmp(fnamecompare,'Jo_20160201_2_cleaned_MA_MarkerData.mat') 
% %             output.t_ext(output.t_ext>82.5&output.t_ext<85.15)=[]
% %             output.t_flx(output.t_flx>82.5&output.t_flx<85.15)=[]
% %             output.t_ext(output.t_ext>104&output.t_ext<104.8)=[]
% %             output.t_flx(output.t_flx>104&output.t_flx<104.8)=[]
% %             output.t_ext(output.t_ext>115.2&output.t_ext<117)=[]
% %             output.t_flx(output.t_flx>115.2&output.t_flx<117)=[]
% %         end
% %         % hand exceptions
% %         if strcmp(fnamecompare,'Jo_20151009_4_cleaned_MA_MarkerData.mat') && j==6 % drift in baseline, do twice and combine
% %             output2 = func_extractPassiveMovements(header.nMarkers(i),passiveType,d)
% %             output.t_ext = cat(1,output.t_ext,output2.t_ext)
% %             output.t_flx = cat(1,output.t_flx,output2.t_flx)
% % 
% %         end
% 
%         passiveTimestamps{j} = output;
%         passiveTimestamps{j}.transformed_joint_pos_SR = MarkerDataMA.MarkerData_pad_SR;
%         
%         %% add to nex file structure
%         % Add passive movement and times to .nex file  
%         t_ext = passiveTimestamps{j}.t_ext;
%         t_flx = passiveTimestamps{j}.t_flx;
%         if isempty(header.passiveCondition);
%             t_ext_name = [passiveTimestamps{j}.joint '_ext'];
%             t_flx_name = [passiveTimestamps{j}.joint '_flx'];
%             joint_pos_name = [passiveTimestamps{j}.joint '_TransformedMarkerPos'];
%         else
%             t_ext_name = [passiveTimestamps{j}.joint '_ext_' passiveCondition{j}];
%             t_flx_name = [passiveTimestamps{j}.joint '_flx_' passiveCondition{j}];
%             joint_pos_name = [passiveTimestamps{j}.joint '_TransformedMarkerPos_' passiveCondition{j}];
%         end
% 
%         
%     	dnex_SingleKluver = nexAddEvent(dnex_SingleKluver, t_ext, t_ext_name);
%     	dnex_SingleKluver = nexAddEvent(dnex_SingleKluver, t_flx, t_flx_name);
%         dnex_SingleKluver = nexAddContinuous(dnex_SingleKluver, 0, MarkerDataMA.MarkerData_pad_SR, output.transformed_joint_pos, joint_pos_name);
% 
%     end
    

%% Right now assuming SingleTargetKluver.m has been run and data saved as for example Jo20151009_6_cleaned_Analyze.mat and copied into appropriate DataDatabase folder
%Future will put that SingleTargetKluver.m analysis into here


%     savefilename = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_cleaned_Analyze']; 
%     temp = load(savefilename)
%     temp.Analyze
      
%% 
%% identify typical start position (semi-manually)

    savefilename = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_SingleTargetKluver_Analyze1']; 

    extract_behavior_timepoints = [];
    if exist([savefilename '.mat']) == 2 || (exist([savefilename1 '.mat']) == 2 )
        fileinfo = [header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_SingleTargetKluver_Analyze1']
        output = overwrite_file(fileinfo)
        if output==1
            disp(['overwriting file:  ' fileinfo])
            extract_behavior_timepoints = 1;
        elseif output ==0
            disp('not overwriting file')
            extract_behavior_timepoints = 0;
        end
    else
       extract_behavior_timepoints = 1; 
    end
    
    if extract_behavior_timepoints == 1
        
        Analyze = [];
        Analyze.SR = MarkerDataMA.MarkerData_pad_SR;
        Analyze.Wrist = MarkerDataMA.MarkerData_pad_d(:,wristix);
        Analyze.Kluver = MarkerDataMA.MarkerData_pad_d(:,kluverix);

        scrsz = get(groot,'ScreenSize');
        figure('Position',[1 1 scrsz(3)/.99 scrsz(4)/.99])
        tempstart = 5000;
        subplot(2,1,1)
        plot(Analyze.Wrist(tempstart:tempstart+10000,:));
        title('SELECT START THEN END OF REACH')
        [x y] = ginput(2)

        subplot(2,1,2)
        plot(Analyze.Kluver(tempstart:tempstart+10000,:));
        title('SELECT KLUVER BOARD START (CLOSED)')
        [x2 y2] = ginput(1)

        startix_wrist = round(x(1))+tempstart;
        endix_wrist =  round(x(2))+tempstart;
        kluver_startix = round(x2(1)+tempstart);

        Analyze = func_SingleTargetKluver_wristpos(Analyze,startix_wrist,endix_wrist,kluver_startix,savefilename)

        save(savefilename,'Analyze')
        
    end
        
        
    

    %%
       if ~(exist([savefilename '.mat']) == 2)||(exist([savefilename1 '.mat']))
    temp = load(savefilename1)
    temp.Analyze.SR=temp.Analyze.SR_MA;
       else
               temp = load(savefilename)
       end
    temp.Analyze
    vel_thresh = 30
    SingleTargetKluverMAData = func_recalculate_kluver_epochs( temp.Analyze, vel_thresh,savefilename);
    
    allgoodix = SingleTargetKluverMAData.goodix_reach==1 & SingleTargetKluverMAData.goodix_return==1;

     for j = 1:length(SingleTargetKluverMAData.ReturnTimeix) %rare instances when return start misidentified
        if  SingleTargetKluverMAData.ReturnTimeix(j)<SingleTargetKluverMAData.start_touch_return_mouth(j,3)
            SingleTargetKluverMAData.ReturnTimeix(j)=SingleTargetKluverMAData.start_touch_return_mouth(j,3);
            j
        end        
     end
    SingleTargetKluverMAData.reaction_time = (SingleTargetKluverMAData.ReachTimeix-SingleTargetKluverMAData.TargetTime)./SingleTargetKluverMAData.SR;
    SingleTargetKluverMAData.reach_time = (SingleTargetKluverMAData.TouchTimeix-SingleTargetKluverMAData.ReachTimeix)./SingleTargetKluverMAData.SR;
    SingleTargetKluverMAData.manipulation_time = (SingleTargetKluverMAData.ReturnTimeix-SingleTargetKluverMAData.TouchTimeix)./SingleTargetKluverMAData.SR;
    SingleTargetKluverMAData.return_time = (SingleTargetKluverMAData.MouthTimeix-SingleTargetKluverMAData.ReturnTimeix)./SingleTargetKluverMAData.SR;
    
    
    TargetAppear = SingleTargetKluverMAData.TargetTime(allgoodix)./SingleTargetKluverMAData.SR;
    ReachStart = SingleTargetKluverMAData.ReachTimeix(allgoodix)./SingleTargetKluverMAData.SR;
    Touch = SingleTargetKluverMAData.TouchTimeix(allgoodix)./SingleTargetKluverMAData.SR;
    ReturnStart = SingleTargetKluverMAData.ReturnTimeix(allgoodix)./SingleTargetKluverMAData.SR;
    Mouth = SingleTargetKluverMAData.MouthTimeix(allgoodix)./SingleTargetKluverMAData.SR;

    dnex_SingleKluver = nexAddEvent(dnex_SingleKluver, TargetAppear, 'TargetAppear');
    dnex_SingleKluver = nexAddEvent(dnex_SingleKluver, ReachStart, 'ReachStart');
    dnex_SingleKluver = nexAddEvent(dnex_SingleKluver, Touch, 'Touch');
    dnex_SingleKluver = nexAddEvent(dnex_SingleKluver, ReturnStart, 'ReturnStart');
    dnex_SingleKluver = nexAddEvent(dnex_SingleKluver, Mouth, 'Mouth');

    dnex_SingleKluver = nexAddContinuous(dnex_SingleKluver, 0, SingleTargetKluverMAData.SR, SingleTargetKluverMAData.Wrist_smooth2(:,1), 'WristPosX');
    dnex_SingleKluver = nexAddContinuous(dnex_SingleKluver, 0, SingleTargetKluverMAData.SR, SingleTargetKluverMAData.Wrist_smooth2(:,2), 'WristPosY');
    dnex_SingleKluver = nexAddContinuous(dnex_SingleKluver, 0, SingleTargetKluverMAData.SR, SingleTargetKluverMAData.Wrist_smooth2(:,3), 'WristPosZ');    
    dnex_SingleKluver = nexAddContinuous(dnex_SingleKluver, 0, SingleTargetKluverMAData.SR, SingleTargetKluverMAData.Wpos_smooth2, 'Transformed2DWristPos');

%% save MA Passive Movement Times in output folder
    % Note that these passive movement times should already be corrected
    % to match with TDT timing

%     %.mat file
    savefilename1 = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_SingleTargetKluver_Analyze2' '.mat']; 
    save(savefilename1, 'SingleTargetKluverMAData')
%     save(savefilename,'passiveTimestamps','passiveOrder','passiveCondition')
    %.nex file
    savefilename = [savedir_block '\' header.MASessionName{files_ix} '_' num2str(header.MAFileName(files_ix)) header.MACleanedFileName{files_ix} '_MA_SingleTargetKluver_Analyze2']; 

    fn = [savefilename '.nex'];
    writeNexFile(dnex_SingleKluver,fn);
    fclose all;

    
    
end


