
function [ TTX, FileInfoTD ] = func_AccessTDTDataTanks( filename )
% Access data tank using ActiveX functions.  This function gives you ACCESS
% to the data tanks, it does not load all the data into RAM

% Change log:
% LAJ: changed input to just be the filename rather than a header

% Note: you will need to have TDT sofware installed
% Go to: http://www.tdt.com/downloads.html
% Install:  -TDT Drivers/RPvdsEx 
%           -OpenEx  (pswd: dragonfly)
%           -OpenDeveloper (pswd: killerbee)

% INPUT:
% filename ; % includes path, e.g.
% 'Y:\Luke\Jo\Recording\Raw\rawTDT\Jo_CR1_DT1_020216\Block-2'


% OUTPUT:
% TTX is the activeX variable that you will use to access data saved in the
% data tank/block
% FileInfoTD is a structure that includes information about the data saved
% in the data tank/block
% e.g:
%           EventName: {'Tick'  'Trig'  'StPd'  'dLFP'  'Stim'  'Blck'  'StmP'  'StmO'  'Luke'  'DBSV'  'plot'}
%     SampleFrequency: [0 0 2.4414e+04 2.4414e+04 2.4414e+04 2.4414e+04 0 0 2.4414e+04 2.4414e+04 762.9395]
%       TankEventType: [257 257 33025 33025 33025 33281 257 257 33025 33025 33025]
%             NumChan: [0 0 1 16 1 1 0 0 100 1 96]
%          EventIndex: [1 2 3 4 5 6 7 8 9 10 11]
%% Define Sesssion, Block etc 

% filename = [header.server_dir header.rawTDTFolder{files_ix} header.TDTSessionName{files_ix} '\Block-' num2str(header.TDTBlockName(files_ix))]

ix1 = strfind(filename,'\')
MyTank = filename(1:(ix1(end)))
MyBlock = ['~' filename((ix1(end)+1):end)]

%% Access the TDT DataTanks
% Get access to the data tank:  
% First instantiate a variable for the ActiveX wrapper interface
h = figure
TTX = actxcontrol('TTank.X', 'Parent', h);
disp('Connect to Server');
TTX.ConnectServer('Local','Me')
disp('Open Tank');
TTX.OpenTank(MyTank,'R')
disp('Access the block (1=accessed, 0=access failed)');
TTX.SelectBlock(MyBlock) 
% Set parameters 
TTX.SetGlobalV('Channel',0);
TTX.SetGlobalStringV('Options','FILTERED');

%  for some reason often you have to run this code twice in a row for it to actually access the tank properly
close(h)
h = figure
TTX = actxcontrol('TTank.X', 'Parent', h);
disp('Connect to Server');
TTX.ConnectServer('Local','Me')
disp('Open Tank');
TTX.OpenTank(MyTank,'R')
disp('Access the block (1=accessed, 0=access failed)');
TTX.SelectBlock(MyBlock) 
TTX.SetGlobalV('Channel',0);
TTX.SetGlobalStringV('Options','FILTERED');
if ans ==0
    error('check block is accessible, in path, and try again')
end

close(h)
h = figure
TTX = actxcontrol('TTank.X', 'Parent', h);
disp('Connect to Server');
TTX.ConnectServer('Local','Me')
disp('Open Tank');
TTX.OpenTank(MyTank,'R')
disp('Access the block (1=accessed, 0=access failed)');
TTX.SelectBlock(MyBlock) 
TTX.SetGlobalV('Channel',0);
TTX.SetGlobalStringV('Options','FILTERED');
if ans ==0
    error('check block is accessible, in path, and try again')
end

% Get file header information
% get str pointer for each data storeName
block_notes=TTX.CurBlockNotes;
ix_StoreName = findstr('NAME=StoreName',block_notes);
if isempty(ix_StoreName)
    error('check block access')
end
% for each hdr item get related info
for i=1:length(ix_StoreName)
    if i==length(ix_StoreName)
        str_event=block_notes(ix_StoreName(i):end);
    else    
        str_event=block_notes(ix_StoreName(i):ix_StoreName(i+1));
    end
    ix_Value = findstr('VALUE=',str_event);
    FileInfoTD.EventName{i}=strtok(str_event((ix_Value(1)+6):end),';');
    ix_SF = findstr('SampleFreq',str_event);
    ix_Value = findstr('VALUE=',str_event(ix_SF:end));
    FileInfoTD.SampleFrequency(i)=str2num(strtok(str_event((ix_SF+ix_Value(1)-1+6):end),';'));
    ix_TET = findstr('TankEvType',str_event);
    ix_Value = findstr('VALUE=',str_event(ix_TET:end));
    FileInfoTD.TankEventType(i)=str2num(strtok(str_event((ix_TET+ix_Value(1)-1+6):end),';'));
    ix_CN = findstr('NumChan',str_event);
    ix_Value = findstr('VALUE=',str_event(ix_CN:end));
    FileInfoTD.NumChan(i)=str2num(strtok(str_event((ix_CN+ix_Value(1)-1+6):end),';'));  
    FileInfoTD.EventIndex(i)=i;
end

    TTX.SetGlobalV('WavesMemLimit', 1024^3); %if long file can run into memory issues

end

