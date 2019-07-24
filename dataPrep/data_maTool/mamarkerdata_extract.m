function [ma_marker, fs] = mamarkerdata_extract(file_matrc, marker_name)
%% extract the x, y  and z coordinates for each marker, nan data were removed
%
%   Example:
%       file_matrc = fullfile('F:','yang7003@umn','NMRC_umn', 'Projects', ...
%           'NWBStandardization','workingfolders','home','data_shared','raw',...
%           'bug','expdata', 'setupchair','bug-190111', 'ma','Bug_20190111_3_cleaned.trc');
% 
%       ma_marker = mamarkerdata_extract(file_matrc, 'Wrist')
%
%   Input:
%       file_matrc:  the full path of cleaned ma .trc file
%       marker_name:  the name of one marker, can be found in .trc file (i.e. 'Shoulder', 'Elbow', 'Wrist', 'Kluver')
%
%
%   Output:
%       ma_marker: the time stamp and x, y, z coordinates of marker, ntimes * 4 (timestamp, x, y, z)
%       fs: the sample rate 

%read the numerical data in the ma .trc file 
numlinestart = 7; % start line of ma marker data 
dataimport = importdata(file_matrc,'\t',numlinestart-1);% reading numeric data starting from line numlinestart 

% dataimport.data (ntimes * ncols), ncols = 2 + 3 * nmarkers 
time = dataimport.data(:,2); %   first column (frame #), second column (time stamps)
madata = dataimport.data(:,3:end); %   third -end columns (x,y,z positions for all markers)

% parse the col number for marker
text_line4 = strsplit(dataimport.textdata{4,1});
col = find(cellfun(@(x) strcmp(x, marker_name), text_line4));
if ~isempty(col)
    col = col -2;
else
    disp(['Did not find the column number of ' marker_name])
    ma_marker = 0;
    fs = 0;
    return
end
ma_marker = madata(:,3*(col-1)+1:col*3); % ma_marker: ntimes * 3
ma_marker = cat(2, time, ma_marker); % ma_marker: ntimes * 4 (time, x, y, z)

% remove the rows containing the nan data
[idx_rows, ~] = ind2sub(size(ma_marker),find(isnan(ma_marker))); % find rows containing nan data
idx_rows = unique(idx_rows);
ma_marker(idx_rows, :) = [];


% extract the sampling rate Fs
fid = fopen(file_matrc);
linenum = 3; % line number for sample rate
C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-2);
varnames = split(C{1}{1});
coli_fs = find(cellfun(@(x) strcmp(x, 'DataRate'), varnames));
if ~isempty(coli_fs)
    %reset the file pointer to the start of the file 
    fseek(fid, 0, 'bof');  
    C = textscan(fid,'%s',1,'delimiter','\n', 'headerlines',linenum-1);
    varvalue = split(C{1}{1});
    fs = str2num(varvalue{coli_fs});
else
    disp(['Can not find the sample rate'])
end
fclose(fid);
end