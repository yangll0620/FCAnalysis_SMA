function ma_joint = majointdata_extract(file_matrc, joint_name)
%% extract the x, y  and z coordinates for each joint, nan data were removed
%
%   Example:
%       file_matrc = fullfile('F:','yang7003@umn','NMRC_umn', 'Projects', ...
%           'NWBStandardization','workingfolders','home','data_shared','raw',...
%           'bug','expdata', 'setupchair','bug-190111', 'ma','Bug_20190111_3_cleaned.trc');
%
%       ma_joint = madata_joint(file_matrc, 'Wrist')
%
%   Input:
%       file_matrc:  the full path of cleaned ma .trc file
%       joint_name:  the name of one joint, can be found in .trc file (i.e. 'Shoulder', 'Elbow', 'Wrist', 'Kluver')
%
%
%   Output:
%       ma_joint: the time stamp and x, y, z coordinates of joint, ntimes * 4 (timestamp, x, y, z)


%read the numerical data in the ma .trc file 
numlinestart = 7; 
dataimport = importdata(file_matrc,'\t',numlinestart-1);% reading numeric data starting from line numlinestart 

% dataimport.data (ntimes * ncols), ncols = 2 + 3 * njoints 
time = dataimport.data(:,2); %   first column (frame #), second column (time stamps)
madata = dataimport.data(:,3:end); %   third -end columns (x,y,z positions for all joints)

% parse the col number for joint
text_line4 = strsplit(dataimport.textdata{4,1});
col = find(cellfun(@(x) strcmp(x, joint_name), text_line4));
if ~isempty(col)
    col = col -2;
else
    disp(['Did not find the column number of ' joint_name])
    return
end
ma_joint = madata(:,3*(col-1)+1:col*3); % ma_joint: ntimes * 3
ma_joint = cat(2, time, ma_joint); % ma_joint: ntimes * 4 (time, x, y, z)

% remove the rows containing the nan data
[idx_rows, ~] = ind2sub(size(ma_joint),find(isnan(ma_joint))); % find rows containing nan data
idx_rows = unique(idx_rows);
ma_joint(idx_rows, :) = [];
end