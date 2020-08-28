function m4_restData_filtered16_18_segDownsample()
%% narrow filtered rest data recorded with grayMatter in frequency [26 28]Hz
%
%   1. narrowed filltered. for stn/gp channel first bipolar then filtered
%
%   2. seg into intervals with same time length (remove bad segments)  
%
%   3. downsample
%
%
%   Input:
%       
%       \m2_restData_selectSeg_M1Power
%   
%       'data_segments', 'fs', 'segsRemain','segsIndex'
%       
%   Output variables:
%   
%
%        lfpsegs_m1, lfpsegs_stn lfpsegs_gp : lfp segments in m1, stn, and gp
%
%          fs: resampled samping rate (default 500Hz)
%       
%          segsRemain, segsIndex: forgot, need to be added
%   

%% folders generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx

% add util path
addpath(genpath(fullfile(codefolder,'util')));

% the corresponding pipeline and the parent folder for this code
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);


%% global variables
% animal
tmp = char(regexp(codefilepath, '/NHP_\w*/', 'match'));
animal = tmp(length('/NHP_')+1:end-1);

segt = 2;

fs_new = 500;


%%  input setup
% band pass frequency
frebp = [16 18];
% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_M1Power');

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['narrowFiltered' num2str(frebp(1)) '_' num2str(frebp(2)) '-segDownsampled-'];

%% starting: narrow filter the lfp data of all the files
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['Narrow Filtering....']);
for filei = 1 : nfiles
    % wait bar
    waitbar(filei/nfiles,f,['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data
    filename = files(filei).name;
    load(fullfile(inputfolder, filename), 'data_segments', 'fs', 'segsRemain','segsIndex');
    
    % band pass filter
    for segi = 1 : size(data_segments,2)
        seg_lfpm1 = data_segments(segi).lfp_m1;
        seg_lfpstn = data_segments(segi).lfp_stn;
        seg_lfpgp = data_segments(segi).lfp_gp;
        
        
        % bipolar stn and gp
        seg_lfpstn = diff(seg_lfpstn, [],2);
        seg_lfpgp = diff(seg_lfpgp, [],2);
        
        % band pass filter
        if ~exist('filteddata_segments')
            filteddata_segments(1).lfp_m1 = bandpass(seg_lfpm1, frebp, fs);
        else
            filteddata_segments(end+1).lfp_m1 = bandpass(seg_lfpm1, frebp, fs);
        end
        % lfp_stn and lfp_gp
        filteddata_segments(end).lfp_stn = bandpass(seg_lfpstn, frebp, fs);
        filteddata_segments(end).lfp_gp = bandpass(seg_lfpgp, frebp, fs);
        
        
        clear seg_lfpm1 seg_lfpGM seg_lfpstn seg_lfpgp
    end
    
    
    % seg
    [lfpsegs_m1, lfpsegs_stn, lfpsegs_gp] = segPerFile(filteddata_segments, segsRemain,segt, fs);
    
    if(isempty(lfpsegs_m1)) % no good segments remained
        continue;
    end
    
    
    % down sample
    lfpsegs_m1 = resampleSeg(lfpsegs_m1, fs, fs_new);
    lfpsegs_stn = resampleSeg(lfpsegs_stn, fs, fs_new);
    lfpsegs_gp = resampleSeg(lfpsegs_gp, fs, fs_new);
    fs = fs_new;
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ...
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    
    save(fullfile(savefolder, savefilename), 'lfpsegs_m1','lfpsegs_stn','lfpsegs_gp', ...
        'fs', 'segsRemain','segsIndex');
    
    
    clear data_segments fs segsRemain  segsIndex
    clear segi idx tmpn savefilename filteddata_segments
    
    
end
close(f);
disp(['narrow filtered lfpdata using grayMatter are saved to ' savefolder])
end


function filteredX = bandpass(X, frebp, fs)
%% band pass for each channel X: ntemporal * nchns

[ntemp, nchns] = size(X);

if ntemp == 1
    X = X';
    [ntemp, nchns] = size(X);
end

filteredX = zeros(ntemp, nchns);
for chi = 1: nchns
    filteredX(:, chi) = filter_bpbutter(X(:,chi),frebp, fs);
end
end



function [lfpsegs_m1, lfpsegs_stn, lfpsegs_gp]  = segPerFile(filteddata_segments, segsRemain, segt, fs)

lfpsegs_m1 =[];
lfpsegs_stn =[];
lfpsegs_gp =[];
for segi = 1: size(filteddata_segments,2)
    
    if segsRemain(segi) ==0 % bad segment, not used
        continue;
    end
    
    seg_lfpm1 = filteddata_segments(segi).lfp_m1;
    seg_lfpstn = filteddata_segments(segi).lfp_stn;
    seg_lfpgp = filteddata_segments(segi).lfp_gp;
    
    
    % seg
    ntemp = size(seg_lfpm1,1);
    nsegtemp = ceil(fs * segt);
    tsegn = floor(ntemp/nsegtemp);
    for i = 1: tsegn
        idx_str = (i-1) * nsegtemp + 1;
        idx_end = i * nsegtemp;
        
        lfpsegs_m1 = cat(3, lfpsegs_m1, seg_lfpm1(idx_str:idx_end,:));
        lfpsegs_stn = cat(3, lfpsegs_stn, seg_lfpstn(idx_str:idx_end,:));
        lfpsegs_gp = cat(3, lfpsegs_gp, seg_lfpgp(idx_str:idx_end,:));
        
        clear idx_str idx_end
    end
    
    clear seg_lfpm1 seg_lfpstn seg_lfpgp
    clear ntemp nsegtemp tsegn
end
end

function resampledX = resampleSeg(X, fs_old, fs_new)
%% downsample along the first dim X: ntemp * nchns * ntrials
%

resampledX = [];
for triali = 1: size(X, 3)
    resamp = resample(X(:,:,triali), round(fs_new), round(fs_old));
    
    resampledX = cat(3, resampledX, resamp);
end

end
