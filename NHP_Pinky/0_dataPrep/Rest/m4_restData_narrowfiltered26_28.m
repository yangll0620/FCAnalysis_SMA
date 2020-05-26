function m4_restData_narrowfiltered26_28()
%% narrow filtered rest data recorded with grayMatter in frequency [26 28]Hz
%   for stn/gp channel first bipolar then filtered
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


%%  input setup
% band pass frequency
frebp = [26 28];
% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm2_restData_selectSeg_M1Power');

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['narrowFiltered' num2str(frebp(1)) '_' num2str(frebp(2)) '-'];

%% starting: narrow filter the lfp data of all the files
files = dir(fullfile(inputfolder, '*.mat'));
nfiles = length(files);
f = waitbar(0, ['Narrow Filtering....']);
for filei = 1 : nfiles
    % wait bar
    waitbar(filei/nfiles,f,['Narrow Filtering [' num2str(frebp(1)) ' ' num2str(frebp(2)) ']Hz lfp data in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data
    filename = files(filei).name;
    load(fullfile(inputfolder, filename), 'data_segments', 'fs', 'segsRemain', 'GMChnAreas','chans_m1','segsIndex');
    
    % band pass filter
    for segi = 1 : size(data_segments,2)
        seg_lfpm1 = data_segments(segi).lfp_m1;
        seg_lfpGM = data_segments(segi).lfp_GM;
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
        % lfp_GM
        filteddata_segments(end).lfp_GM = bandpass(seg_lfpGM, frebp, fs);
        
        
        clear seg_lfpm1 seg_lfpGM seg_lfpstn seg_lfpgp
    end
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'filteddata_segments','fs', 'segsRemain', 'GMChnAreas','chans_m1','segsIndex');
    
    
    clear data_segments fs segsRemain GMChnAreas chans_m1 segsIndex
    clear segi idx tmpn savefilename filteddata_segments
    
    
end
close(f);
disp(['narrow filtered lfpdata using grayMatter are saved to ' savefolder])


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



