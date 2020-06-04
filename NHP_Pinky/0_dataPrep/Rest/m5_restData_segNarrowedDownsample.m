function m5_restData_segNarrowedDownsample()
%% seg into intervals with same time length (remove bad segments) and downsample
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
segt = 2;

fs_new = 500;

% input folder: extracted raw rest data with grayMatter 
inputfolder = fullfile(codecorresParentfolder, 'm4_restData_narrowfiltered26_28');

%% save setup
savefolder = codecorresfolder;
savefilename_addstr = ['seged'];



%% Start Here

files = dir(fullfile(inputfolder, '*moderate*.mat'));
nfiles = length(files);
f = waitbar(0, ['Seging Filtered Rest Data....']);


for filei = 1 : nfiles
    disp(files(filei).name)
    
    % wait bar
    waitbar(filei/nfiles,f,['Seging and downsample Filtered Rest Data in file ' num2str(filei) '/' num2str(nfiles)]);
    
    % load data
    filename = files(filei).name;
    load(fullfile(inputfolder, filename), 'filteddata_segments', 'fs', 'segsRemain', 'GMChnAreas','chans_m1','segsIndex');
    
    
    [lfpsegs_m1, lfpsegs_GM, lfpsegs_stn, lfpsegs_gp] = segPerFile(filteddata_segments, segsRemain,segt, fs);
    
    if(isempty(lfpsegs_m1)) % no good segments remained
        continue;
    end
    
    % down sample
    lfpsegs_m1 = resampleSeg(lfpsegs_m1, fs, fs_new);
    lfpsegs_GM = resampleSeg(lfpsegs_GM, fs, fs_new);
    lfpsegs_stn = resampleSeg(lfpsegs_stn, fs, fs_new);
    lfpsegs_gp = resampleSeg(lfpsegs_gp, fs, fs_new);
    fs = fs_new;
    
    % save
    idx = strfind(filename, [animal '_']);
    tmpn = length([animal '_']);
    savefilename = [filename(idx:idx+tmpn-1) savefilename_addstr ... 
        upper(filename(idx+tmpn)) filename(idx+tmpn+1:end)];
    save(fullfile(savefolder, savefilename), 'lfpsegs_m1','lfpsegs_GM','lfpsegs_stn','lfpsegs_gp', ...
        'fs', 'segsRemain', 'GMChnAreas','chans_m1','segsIndex');
    
    
end
end


function [lfpsegs_m1, lfpsegs_GM, lfpsegs_stn, lfpsegs_gp]  = segPerFile(filteddata_segments, segsRemain, segt, fs)

lfpsegs_m1 =[];
lfpsegs_GM =[];
lfpsegs_stn =[];
lfpsegs_gp =[];
for segi = 1: size(filteddata_segments,2)
    
    if segsRemain(segi) ==0 % bad segment, not used
        continue;
    end
    
    seg_lfpm1 = filteddata_segments(segi).lfp_m1;
    seg_lfpGM = filteddata_segments(segi).lfp_GM;
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
        lfpsegs_GM = cat(3, lfpsegs_GM, seg_lfpGM(idx_str:idx_end,:));
        lfpsegs_stn = cat(3, lfpsegs_stn, seg_lfpstn(idx_str:idx_end,:));
        lfpsegs_gp = cat(3, lfpsegs_gp, seg_lfpgp(idx_str:idx_end,:));
        
        clear idx_str idx_end
    end
    
    clear seg_lfpm1 seg_lfpGM seg_lfpstn seg_lfpgp
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

