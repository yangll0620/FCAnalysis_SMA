function fig5_imCohChanges_Reachfreeze2ReachPhases(varargin)
%   
%   Usage:
%       fig5_imCohChanges_Reachfreeze2ReachPhases('pos_ifig', [50 50 400 200])
%
%   Inputs:
%
%       Name-Value:
%           'pos_ifig' - position and size of the reachtime statistical figure [left bottom fig_width fig_height], default [50 50 400 200]
%  


% parse params
p = inputParser;
addParameter(p, 'pos_ifig', [50 50 400 200], @(x) assert(isvector(x) && isnumeric(x) && length(x)==4));


parse(p,varargin{:});
pos_ifig = p.Results.pos_ifig;



codefilepath = mfilename('fullpath');


% find the codefolder
tmp = regexp(codefilepath, '.*\code', 'match');
if length(tmp) ~= 1
    disp('can not find code path correctly.')
    return;
end
codefolder = tmp{1};
clear tmp

% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));

%% Input & save
[~, ~, pipelinefolder, outputfolder] = exp_subfolders();
[~, funcname, ~]= fileparts(codefilepath);

input_folder = fullfile(pipelinefolder, 'NHPs', 'Kitty', '0_dataPrep', 'SKT', 'fs500Hz', 'm5_imCohChanges_Freeze2ReachPhases');
ciCoh_Changes_file = fullfile(input_folder, 'ciCohs-Changes-Freeze2reachPhases.mat');

savefolder = fullfile(outputfolder, 'results', 'figures', funcname);
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

aisavefolder = fullfile(outputfolder,'results','figures', 'Illustrator', funcname);
if ~exist(aisavefolder, 'dir')
    mkdir(aisavefolder)
end
copy2folder = aisavefolder;

savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, "dir")
    rmdir(savecodefolder, 's');
end
copyfile2folder(codefilepath, savecodefolder);


savefilename = funcname;

%% Code start here
disp(['running ' funcname]);
load(ciCoh_Changes_file, 'ciCohChanges','psedociCohChanges', 'f_selected',  'T_chnsarea')

eBasePhases = {'preMove'; 'earlyReach'};
for ebi = 1 : length(eBasePhases)
    eBasePhase = eBasePhases{ebi};

    show_xticklabels = false;
    show_xlabel = false;
    if ebi == length(eBasePhases)
        show_xticklabels = true;
        show_xlabel = true;
    end


    reachfreezeTypes = fieldnames(ciCohChanges.(eBasePhase));
    for fri = 1 : length(reachfreezeTypes)
        subfreezeType = reachfreezeTypes{fri};
    end
end
