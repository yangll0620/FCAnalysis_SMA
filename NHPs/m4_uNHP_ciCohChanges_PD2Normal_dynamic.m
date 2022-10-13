function m4_uNHP_ciCohChanges_PD2Normal_dynamic(animal, varargin)
% plot cicoh changes for compared conditions
%
%   Example Usage:
%           m4_uNHP_ciCohChanges_PD2Normal_dynamic('Kitty', 'newRun', false, 'shuffleN_psedoTest', 50)
%   
%   Input:
%       animal
%
%       Name-Value:
%           ci_str - condition start index
%           ci_end - condition end index
%           newRun - true or false(default), for running new or not
%           shuffleN_psedoTest -  default 500

codefilepath = mfilename('fullpath');

% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add path
addpath(genpath(fullfile(codefolder,'util')));
addpath(genpath(fullfile(codefolder,'NHPs')));
addpath(genpath(fullfile(codefolder,'connAnalyTool')));
addpath(genpath(fullfile(codefolder,'toolbox')));


% parse params
p = inputParser;
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;

if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end


[~,codefilename,~] = fileparts(codefilepath); 
if strcmpi(animal, 'Kitty')
    NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'SKT','fs500Hz', 'longerTrials', codefilename);
end
if strcmpi(animal, 'Jo')
    NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , 'SKT','fs500Hz',  codefilename);
end
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);


%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

ciCohChangesfile_prefix =[animal '-ciCohChanges'];
saveimg_prefix = [animal '-ciCohChanges'];
imgtitle_prefix = [animal '-FCChanges'];

%%  input setup
ciCohfolder = fullfile(codecorresParentfolder, 'm3_uNHP_ciCOH_dynamic');

inputfile_cicoh_prefix = [animal '-ciCoh_Dynamic_'];
image_type = 'tif';



%% Code start here
baseCond = 'normal';
compCond = 'PD';

    
ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix  '_' compCond '2' baseCond '.mat']);
disp([animal '-' compCond '-based ' baseCond])

%%% ----------- case no ciCohPhasefile or new Run -------- %%%
if(~exist(ciCohChangesfile, 'file') || newRun)
    
    ciCohfile_base = fullfile(ciCohfolder, [inputfile_cicoh_prefix  '_' baseCond '.mat']);
    load(ciCohfile_base, 'f_selected', 'T_chnsarea', 'fs', 'f_AOI', 't_selected', 't_trialdur')
    load(ciCohfile_base, 'ciCoh',  'lfptrials_win','lfptrials', 'psedociCohs');
    ciCoh_base = ciCoh;
    lfptrials_win_base = lfptrials_win;
    lfptrials_base = lfptrials;
    psedociCohs_base = psedociCohs;
    clear ciCohfile_base ciCoh lfptrials_win lfptrials psedociCohs
    
    
    ciCohfile_comp = fullfile(ciCohfolder, [inputfile_cicoh_prefix  '_' compCond '.mat']);
    load(ciCohfile_comp, 'ciCoh',  'lfptrials_win','lfptrials', 'psedociCohs');
    ciCoh_comp = ciCoh;
    lfptrials_win_comp = lfptrials_win;
    lfptrials_comp = lfptrials;
    psedociCohs_comp = psedociCohs;
    clear ciCohfile_comp ciCoh lfptrials_win lfptrials psedociCohs
    

    ciCohs.(baseCond) = ciCoh_base;
    lfptrials_wins.(baseCond) = lfptrials_win_base;
    lfptrialss.(baseCond) = lfptrials_base;
    psedociCohs.(baseCond) = psedociCohs_base;
    
    ciCohs.(compCond) = ciCoh_comp;
    lfptrials_wins.(compCond) = lfptrials_win_comp;
    lfptrialss.(compCond) = lfptrials_comp;
    psedociCohss.(compCond) = psedociCohs_comp;
    
    
    % save
    save(ciCohChangesfile, 'baseCond', 'compCond');
    save(ciCohChangesfile, 'f_selected', 'T_chnsarea', 'fs', 'f_AOI', 't_selected', 't_trialdur', '-append');
    save(ciCohChangesfile, 'ciCohs', 'lfptrials_wins', 'lfptrialss', 'psedociCohss', '-append');
    
    ciCohChanges = abs(ciCoh_comp) - abs(ciCoh_base);
    save(ciCohChangesfile, 'ciCohChanges', '-append');
    
    clear ciCoh_base lfptrials_win_base lfptrials_base psedociCohs_base
    clear ciCoh_comp lfptrials_win_comp lfptrials_comp psedociCohs_comp
    clear('f_selected', 'T_chnsarea', 'fs', 'f_AOI', 't_selected', 't_trialdur');
    clear('ciCohs', 'lfptrials_wins', 'lfptrialss', 'psedociCohss', '-append');  
end

%%%----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
load(ciCohChangesfile, 'lfptrials_wins', 'f_AOI', 'fs', 'baseCond', 'compCond');
lfptrials_win_base = lfptrials_wins.(baseCond); % nchns * ntemp * ntrials * nts
lfptrials_win_comp = lfptrials_wins.(compCond);
for ti = 1 : size(lfptrials_win_base, 4)
    lfptrials_base = squeeze(lfptrials_win_base(:, :, :, ti));
    lfptrials_comp = squeeze(lfptrials_win_comp(:, :, :, ti));
    psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials_comp, lfptrials_base, fs, f_AOI, ciCohChangesfile, ti, 'codesavefolder', savecodefolder);
    clear lfptrials_base lfptrials_comp
end
clear('lfptrials_wins', 'f_AOI', 'fs', 'baseCond', 'compCond');



%%% -- plot and save section --- %%% 

% ciCohChanges: nchns * nchns * nf * nt
load(ciCohChangesfile, 'ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea', 't_selected')

time0name = 'reachonset';
nchns = size(ciCohChanges, 1);
chnsNames = T_chnsarea.brainarea;
chnsNames{cellfun(@(x) contains(x, 'stn'), chnsNames)} = 'STN';
chnsNames{cellfun(@(x) contains(x, 'gp'), chnsNames)} = 'GP';

t_AOI = [-1 1];
mask_tAOI = (t_selected >= t_AOI(1) & t_selected <= t_AOI(2));
ciCohChanges_tAOI = ciCohChanges(:,:, :, mask_tAOI);
t_selected_AOI = t_selected(mask_tAOI);
clear mask_tAOI

% plot original ciCohChanges
title_prefix = 'original-ciCoh-Changes-';
histClim = [-1 1];
for chi = 1 : nchns -1
    chnnamei = chnsNames{chi};
    for chj = chi + 1 : nchns
        chnnamej = chnsNames{chj};
        ciCohChange_1pair = squeeze(ciCohChanges_tAOI(chi, chj, :, :));
        titlename = [title_prefix chnnamei '-' chnnamej];
        
        plot_ciCohHist_dynamic(ciCohChange_1pair, f_selected, t_selected_AOI, titlename, time0name, histClim);
        
        saveimgname = titlename;
        saveas(gcf, fullfile(savefolder, saveimgname), 'tif');
        close(gcf)
        
        clear chnnamej ciCohChange_1pair titlename saveimgname
    end
    clear chnnamei
end
clear title_prefix histClim


function psedociCohChanges_extract_save(suffi_end, lfptrials_comp, lfptrials_base, fs, f_AOI, ciCohChangesfile, ti, varargin)
%
%   Inputs:
%       suffi_end
%       lfptrials: nchns * ntemp * ntrials(nsegs)
%       fs
%       f_AOI
%       ciCohPhasefile:
%       
%       Name-Value: 
%           'codesavefolder' - code saved folder  
%   
% save psedociCohChanges: nchns * nchns * nf * nshuffle, saved to ciCohChangesfile


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


psedoVar = whos('-file',ciCohChangesfile, 'psedociCohChanges');
if(isempty(psedoVar))
    psedociCohChanges = struct();
else
    load(ciCohChangesfile, 'psedociCohChanges');
end
clear psedoVar


fiename = ['t' num2str(ti)];
if ~isfield(psedociCohChanges, fiename)
    psedociCohChanges.(fiename) = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChanges.(fiename), 4) + 1;
end

if shuffi_str > suffi_end
    return;
end

disp(['ti = ' num2str(ti) ' psedo changes test start at ' num2str(shuffi_str) ' times'])
lfp_combined = cat(3, lfptrials_base, lfptrials_comp);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials_comp, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
   
    psedolfp_base = lfp_combined(:, :, masksBase);
    psedolfp_comp = lfp_combined(:, :, ~masksBase);

    
    [psedociCoh_base, ~] = ciCohSKTAllchns_FFT_NoAmp(psedolfp_base, fs, f_AOI, 'codesavefolder', codesavefolder);
    [psedociCoh_comp, ~] = ciCohSKTAllchns_FFT_NoAmp(psedolfp_comp, fs, f_AOI, 'codesavefolder', codesavefolder);
    
    psedociCohChanges.(fiename) = cat(4, psedociCohChanges.(fiename), psedociCoh_comp - psedociCoh_base);
    
    if(mod(si, 10) == 0)
        disp(['ti = ' num2str(ti) ' pesdo changes test finished ' num2str(si)])
        save(ciCohChangesfile, 'psedociCohChanges', '-append');
    end
    
    clear randomSKTInds masksBase psedolfp_base psedolfp_comp
    clear psedociCoh_base psedociCoh_comp
end


