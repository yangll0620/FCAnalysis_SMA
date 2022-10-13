function m4_uNHP_ciCohChanges_compEvent_Anti(animal, varargin)
% plot cicoh changes for compared conditions
%
%   Example Usage:
%           m4_uNHP_ciCohChanges_compEvent_Anti('Jo', 'compEi_str', 1, 'ci_str', 1, 'shuffleN_psedoTest', 100)
%           m4_uNHP_ciCohChanges_compEvent_Anti('Kitty', 'compEi_str', 1, 'ci_str', 1, 'shuffleN_psedoTest', 100)
%   
%   Input:
%       animal
%
%       Name-Value:
%           compEi_str - compare events start index
%           compEi_end - compare events end index
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
addParameter(p, 'compEi_str', 1, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'compEi_end', [], @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', [], @isscalar);
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
compEi_str = p.Results.compEi_str;
compEi_end = p.Results.compEi_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;

if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end

[~, ~, pipelinefolder, ~] = exp_subfolders();

%% save setup
[~,codefilename,~] = fileparts(codefilepath); 
savecodefolder = fullfile(codefolder, 'NHPs', animal, '0_dataPrep', 'SKT', 'fs500Hz', codefilename);
[savefolder, ~] = code_corresfolder(savecodefolder, false, false);

savecodefolder = fullfile(savefolder, 'code');
copyfile2folder(codefilepath, savecodefolder);

ciCohChangesfile_prefix =[animal '-ciCohChanges'];
saveimg_prefix = [animal '-ciCohChanges'];
imgtitle_prefix = [animal '-FCChanges'];

%%  input setup
ciCohfolder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'SKT', 'fs500Hz', 'm3_uNHP_ciCOH_Anti');
inputfile_cicoh_prefix = [animal '-ciCoh'];

image_type = 'tif';



%% Code start here
cond_cell = cond_cell_extract(animal);
if isempty(ci_end)
    ci_end = length(cond_cell);
end

baseEvents = {'Anticipated'};
compEvents = {'earlyReach'};

if isempty(compEi_end)
    compEi_end = length(compEvents);
end

for baseEi = 1 : length(baseEvents)
    baseEvent = baseEvents{baseEi};
    for compEi = compEi_str: compEi_end
        compEvent = compEvents{compEi};
        
        if strcmpi(baseEvent, compEvent)
            clear compEvent
            continue;
        end
        
        for ci = ci_str : ci_end
            pdcond = cond_cell{ci};
            
            ciCohChangesfile = fullfile(savefolder, [ciCohChangesfile_prefix '_' pdcond '_' compEvent '2' baseEvent '.mat']);
            
            disp([animal '-' compEvent '-based ' baseEvent '-' pdcond])
            
            %%% ----------- case no ciCohPhasefile or new Run -------- %%%
            if(~exist(ciCohChangesfile, 'file') || newRun)
                
                ciCohfile_base = fullfile(ciCohfolder, [inputfile_cicoh_prefix '_' pdcond '_' baseEvent '.mat']);
                load(ciCohfile_base, 'ciCoh',  'lfptrials', 'f_selected', 'T_chnsarea', 'fs', 'f_AOI');
                ciCoh_base = ciCoh;
                lfptrials_base = lfptrials;
                clear ciCohfile_base ciCoh lfptrials
                
                ciCohfile_comp = fullfile(ciCohfolder, [inputfile_cicoh_prefix '_' pdcond '_' compEvent '.mat']);
                load(ciCohfile_comp, 'ciCoh', 'lfptrials');
                ciCoh_comp = ciCoh;
                lfptrials_comp = lfptrials;
                clear ciCohfile_comp ciCoh lfptrials
                
                save(ciCohChangesfile, 'ciCoh_base', 'ciCoh_comp', 'lfptrials_base', 'lfptrials_comp', 'f_selected', 'T_chnsarea', 'fs', 'f_AOI')
                
                ciCohChanges = abs(ciCoh_comp) - abs(ciCoh_base);
                save(ciCohChangesfile, 'ciCohChanges', '-append');
                
                clear ciCohChanges
                clear('ciCoh_base', 'ciCoh_comp', 'lfptrials_base', 'lfptrials_comp', 'f_selected', 'T_chnsarea', 'fs', 'f_AOI')
            end
            
            %%%----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
            psedoVar = whos('-file',ciCohChangesfile, 'psedociCohChanges');
            if isempty(psedoVar) || psedoVar.size(4)< shuffleN_psedoTest
                load(ciCohChangesfile, 'lfptrials_base', 'lfptrials_comp', 'f_AOI', 'fs');
                
                psedociCohChanges_extract_save(shuffleN_psedoTest, lfptrials_comp, lfptrials_base, fs, f_AOI, ciCohChangesfile, 'codesavefolder', savecodefolder);
                
                clear('lfptrials_base', 'lfptrials_comp', 'f_AOI', 'fs');
            end
            clear psedoVar
            
            
            %%% -- plot and save section --- %%%
            load(ciCohChangesfile, 'ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea')
            
            % sig and flatten
            [sigciCohChanges]= sigciCoh_extract(psedociCohChanges, ciCohChanges);
            [ciCohChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigciCohChanges, T_chnsarea);
            
            % plot and save ciCoh Histogram image
            titlename = [imgtitle_prefix ':' pdcond '-' compEvent '/' baseEvent];
            plot_ciCohHistogram(ciCohChanges_flatten, chnPairNames, f_selected, titlename, 'histClim', [-1 1],...
                'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1]);
            saveimgname = [saveimg_prefix '-' pdcond '-' compEvent  '2' baseEvent '.' image_type];
            saveas(gcf, fullfile(savefolder, saveimgname), image_type);
            close(gcf)
            
            %%% clear
            clear pdcond ciCohChangesfile
            clear('ciCohChanges', 'psedociCohChanges', 'f_selected',  'T_chnsarea')
            clear sigciCohChanges ciCohChanges_flatten chnPairNames
            clear titlename saveimgname
        end
        
        clear ci compEvent 
    end
    clear baseEvent   compEi
end

function psedociCohChanges_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ciCohChangesfile, varargin)
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


load(ciCohChangesfile, 'psedociCohChanges');

if(~exist('psedociCohChanges', 'var'))
    psedociCohChanges = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedociCohChanges, 4) + 1;
end

lfp_combined = cat(3, lfptrials_base, lfptrials);
ntotal = size(lfp_combined, 3);
ntrials = size(lfptrials, 3);
for si = shuffi_str : suffi_end
    randomSKTInds =  randsample(ntotal,ntrials);
    
    masksBase = logical([1: ntotal]);
    masksBase(randomSKTInds) = 0;
   
    psedolfp_base = lfp_combined(:, :, masksBase);
    psedolfp_comp = lfp_combined(:, :, ~masksBase);

    
    [psedociCoh_base, ~] = ciCohSKTAllchns_FFT_NoAmp(psedolfp_base, fs, f_AOI, 'codesavefolder', codesavefolder);
    [psedociCoh_comp, ~] = ciCohSKTAllchns_FFT_NoAmp(psedolfp_comp, fs, f_AOI, 'codesavefolder', codesavefolder);
    
    psedociCohChanges = cat(4, psedociCohChanges, psedociCoh_comp - psedociCoh_base);
    
    if(mod(si, 10) == 0)
        disp(['pesdo test finished ' num2str(si)])
        save(ciCohChangesfile, 'psedociCohChanges', '-append');
    end
    
    clear randomSKTInds masksBase psedolfp_base psedolfp_comp
    clear psedociCoh_base psedociCoh_comp
end


