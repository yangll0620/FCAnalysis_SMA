function m4_uNHP_ccAmpChanges_compEvent(animal, varargin)
% plot cicoh changes for compared conditions
%
%   Example Usage:
%           m4_uNHP_ccAmpChanges_compEvent('Kitty', 'compEi_str', 1, 'ci_str', 1)
%           m4_uNHP_ccAmpChanges_compEvent('Kitty', 'compEi_str', 1, 'ci_str', 1, 'ci_end', 1, 'shuffleN_psedoTest', 500, 'newRun', true)
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

ccAmpChangesfile_prefix =[animal '-ccAmpChanges'];
saveimg_prefix = [animal '-ciCohChanges'];
imgtitle_prefix = [animal '-FCChanges'];

%%  input setup
ccAmpfolder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'SKT', 'fs500Hz', 'm3_uNHP_ccAmp');
inputfile_cicoh_prefix = [animal '-ccAmp_'];

image_type = 'tif';



%% Code start here
cond_cell = cond_cell_extract(animal);
if isempty(ci_end)
    ci_end = length(cond_cell);
end
[~, tbl_compEvents]= compCondEvents_extract(animal);
if isempty(compEi_end)
    compEi_end = height(tbl_compEvents);
end

for compEi = compEi_str: compEi_end
    baseEvent = tbl_compEvents.baseEvent{compEi};
    compEvent = tbl_compEvents.compEvent{compEi};
    
    if ~strcmpi(baseEvent, 'preMove')
        clear compEvent baseEvent  
        continue;
    end

    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};
        
        [~, ~, align2name_base] = SKT_EventPhase_align2_tAOI_extract(baseEvent, animal, 'pdcond', pdcond, 'codesavefolder', savecodefolder);
        [~, ~, align2name_comp] = SKT_EventPhase_align2_tAOI_extract(compEvent, animal, 'pdcond', pdcond, 'codesavefolder', savecodefolder);
        
        ccAmpChangesfile = fullfile(savefolder, [ccAmpChangesfile_prefix '_' pdcond '_' compEvent '2' baseEvent '.mat']);

        disp([animal '-' compEvent '-based ' baseEvent '-' pdcond])

        %%% ----------- case no ciCohPhasefile or new Run -------- %%%
        if(~exist(ccAmpChangesfile, 'file') || newRun)
            
            ccAmpfile_base = fullfile(ccAmpfolder, [inputfile_cicoh_prefix pdcond '_' baseEvent '_align2' align2name_base '.mat']);
            load(ccAmpfile_base, 'ccAmp',  'lfptrials', 'f_selected', 'T_chnsarea', 'fs', 'f_AOI');
            ccAmp_base = ccAmp; 
            lfptrials_base = lfptrials;
            clear ccAmpfile_base ccAmp lfptrials
            
            ccAmpfile_comp = fullfile(ccAmpfolder, [inputfile_cicoh_prefix pdcond '_' compEvent '_align2' align2name_comp '.mat']);
            load(ccAmpfile_comp, 'ccAmp', 'lfptrials');
            ccAmp_comp = ccAmp; 
            lfptrials_comp = lfptrials;
            clear ccAmpfile_comp ccAmp lfptrials
            
            save(ccAmpChangesfile, 'ccAmp_base', 'ccAmp_comp', 'lfptrials_base', 'lfptrials_comp', 'f_selected', 'T_chnsarea', 'fs', 'f_AOI')
            
            ccAmpChanges = abs(ccAmp_comp) - abs(ccAmp_base);
            save(ccAmpChangesfile, 'ccAmpChanges', '-append');

            clear ccAmpChanges
            clear('ccAmp_base', 'ccAmp_comp', 'lfptrials_base', 'lfptrials_comp', 'f_selected', 'T_chnsarea', 'fs', 'f_AOI')
        end
        
        %%%----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
        psedoVar = whos('-file',ccAmpChangesfile, 'psedoccAmpChanges');
        if isempty(psedoVar) || psedoVar.size(4)< shuffleN_psedoTest
            load(ccAmpChangesfile, 'lfptrials_base', 'lfptrials_comp', 'f_AOI', 'fs');
            
            psedoccAmpChanges_extract_save(shuffleN_psedoTest, lfptrials_comp, lfptrials_base, fs, f_AOI, ccAmpChangesfile, 'codesavefolder', savecodefolder);
            
            clear('lfptrials_base', 'lfptrials_comp', 'f_AOI', 'fs');
        end
        clear psedoVar
            

        %%% -- plot and save section --- %%%
        load(ccAmpChangesfile, 'ccAmpChanges', 'psedoccAmpChanges', 'f_selected',  'T_chnsarea')
        
        % sig and flatten
        [sigccAmpChanges]= sigciCoh_extract(psedoccAmpChanges, ccAmpChanges);
        [ccAmpChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigccAmpChanges, T_chnsarea);
        
        % plot and save ciCoh Histogram image
        titlename = [imgtitle_prefix ':' pdcond '-' compEvent '/' baseEvent];
        plot_ciCohHistogram(ccAmpChanges_flatten, chnPairNames, f_selected, titlename, 'histClim', [-1 1],...
            'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-1 0 1]);
        saveimgname = [saveimg_prefix '-' pdcond '-' compEvent  '2' baseEvent '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        close(gcf)
        
        %%% clear
        clear pdcond ccAmpChangesfile align2name_base align2name_comp
        clear('ccAmpChanges', 'psedoccAmpChanges', 'f_selected',  'T_chnsarea')
        clear sigccAmpChanges ccAmpChanges_flatten chnPairNames
        clear titlename saveimgname
    end
    
    clear ci compEvent baseEvent     
end


function psedoccAmpChanges_extract_save(suffi_end, lfptrials, lfptrials_base, fs, f_AOI, ccAmpChangesfile, varargin)
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
% save psedoccAmpChanges: nchns * nchns * nf * nshuffle, saved to ccAmpChangesfile


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


load(ccAmpChangesfile, 'psedoccAmpChanges');

if(~exist('psedoccAmpChanges', 'var'))
    psedoccAmpChanges = [];
    shuffi_str = 1;
else
    shuffi_str = size(psedoccAmpChanges, 4) + 1;
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

    
    [psedoccAmp_base, ~] = ccAmp_SKTAllchns_FFT(psedolfp_base, fs, f_AOI, 'codesavefolder', codesavefolder);
    [psedoccAmp_comp, ~] = ccAmp_SKTAllchns_FFT(psedolfp_comp, fs, f_AOI, 'codesavefolder', codesavefolder);
    
    psedoccAmpChanges = cat(4, psedoccAmpChanges, psedoccAmp_comp - psedoccAmp_base);
    
    if(mod(si, 10) == 0)
        disp(['pesdo test finished ' num2str(si)])
        save(ccAmpChangesfile, 'psedoccAmpChanges', '-append');
    end
    
    clear randomSKTInds masksBase psedolfp_base psedolfp_comp
    clear psedoccAmp_base psedoccAmp_comp
end


