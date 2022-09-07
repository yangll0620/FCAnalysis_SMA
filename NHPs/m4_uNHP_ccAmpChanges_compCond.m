function m4_uNHP_ccAmpChanges_compCond(animal, varargin)
% plot cicoh changes for compared conditions
%
%   Example Usage:
%           m4_uNHP_ccAmpChanges_compCond('Jo', 'compCi_str', 1, 'ei_str', 3, 'ei_end', 3)
%           m4_uNHP_ccAmpChanges_compCond('Kitty', 'compCi_str', 1, 'ei_str', 1, 'shuffleN_psedoTest', 10, 'newRun', true)
%   
%   Input:
%       animal
%
%       Name-Value:
%           compCi_str - compare cond start index
%           compCi_end - compare cond end index
%           ei_str - event start index
%           ei_end - event end index
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
addParameter(p, 'compCi_str', 1, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'compCi_end', [], @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', [], @isscalar);
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

parse(p,varargin{:});
compCi_str = p.Results.compCi_str;
compCi_end = p.Results.compCi_end;
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
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
saveimg_prefix = [animal '-ccAmpChanges'];
imgtitle_prefix = [animal '-FCChanges'];

%%  input setup
ccAmpfolder = fullfile(pipelinefolder, 'NHPs', animal, '0_dataPrep' , 'SKT', 'fs500Hz', 'm3_uNHP_ccAmp');
inputfile_cicoh_prefix = [animal '-ccAmp_'];


image_type = 'tif';


%% Code start here

EventPhases = SKT_eventPhases_extract(animal, 'codesavefolder', savecodefolder);

if isempty(ei_end)
    ei_end = length(EventPhases);
end
[tbl_compConds, ~]= compCondEvents_extract(animal);
if isempty(compCi_end)
    compCi_end = height(tbl_compConds);
end

for compCi = compCi_str: compCi_end
    basepd = tbl_compConds.baseCond{compCi};
    comppd = tbl_compConds.compCond{compCi};
    
   
    for ei = ei_str : ei_end
        event = EventPhases{ei};
        [~, ~, align2name_base] = SKT_EventPhase_align2_tAOI_extract(event, animal, 'pdcond', basepd, 'codesavefolder', savecodefolder);
        [~, ~, align2name_comp] = SKT_EventPhase_align2_tAOI_extract(event, animal, 'pdcond', comppd, 'codesavefolder', savecodefolder);
        
        ccAmpChangesfile = fullfile(savefolder, [ccAmpChangesfile_prefix  '_' comppd '2' basepd '_' event '.mat']);
        
        disp([animal '-' comppd '2 ' basepd '-' event])
        
        %%% ----------- case no ciCohPhasefile or new Run -------- %%%
        if(~exist(ccAmpChangesfile, 'file') || newRun)

            ccAmpfile_base = fullfile(ccAmpfolder, [inputfile_cicoh_prefix basepd '_' event '_align2' align2name_base '.mat']);
            load(ccAmpfile_base, 'ccAmp',  'lfptrials', 'f_selected', 'T_chnsarea', 'fs', 'f_AOI');
            ccAmp_base = ccAmp; 
            lfptrials_base = lfptrials;
            clear ccAmpfile_base ccAmp lfptrials
            
            ccAmpfile_comp = fullfile(ccAmpfolder, [inputfile_cicoh_prefix comppd '_' event '_align2' align2name_comp '.mat']);
            load(ccAmpfile_comp, 'ccAmp',  'lfptrials');
            ccAmp_comp = ccAmp; 
            lfptrials_comp = lfptrials;
            clear ciCohfile_comp ccAmp lfptrials
            
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
        load(ccAmpChangesfile, 'ccAmpChanges', 'psedoccAmpChanges', 'f_selected', 'T_chnsarea')
        
        % sig and flatten
        [sigccAmpChanges]= sigciCoh_extract(psedoccAmpChanges, ccAmpChanges);
        [ccAmpChanges_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(sigccAmpChanges, T_chnsarea);
        
        % plot and save ciCoh Histogram image
        titlename = [imgtitle_prefix ':' comppd '/' basepd '-'  event];
        plot_ciCohHistogram(ccAmpChanges_flatten, chnPairNames, f_selected, titlename, 'histClim', [-1 1],...
            'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder, 'cbarStr', 'ciCohChange', 'cbarTicks', [-0.5 0 0.5]);
        saveimgname = [saveimg_prefix '-' event '-' comppd  '2' basepd '.' image_type];
        saveas(gcf, fullfile(savefolder, saveimgname), image_type);
        close(gcf)
        
        %%% clear
        clear pdcond ciCohChangesfile
        clear ccAmpChanges f_selected T_chnsarea
        clear ccAmpChanges_flatten chnPairNames
        clear titlename saveimgname nshuffle
        clear align2_comp align2name_comp align2name_base align2_base
    end
    
    clear basepd comppd 
    clear files_base files_comp
    
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
