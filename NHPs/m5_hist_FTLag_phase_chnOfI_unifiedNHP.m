function m5_hist_FTLag_phase_chnOfI_unifiedNHP(animal, varargin)
% plot cicoh Histogram, frequency time lag and phase of interested channels
%   Input:
%       Name-Value: 
%           animal
%           ei_str - event start index
%           ci_str - condition start index

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


% find animal corresponding folder
[~, codefilename]= fileparts(codefilepath);
SKTSubfolder = 'SKT';
if strcmpi(animal, 'Kitty')
    SKTSubfolder = 'SKT_SegV';
end
NHPCodefilepath = fullfile(codefolder, 'NHPs', animal, '0_dataPrep' , SKTSubfolder, codefilename);
[codecorresfolder, codecorresParentfolder] = code_corresfolder(NHPCodefilepath, true, false);

%% Input setup
inputfolder = fullfile(codecorresParentfolder, 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');
if strcmpi(animal, 'Kitty')
    ylimit_ftLag = [0 20];
end
if strcmpi(animal, 'Jo')
    ylimit_ftLag = [0 30];
end

%% save setup
savefolder = codecorresfolder;
copyfile2folder(codefilepath, fullfile(savefolder, 'code'));


%% Code start here
cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();
chnsOfI = chnOfInterest_extract(animal, 'codesavefolder', fullfile(savefolder, 'code'));

files = dir(fullfile(inputfolder, '*.mat'));
for fi = 1 : length(files)
    
    filename = files(fi).name;
    pdcond = cond_cell{cellfun(@(x) contains(files(fi).name, x), cond_cell)};
    ephase = EventPhases{cellfun(@(x) contains(files(fi).name, x), EventPhases)};
    
    load(fullfile(inputfolder, filename), 'ciCoh', 'deltaphis_allChnsTrials', 'f_selected', 'psedociCohs', 'T_chnsarea');
    
    % select the data of chnOfI
    mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_chnOfI, :);
    ciCoh = ciCoh(mask_chnOfI, mask_chnOfI, :);
    deltaphis_allChnsTrials = deltaphis_allChnsTrials(mask_chnOfI, mask_chnOfI, :, :);
    psedociCohs = psedociCohs(mask_chnOfI, mask_chnOfI, :, :);
    
    sigciCohs= sigciCoh_extract(psedociCohs, ciCoh);
    
    % plot
    freqtimeLag_Plot(T_chnsarea, sigciCohs, deltaphis_allChnsTrials, f_selected, ylimit_ftLag, animal, ephase, pdcond, savefolder)
    
    clear filename pdcond ephase 
    clear sigciCohs
    clear('ciCoh', 'deltaphis_allChnsTrials', 'f_selected', 'psedociCohs', 'T_chnsarea');
end

function freqtimeLag_Plot(T_chnsarea, sigciCohs, deltaphis_allChnsTrials, f_selected, ylimit, animal, ephase, pdcond, savefolder)
for bi = 1 : length(T_chnsarea.brainarea)-1
    site1 = T_chnsarea.brainarea{bi};
    for bj = bi+1 : length(T_chnsarea.brainarea)
        site2 = T_chnsarea.brainarea{bj};
        if (contains(site1, 'stn')&& contains(site2, 'stn')) || (contains(site1, 'gp')&& contains(site2, 'gp'))
            clear site2
            continue;
        end
        
        contInds = contiFreqrange(sigciCohs(bi, bj, :)); % get continuous sig inds
        deltaphis_trials = squeeze(deltaphis_allChnsTrials(bi, bj, :, :));
        sigcicoh_1pair = squeeze(sigciCohs(bi, bj, :));
        
        for ci = 1 : size(contInds, 1)
            idx_str = contInds(ci, 1);
            idx_end = contInds(ci, 2);
            
            figTitle_prefix = [animal '-' ephase ' ' site1 '-' site2 '-' pdcond ' time-frequency lag Histogram'];
            saveimgname_prefix = [animal  '_' ephase '_freqtimeLagPlot_' site1 '_' site2 '-' pdcond];
            ContFreq_TFLag_Plot(f_selected(idx_str:idx_end), deltaphis_trials(idx_str:idx_end, :), sigcicoh_1pair(idx_str:idx_end), ...
                figTitle_prefix, savefolder, saveimgname_prefix, ylimit)
            
            clear idx_str idx_end figTitle_prefix saveimgname_prefix
        end
        
        clear site2
        clear contInds deltaphis_trials sigcicoh_1pair
    end
    clear site1
end

function contInds = contiFreqrange(data)
% data : vector n * 1
conTag = false;
contInds = [];
for di = 1 : length(data)
    if(data(di) > 0 && conTag == false)
        idx_str = di;
        conTag = true;
    end
    if(data(di) == 0 && conTag == true)
        contInds = cat(1, contInds, [idx_str di-1]);
        conTag = false;
    end
    if(di == length(data) && data(di) >0 && conTag == true)
        contInds = cat(1, contInds, [idx_str di]);
    end
end


function ContFreq_TFLag_Plot(freqs, deltaphis_trials, sigcicoh, figTitle_prefix, savefolder, saveimgname_prefix, ylimit)
%   plot Time frequency lag for length(freqs) > 2
%
%   Input:
%       freqs: vector nf * 1
%       deltaphis_trials: delta phis across trials nf * ntrials
%       sigcicoh: vector nf * 1

nfs = length(freqs);

image_type = 'tif';

% figure setup
leftmargin = 0.03;
rightmargin = 0.03;
uppermargin = 0.1;
lowermargin = 0.1;

deltax = 0.05;
deltay = 0.05;
figpos_base = [250 250 1600 420];
histedge = [-0.05:0.0025:0.05 ];
ncols = 8;


nrows = ceil(nfs / ncols);
figpos = figpos_base;
figpos(4) = figpos(4) * nrows;
width = ((1-leftmargin - rightmargin) - (ncols -1) * deltax)/ncols;
height = ((1- uppermargin - lowermargin)- (nrows -1)* deltay)/nrows;
fig = figure('Position', figpos);
for fi = 1 : nfs
    f = freqs(fi);
    deltaphis = deltaphis_trials(fi, :);
    
    
    % convert to [-pi pi]
    for di = 1 : length(deltaphis)
        if(deltaphis(di) <= -pi)
            deltaphis(di) = deltaphis(di) + 2 * pi;
        else
            if (deltaphis(di) > pi)
                deltaphis(di) = deltaphis(di) - 2 * pi;
            end
        end
    end
    
    deltat = deltaphis / (2 * pi * f);
    
    coli = mod(fi, ncols);
    rowi = ceil(fi/ncols);
    if coli == 0
        coli = ncols;
    end

   
    x = leftmargin + (coli -1) * (deltax + width);
    y = 1-uppermargin - height - (rowi -1) * (deltay + height);
    ax = axes(fig, 'Position', [x y  width height]);
    histogram(ax, deltat, histedge)
    
    ylabel([num2str(f) 'Hz'])
    ylim(ylimit)
    if coli == 1
        xlabel('\Deltat /s')
    end
    view(90,90)
    annotation(fig,'textbox',...
    [x+width-0.05 y+height-0.05 0.06 0.05],...
    'String',{['cicoh = ' num2str(sigcicoh(fi))]},'LineStyle','none','FitBoxToText','off');
    
    if fi == nfs       
        % Create textbox
        annotation(fig,'textbox',[0.04 0.91 0.5 0.08],...
            'String',[figTitle_prefix '-freq' num2str(round(freqs(1))) '-' num2str(round(freqs(end))) 'Hz'],...
            'LineStyle','none', 'FitBoxToText','off'); 
        saveas(fig, fullfile(savefolder, [saveimgname_prefix '_f' num2str(round(freqs(1))) '_' num2str(round(freqs(end))) 'Hz']), image_type);
        close(fig);
        clear f_end
    end
end