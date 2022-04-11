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

if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end


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
inputfolder = fullfile(codecorresParentfolder, 'fs500Hz', 'm4_imCohPhaseUsingFFT_EventPhase_unifiedNHP');
if strcmpi(animal, 'Kitty')
    ylimit_ftLag = [0 20];
end
if strcmpi(animal, 'Jo')
    ylimit_ftLag = [0 30];
end

histClim = [0 1];
roseRLim = [0 0.3];

image_type = 'tif';

%% save setup
savefolder = codecorresfolder;
copyfile2folder(codefilepath, fullfile(savefolder, 'code'));
savecodefolder = fullfile(savefolder, 'code');


%% Code start here
runTFLag = false;
runCicohHist = true;
runRosePlot = false;

cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);

files = dir(fullfile(inputfolder, '*.mat'));
for fi = 1 : length(files)
    
    filename = files(fi).name;
    pdcond = cond_cell{cellfun(@(x) contains(files(fi).name, x), cond_cell)};
    idx_phase = find(cellfun(@(x) contains(files(fi).name, x), EventPhases));
    if isempty(idx_phase)
        continue;
    end
    ephase = EventPhases{idx_phase};
    
    load(fullfile(inputfolder, filename), 'ciCoh', 'deltaphis_allChnsTrials', 'f_selected', 'psedociCohs', 'T_chnsarea');
    
    % select the data of chnOfI
    mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
    T_chnsarea = T_chnsarea(mask_chnOfI, :);
    ciCoh = ciCoh(mask_chnOfI, mask_chnOfI, :);
    deltaphis_allChnsTrials = deltaphis_allChnsTrials(mask_chnOfI, mask_chnOfI, :, :);
    psedociCohs = psedociCohs(mask_chnOfI, mask_chnOfI, :, :);
    
    sigciCohs= sigciCoh_extract(psedociCohs, ciCoh);
    
    % plot
    if runTFLag
        TFLagsavefolder = fullfile(savefolder, 'TFLag', ephase);
        if ~exist(TFLagsavefolder, 'dir')
            mkdir(TFLagsavefolder);
        end
        freqtimeLag_Plot(T_chnsarea, sigciCohs, deltaphis_allChnsTrials, f_selected, ylimit_ftLag, animal, ephase, pdcond, TFLagsavefolder);
        
        clear TFLagsavefolder
    end
    
    if runCicohHist || runRosePlot
        [sigciCohs_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(sigciCohs, deltaphis_allChnsTrials, T_chnsarea, 'codesavefolder', savecodefolder);
        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(ephase, animal, pdcond, 'codesavefolder', savecodefolder);
        ntrials = size(deltaphis_allChnsTrials, 4);
        
        if runCicohHist
            nshuffle = size(psedociCohs, 4);
            titlename = [animal '-'  pdcond '-'  ephase '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
            
            plot_ciCohHistogram(sigciCohs_flatten, chnPairNames, f_selected, titlename, histClim, ...
                'codesavefolder', savecodefolder, 'fig_width', 500, 'fig_height', 200);
            
            saveimgname = [animal '_' ephase '_' pdcond '_align2' char(align2) '.' image_type];
            histsavefolder = fullfile(savefolder, 'ciCohHist');
            if ~exist(histsavefolder, 'dir')
                mkdir(histsavefolder)
            end
            saveas(gcf, fullfile(histsavefolder, saveimgname), image_type);
            close(gcf)
            clear nshuffle titlename saveimgname histsavefolder
        end
        
        % rose histogram of deltaphis_allChnsTrials
        if runRosePlot
            titlename_prefix = [animal '-'  pdcond '-'  ephase];
            subtitlename = [ephase '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s, align2 = ' char(align2) ', ntrials = ' num2str(ntrials)];
            savefile_prefix = [animal 'trialPhaseDiff'];
            savefile_suffix = [pdcond '_' ephase   '_align2' char(align2)];
            rosePlotsavefolder = fullfile(savefolder, 'rosePlot', ephase);
            if ~exist(rosePlotsavefolder, 'dir')
                mkdir(rosePlotsavefolder);
            end
            
            plotsave_deltaphirose(deltaphis_flatten, sigciCohs_flatten, chnPairNames, f_selected, titlename_prefix, subtitlename, rosePlotsavefolder, savefile_prefix, savefile_suffix, image_type,...
                'codesavefolder', savecodefolder, 'roseRLim', roseRLim);
            
            clear titlename_prefix subtitlename savefile_prefix savefile_suffix
        end
        
        clear sigciCohs_flatten deltaphis_flatten chnPairNames
        clear align2 t_AOI align2name ntrials
    end
    
    clear filename pdcond ephase 
    clear sigciCohs mask_chnOfI
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
histedge = [-0.05:0.0025:0.05 ];

width_outsubplot = 180;
height_outsubplot = 410;
deltaWidth = 20;
dispix_outinpos = [30 20 20 50];
width_insubplot = width_outsubplot - dispix_outinpos(1) - dispix_outinpos(3);
height_insubplot = height_outsubplot - dispix_outinpos(2) - dispix_outinpos(4);


x_leftmargin = 20;
x_rightmargin = 0;
y_topmargin = 10;
y_bottommargin = 20;

fig_width = nfs* width_outsubplot + (nfs-1)* deltaWidth + x_leftmargin + x_rightmargin;
fig_height = height_outsubplot + y_bottommargin + y_topmargin;
figpos = [250 250 fig_width fig_height];
fig = figure('Units', 'pixels', 'Position', figpos);

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
    
   
    x_out = x_leftmargin + (fi-1) * width_outsubplot;
    x = x_out + dispix_outinpos(1);
    y = y_bottommargin + dispix_outinpos(2); 
    pos = [x y width_insubplot height_insubplot];
    ax = axes(fig, 'Units', 'pixels', 'Position', pos);
    histogram(ax, deltat, histedge)
    
    ylabel([num2str(f) 'Hz'])
    ylim(ylimit)
    view(90,90)
    
    % annotation for cicoh and median
    annotation(fig,'textbox',[x+width-0.05 y+height-0.05 0.06 0.05],...
    'String',{['cicoh = ' num2str(sigcicoh(fi))]},'LineStyle','none','FitBoxToText','off');

    
    % annotation for median deltat
    med = round(median(abs(deltat)),3)*1000;
    T = round(1/f,3)*1000;
    med_4T = [med med+T  med+2*T med+3*T];

    annotation(fig,'textbox',[x+0.005 y+0.07 0.07 0.05],'String',{['median \Deltat=[' num2str(med_4T) ']ms']},'LineStyle','none','FitBoxToText','off');

    
    if fi == nfs       
        % Create textbox
        annotation(fig,'textbox',[0.05 0.88 0.9 0.08],...
            'String',[figTitle_prefix '-freq' num2str(round(freqs(1))) '-' num2str(round(freqs(end))) 'Hz'],...
            'LineStyle','none', 'FitBoxToText','off'); 
        saveas(fig, fullfile(savefolder, [saveimgname_prefix '_f' num2str(round(freqs(1))) '_' num2str(round(freqs(end))) 'Hz']), image_type);
        close(fig);
        clear f_end
    end
end