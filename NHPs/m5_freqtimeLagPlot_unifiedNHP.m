function m5_freqtimeLagPlot_unifiedNHP(animal, varargin)

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


%% save setup
savefolder = codecorresfolder;
copyfile2folder(codefilepath, fullfile(savefolder, 'code'));

%% Code start here
cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();

files = dir(fullfile(inputfolder, '*.mat'));
for fi = 1 : length(files)
    
    filename = files(fi).name;
    pdcond = cond_cell{cellfun(@(x) contains(files(fi).name, x), cond_cell)};
    ephase = EventPhases{cellfun(@(x) contains(files(fi).name, x), EventPhases)};
    
    load(fullfile(inputfolder, filename), 'ciCoh', 'deltaphis_allChnsTrials', 'f_selected', 'psedociCohs', 'T_chnsarea');
    
    sigciCoh= sigciCoh_extract(psedociCohs, ciCoh);
    
    for bi = 1 : length(T_chnsarea.brainarea)-1
        site1 = T_chnsarea.brainarea{bi};
        for bj = bi+1 : length(T_chnsarea.brainarea)
            site2 = T_chnsarea.brainarea{bj};
            if (contains(site1, 'stn')&& contains(site2, 'stn')) || (contains(site1, 'gp')&& contains(site2, 'gp'))
                continue;
            end
            
            contInds = contiFreqrange(sigciCoh(bi, bj, :)); % get continuous sig inds
            deltaphis_trials = squeeze(deltaphis_allChnsTrials(bi, bj, :, :));
            
            for ci = 1 : size(contInds, 1)
                idx_str = contInds(ci, 1);
                idx_end = contInds(ci, 2);
                
                figTitle_prefix = [animal '-' pdcond '-' ephase ' ' site1 '-' site2 ' time-frequency lag Histogram'];
                saveimgname_prefix = [animal '_' pdcond '_' ephase '_freqtimeLagPlot_' site1 '_' site2];
                plotTimefreqLag(f_selected(idx_str:idx_end), deltaphis_trials, figTitle_prefix, savefolder, saveimgname_prefix)
                clear idx_str idx_end
            end
        end
    end
    
    clear('filename', 'ciCoh', 'deltaphis_allChnsTrials', 'f_selected', 'psedociCohs', 'T_chnsarea');
    
    
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


function plotTimefreqLag(freqs, deltaphis_trials, figTitle_prefix, savefolder, saveimgname_prefix)
%   plot Time frequency lag for length(freqs) > 2
%
%   Input:
%       freqs: vector nf * 1
%       deltaphis_trials: delta phis across trials nf * ntrials

nfs = length(freqs);
if nfs <=2
    return
end

image_type = 'tif';

% figure setup
x0 = 0.03;
y0 = 0.1;
width = 0.07;
height = 0.8;
deltax = 0.05;
figpos = [250 250 1600 420];
histedge = [-0.05:0.005:0.05 ];
nPerFig = 8;


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
    
    if mod(fi, nPerFig) == 1
        fig = figure('Position', figpos);
        f_str = freqs(fi);
    end
    ax = axes(fig, 'Position', [x0+(mod(fi, nPerFig)-1)*(deltax + width) y0  width height]);
    histogram(ax, deltat, histedge)
    
    ylabel([num2str(f) 'Hz'])
    ylim([0 40])
    if mod(fi, nPerFig) == 1
        xlabel('\Deltat')
    end
    view(90,90)
    
    if fi == nfs || mod(fi, nPerFig) == 0
        f_end = freqs(fi);
        
        % Create textbox
        annotation(fig,'textbox',[0.34 0.92 0.3 0.07],...
            'String',[figTitle_prefix '-freq' num2str(round(f_str)) '-' num2str(round(f_end)) 'Hz'],...
            'LineStyle','none', 'FitBoxToText','off');
        
        
        saveas(fig, fullfile(savefolder, [saveimgname_prefix '_f' num2str(round(f_str)) '_' num2str(round(f_end)) 'Hz']), image_type);
        close(fig);
        clear f_end
    end
end

