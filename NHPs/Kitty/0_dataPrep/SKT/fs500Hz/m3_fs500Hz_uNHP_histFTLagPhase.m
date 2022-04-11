function m3_fs500Hz_uNHP_histFTLagPhase(animal, varargin)
% plot cicoh Histogram, frequency time lag and phase of interested channels
%
%   Input:
%       animal
%
%       Name-Value:
%           ei_str - event start index
%           ei_end - event end index
%           ci_str - event start index
%           ci_end - condition end index
%           runCicohHist -  true (default) or false(default)
%           runRosePlot -  true or false(default)
%           newRun - true or false(default), for running new or not

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
addParameter(p, 'ei_str', 1, @isscalar);
addParameter(p, 'ei_end', [], @isscalar);
addParameter(p, 'ci_str', 1, @isscalar);
addParameter(p, 'ci_end', [], @isscalar);
addParameter(p, 'runCicohHist', true, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'runRosePlot', false, @(x) assert(islogical(x) && isscalar(x)));
addParameter(p, 'shuffleN_psedoTest', 500, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'newRun', false, @(x) assert(islogical(x) && isscalar(x)));

if ~exist('animal', 'var')
    animal = input('animal = ', 's');
end

parse(p,varargin{:});
ei_str = p.Results.ei_str;
ei_end = p.Results.ei_end;
ci_str = p.Results.ci_str;
ci_end = p.Results.ci_end;
runCicohHist = p.Results.runCicohHist;
runRosePlot = p.Results.runRosePlot;
shuffleN_psedoTest = p.Results.shuffleN_psedoTest;
newRun = p.Results.newRun;
disp('p.Results =  ' )
p.Results


% codecorresfolder, codecorresParentfolder
[codecorresfolder, codecorresParentfolder] = code_corresfolder(codefilepath, true, false);

%% save setup
savefolder = codecorresfolder;
savecodefolder = fullfile(savefolder, 'code');
if exist(savecodefolder, 'dir')
    rmdir(savecodefolder,'s');
end
copyfile2folder(codefilepath, savecodefolder);

phsubfolder = fullfile(savefolder, 'phases');
if ~exist(phsubfolder, 'dir')
    mkdir(phsubfolder);
end


ciCohPhasefile_prefix =[animal ' ciCohPhasefile'];

%%  input setup
inputfolder = fullfile(codecorresParentfolder, 'm2_segSKTData_SelectTrials_chnOfI');

image_type = 'tif';
f_AOI = [8 40];

if runCicohHist
    histClim = [0 1];
    
    histsavefolder = fullfile(savefolder, 'ciCohHist');
    if ~exist(histsavefolder, 'dir')
        mkdir(histsavefolder)
    end
end
if runRosePlot
    roseRLim = [0 0.3];
    
    savefile_prefix = [animal 'trialPhaseDiff'];
end


%% Code start here
cond_cell = cond_cell_extract(animal);
EventPhases = SKT_eventPhases_extract();
if isempty(ei_end)
    ei_end = length(EventPhases);
end
if isempty(ci_end)
    ci_end = length(cond_cell);
end
chnsOfI = chnsOfInterest_extract(animal, 'codesavefolder', savecodefolder);

for ei = ei_str: ei_end
    event = EventPhases{ei};
    
    subphsavefolder = fullfile(phsubfolder, event);
    if ~exist(subphsavefolder, 'dir')
        mkdir(subphsavefolder);
    end
    
    for ci = ci_str : ci_end
        pdcond = cond_cell{ci};

        disp([ animal ' ' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz ' event '-' pdcond])

        [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond, 'codesavefolder', savecodefolder);

        % load(and extract) ciCohPhasefile
        ciCohPhasefile = fullfile(savefolder, [ciCohPhasefile_prefix  '_' num2str(f_AOI(1)) '-' num2str(f_AOI(2)) 'Hz_' pdcond '_' event '_align2' align2name '.mat']);
        
        %%% ----------- case no ciCohPhasefile or new Run -------- %%%
        if(~exist(ciCohPhasefile, 'file') || newRun)

            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end

            % lfpseg align2
            [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2,t_AOI, 'codesavefolder', savecodefolder);

            % extract data of chns of AOI
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            lfptrials = lfptrials(mask_chnOfI, :, :);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);

            %  extract and save deltaphis_allChnsTrials and cicoh
            [deltaphis_allChnsTrials, ciCoh, f_selected]= ciCoh_trialDeltaPhi(lfptrials, fs, f_AOI, 'codesavefolder', savecodefolder);
            ntrials = size(lfptrials, 3);
            save(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');

            %  extract and save psedociCohs
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);


            clear files lfptrials fs T_chnsarea mask_chnOfI
            clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected');
        end
        
        %%% ----------- case no psedociCohs variable or psedociCohs nshuffle < shuffleN_psedoTest -------- %%%
        pseciCohsVar = whos('-file',ciCohPhasefile, 'psedociCohs');
        if isempty(pseciCohsVar) || pseciCohsVar.size(4)< shuffleN_psedoTest
            files = dir(fullfile(inputfolder, ['*_' pdcond '_*.mat']));
            if isempty(files)
                continue;
            end
            
            % lfpseg align2
            [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2,t_AOI, 'codesavefolder', savecodefolder);

            % extract data of chns of AOI
            mask_chnOfI = cellfun(@(x) any(strcmp(chnsOfI, x)), T_chnsarea.brainarea);
            lfptrials = lfptrials(mask_chnOfI, :, :);
            T_chnsarea = T_chnsarea(mask_chnOfI, :);

            %  extract and save psedociCohs
            psedociCoh_extract_save(shuffleN_psedoTest, lfptrials, fs, f_AOI, ciCohPhasefile, 'codesavefolder', savecodefolder);
            
            clear files lfptrials fs T_chnsarea mask_chnOfI
        end
        clear pseciCohsVar


        %%% -- plot section --- %%%
        load(ciCohPhasefile, 'deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');
        
        % extract sigciCoh
        [sigciCoh]= sigciCoh_extract(psedociCohs, ciCoh, 'codesavefolder', savecodefolder);


        % extract ciCoh_flatten and chnPairNames, such as M1-stn0-1
        [sigciCoh_flatten, deltaphis_flatten, chnPairNames] = ciCohDephiFlatten_chnPairNames_extract(sigciCoh, deltaphis_allChnsTrials, T_chnsarea, 'codesavefolder', savecodefolder);
        [ciCoh_flatten_used, deltaphis_flatten_used, chnPairNames_used]= ciCoh_deltaPhi_Used(chnPairNames, sigciCoh_flatten, deltaphis_flatten, [], 'codesavefolder', savecodefolder);


        % plot and save ciCoh Histogram image
        if runCicohHist
            nshuffle = size(psedociCohs, 4);
            titlename = [animal '-'  pdcond '-'  event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s,' ' align2 = ' align2name ', ntrials = ' num2str(ntrials) ' nshuffle= ' num2str(nshuffle)];
            
            plot_ciCohHistogram(ciCoh_flatten_used, chnPairNames_used, f_selected, titlename, ...
                'fig_width', 500, 'fig_height', 200, 'codesavefolder', savecodefolder);
            
            saveimgname = [animal '_' event '_' pdcond '_align2' char(align2) '.' image_type];
            saveas(gcf, fullfile(histsavefolder, saveimgname), image_type);
            close gcf
            clear titlename  saveimgname nshuffle
        end


        % rose histogram of deltaphis_allChnsTrials
        if runRosePlot
            titlename_prefix = [animal '-'  pdcond '-'  event];
            subtitlename = [event '['  num2str(t_AOI(1)) ' ' num2str(t_AOI(2))   ']s, align2 = ' char(align2) ', ntrials = ' num2str(ntrials)];
            savefile_suffix = [event '_' pdcond '_align2' char(align2)];
            rosePlotsavefolder = fullfile(savefolder, 'rosePlot', ephase);
            if ~exist(rosePlotsavefolder, 'dir')
                mkdir(rosePlotsavefolder);
            end
            
            plotsave_deltaphirose(deltaphis_flatten_used, ciCoh_flatten_used, chnPairNames_used, f_selected, titlename_prefix, subtitlename, rosePlotsavefolder, savefile_prefix, savefile_suffix, image_type,...
                'codesavefolder', savecodefolder, 'roseRLim', roseRLim);
            close gcf
            clear titlename_prefix subtitlename savefile_prefix savefile_suffix
        end


        clear pdcond subpdsavefolder align2 t_AOI align2name ciCohPhasefile
        clear('deltaphis_allChnsTrials', 'ciCoh', 'T_chnsarea', 'ntrials', 'f_selected', 'psedociCohs');

        close all
    end
end
end

function plot_ciCohHistogram(ciCoh_flatten, chnPairNames, f_selected, titlename, varargin)
%
%   Inputs:
%       ciCoh_flatten:
%       chnPairNames
%       f_selected
%       titlename
%       histClim
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
%           'histClim' - ciCoh histogram clim (default [0 1])
%           'fig_left' - figure position left (default 50)
%           'fig_bottom' - figure position bottom (default 50)
%           'fig_width' - figure position width (default 1200)
%           'fig_height' - figure position height (default 60)
%           'cbarTicks' - vector, default [0 0.5 1]


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
addParameter(p, 'histClim', [0 1], @(x) assert(isnumeric(x) && isvector(x) && length(x)==2));
addParameter(p, 'fig_left', 50, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_bottom', 50, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_width', 1000, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'fig_height', 250, @(x) assert(isnumeric(x) && isscalar(x)));
addParameter(p, 'cbarTicks', [0 0.5 1], @(x) assert(isvector(x) && isnumeric(x)));

parse(p,varargin{:});
codesavefolder = p.Results.codesavefolder;
histClim = p.Results.histClim;
fig_left = p.Results.fig_left;
fig_bottom = p.Results.fig_bottom;
fig_width = p.Results.fig_width;
fig_height = p.Results.fig_height;
cbarTicks = p.Results.cbarTicks;


% copy code to savefolder if not empty
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



% plot
dispix_outinpos = [80 30 15 40]; % left, top, right, bottom of distance between outer and inner position
margin_outpos = [5 5 5 5]; % left, top, right, bottom margin of outer position
outpos = [margin_outpos(1) margin_outpos(4) fig_width-margin_outpos(1)-margin_outpos(3) fig_height-margin_outpos(2)-margin_outpos(4)];
innerpos = [outpos(1)+dispix_outinpos(1) outpos(2)+dispix_outinpos(4) outpos(3)-dispix_outinpos(1)-dispix_outinpos(3) outpos(4)-dispix_outinpos(2)-dispix_outinpos(4)];
figure;
set(gcf, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
imagesc(ciCoh_flatten)
colormap(jet)
set(gca, 'Units', 'pixels');
set(gca, 'OuterPosition', outpos, 'Position', innerpos)
[npairs, nf] = size(ciCoh_flatten);
xticks([1:nf])
xticklabels(round(f_selected))
yticks([1:npairs]);
set(gca,'YTickLabel',chnPairNames,'fontsize',10,'FontWeight','normal')
xlabel('freqs/Hz')
title(titlename, 'FontSize', 10, 'FontWeight', 'normal')
set(gca,'CLim', histClim)
c = colorbar;
c.Label.String = 'ciCoh';
if ~isempty(cbarTicks)
    c.Ticks = cbarTicks;
end

chnPair_prev = '';
for ci = 1: length(chnPairNames)
    chnPair = chnPairNames{ci};
    
    % replace M1-stn0-1 to M1-STN
    s_stn = regexp(chnPair, 'stn[0-9]*-[0-9]*', 'match');
    if ~isempty(s_stn)
        for si = 1 : length(s_stn)
            chnPair = strrep(chnPair, s_stn{si}, 'STN');
        end
    end
    % replace M1-stn0-1 to M1-STN
    s_gp = regexp(chnPair, 'gp[0-9]*-[0-9]*', 'match');
    if ~isempty(s_gp)
        for si = 1 : length(s_gp)
            chnPair = strrep(chnPair, s_gp{si}, 'GP');
        end
    end
    
    if ~strcmp(chnPair_prev, '') && ~strcmp(chnPair_prev, chnPair) % a new site pairs
        hold on; plot(gca, xlim, [(ci + ci -1)/2 (ci + ci -1)/2], 'w--')
        % Create line
    end
    chnPair_prev = chnPair;
    
    clear s_stn s_gp chnPair
end
end


function [lfptrials, fs_lfp, T_chnsarea] = lfpseg_selectedTrials_align2(files, align2, tdur_trial, varargin)
% extract lfp seg data respect to targetonset, reachonset, reach and returnonset separately
% [lfptrials, fs, T_chnsarea] = lfpseg_selectedTrials_align2PeakV(files, [t_AOI(1) t_AOI(2)], 'codesavefolder', savecodefolder);
%
%   not include trials with t_reach <0.2s
% 
%         Args:
%             align2: the event to be aligned 
% 
%             tdur_trial: the duration of extracted trials respected to event(e.g. [-0.5 0.6])
%             
%
%       Name-Value: 
%           'codesavefolder' - code saved folder
% 
%         return:
%             lfptrials: nchns * ntemp * ntrials
% 
%             chnAreas:
% 
%             fs:


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end

coli_align2 = uint32(align2);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);

t_minmax_reach = 0.2;

load(fullfile(files(1).folder, files(1).name),  'fs_lfp', 'T_chnsarea');

nfiles = length(files);
lfptrials = [];
for filei = 1 : nfiles
    
    % load data, lfpdata: [nchns, ntemps, ntrials]
    filename = files(filei).name;
    load(fullfile(files(filei).folder, filename), 'lfpdata', 'T_idxevent_lfp', 'selectedTrials', 'T_idxevent_ma', 'smoothWspeed_trial', 'fs_ma');
    
    if(height(T_idxevent_lfp) == 1)
        disp([filename ' has only 1 trial, skip!']);
        continue;
    end
    
    
    ntrials = length(lfpdata);
    for tri = 1: ntrials
        
        % ignore trials marked with 0
        if ~selectedTrials(tri)
            continue
        end
        
        % select trials based on reach duration
        t_reach = (T_idxevent_lfp{tri, coli_reach} - T_idxevent_lfp{tri, coli_reachonset}) / fs_lfp;
        if t_reach < t_minmax_reach 
            clear t_reach
            continue
        end
        
        if align2 == SKTEvent.PeakV
            % find peakV and its timepoint
            idx_reachonset_ma = T_idxevent_ma{tri, coli_reachonset};
            idx_reach_ma = T_idxevent_ma{tri, coli_reach};
            [~, idx] = max(smoothWspeed_trial{tri}(idx_reachonset_ma: idx_reach_ma, 1));
            idx_peakV_ma = idx + idx_reachonset_ma -1;
            t_reachonset2peakV = (idx_peakV_ma - idx_reachonset_ma)/ fs_ma;
            t_peakV2reach = (idx_reach_ma - idx_peakV_ma)/ fs_ma;
            
            if t_reachonset2peakV < abs(tdur_trial(1)) || t_peakV2reach < tdur_trial(2)
                clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
                clear t_reachonset2peakV  t_peakV2reach
                continue;
            end
            
            % extract trial with t_dur
            idx_peakV_lfp = round(idx_peakV_ma / fs_ma * fs_lfp);
            idx_time0 = idx_peakV_lfp;
            
            clear idx idx_reachonset_ma idx_reach_ma idx_peakV_ma
            clear t_reachonset2peakV  t_peakV2reach
            clear idx_peakV_lfp
        else
            idx_time0 = T_idxevent_lfp{tri, coli_align2}; 
        end
        
        
        % extract phase for 1 trial
        lfp_1trial = lfpdata{tri};
        idxdur = round(tdur_trial * fs_lfp) + idx_time0;
        if idxdur(1) == 0
            idxdur(1) = 1;
        else
            idxdur(1) = idxdur(1) + 1;
        end
        lfp_phase_1trial = lfp_1trial(:,idxdur(1) :idxdur(2));
           
        % cat into lfptrials
        lfptrials = cat(3, lfptrials, lfp_phase_1trial);
        
        clear t_reach idxdur lfp_phase_1trial lfp_1trial
    end
end
end