function goodTrials = check_spectrogram_oneGroup(lfpdata, T_idxevent, T_chnsarea, fs, goodTrials, filename)
% lfpdata: nchns * ntemp * ntrials


% global parameters


coli_targetonset = uint32(SKTEvent.TargetOnset);
coli_reachonset = uint32(SKTEvent.ReachOnset);
coli_reach = uint32(SKTEvent.Reach);
coli_returnonset = uint32(SKTEvent.ReturnOnset);
coli_mouth = uint32(SKTEvent.Mouth);
% align to event
coli_align2 = coli_reachonset;

% subplot/Figure parameters
fig_left = 50;
fig_bottom = 50;
fig_width = 1800;
fig_height = 900;





% position of Next and previous button
btn_width = 50;
btn_height = 30;
pos_btn_next_left = (subp_endLeft) * fig_width;
pos_btn_next_bottom = (0.05) * fig_height;
pos_btn_prev_left = pos_btn_next_left;
pos_btn_prev_bottom = pos_btn_next_bottom + btn_height + 10;
% Finish Button
pos_btn_finish_left = (subp_endLeft) * fig_width;
pos_btn_finish_bottom = 0.5 * fig_height;


% trial number checkbox parameters
cb_width = 150;
cb_height = 50;


[nchns, ~, ntrials] = size(lfpdata);
ngs = ceil(ntrials / subp_ntrials);

fig = figure(); set(fig, 'PaperUnits', 'points',  'Position', [fig_left fig_bottom fig_width fig_height]);
hcbs = zeros(ntrials, 1);
checkedAllGs = zeros(ngs, 1);


gi = 0;
clf(fig);
c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [pos_btn_next_left pos_btn_next_bottom btn_width btn_height]);
c_next.Callback = @nextTrials;
title(['gi = ' num2str(gi + 1) '/' num2str(ngs)])
tri_str = 1;
tri_end = subp_ntrials;
plot_spectrogram()
checkedAllGs(gi + 1) = 1;

uiwait(fig)


    function showNewTrials()
        clf(fig);
        
        if(gi < ngs-1)
            c_next = uicontrol(fig, 'Style','pushbutton', 'String', 'Next', 'Position', [pos_btn_next_left pos_btn_next_bottom btn_width btn_height]);
            c_next.Callback = @nextTrials;
        end
        
        if(gi > 0)
            c_prev = uicontrol(fig, 'Style','pushbutton', 'String', 'Previous', 'Position', [pos_btn_prev_left pos_btn_prev_bottom btn_width btn_height]);
            c_prev.Callback = @prevTrials;
        end
        
        % start and end trial number for next graph
        tri_str = gi * subp_ntrials + 1;
        tri_end = (gi + 1) * subp_ntrials;
        if tri_end > ntrials
            tri_end = ntrials;
        end
        plot_spectrogram()
        
        checkedAllGs(gi + 1) = 1;
        if all(checkedAllGs)
            c_Finish = uicontrol(fig, 'Style','pushbutton', 'String', 'Finish', 'Position', [pos_btn_finish_left pos_btn_finish_bottom btn_width btn_height]);
            c_Finish.Callback = @finishCheck;
        end
        
    end

    function finishCheck(src,event)
        finishedCheck = true;
        close(fig)
    end % finishCheck

    function nextTrials(src,event)
        gi = mod(gi + 1, ngs);
        showNewTrials()
    end % nextTrials

    function prevTrials(src,event)
        gi = mod(gi - 1, ngs);
        showNewTrials()
    end % prevTrials

    function box_value(hObj,event)
        % Called when boxes are used
        
        % find the used box
        idx = find(hcbs==hObj);
        
        goodTrials(idx) =  get(hObj,'Value');
    end

    function plot_spectrogram_trialbyTrial(lfpdata, T_idxevent, T_chnsarea, fs, tri_str, tri_end, coli_align2, coli_mouth,filename, fig_width, fig_height, cb_width, cb_height)
        % plot lfpdata_1group of all the channels: nchns * ntemp * ntrial
        
        twin = 0.2;
        toverlap = 0.15;
        f_AOI = [5 100];
        
        
        areaName_left = 0.003;
        
        % spectrogram subplot parameters
        subp_startLeft = 0.05;
        subp_endLeft = 0.97;
        subp_startTop = 0.95;
        subp_width = 0.2;
        subp_height = 0.12;
        supb_deltaX = 0.01;
        supb_deltaY = 0.01;
        subp_nchns = floor((subp_startTop - subp_height)/(subp_height + supb_deltaY )) + 1;
        subp_ntrials = floor((subp_endLeft - subp_startLeft - subp_width)/(subp_width + supb_deltaX )) + 1;
        
        
        
        [nchns, ~, ntrials] = size(lfpdata);
        
        annotation(gcf,'textbox',...
            [subp_startLeft subp_startTop 1 0.03],...
            'String', {filename}, ...
            'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
        
        for tri = tri_str: tri_end
            coli = mod(tri,subp_ntrials);
            if coli == 0
                coli = subp_ntrials;
            end
            
            % left postion for tri
            subp_left = (coli -1) * (subp_width + supb_deltaX )+ subp_startLeft;
            
            % plot trial number Name and Checkbox
            pos_cb_trialN_left = (subp_left + subp_width/2) * fig_width;
            pos_cb_trialN_bottom = (subp_startTop) * fig_height;
            hcbs(tri) = uicontrol('Style','checkbox','Value', goodTrials(tri),...
                'Position', [pos_cb_trialN_left pos_cb_trialN_bottom cb_width cb_height],'String', ['triali = ' num2str(tri) '/' num2str(ntrials)]);
            set(hcbs(tri),'Callback',{@box_value});

            
            for chi = 1 : nchns
                
                x = lfpdata(chi, 1:T_idxevent{tri, coli_mouth}, tri);
                
                % spectrogram
                [~, freqs, times, psd] = spectrogram(x, round(twin * fs), round(toverlap * fs),[],fs); % psd: nf * nt
                
                % select freqs_plot and corresponding psd
                idx_f = (freqs >= f_AOI(1) & freqs <= f_AOI(2));
                freqs_plot =  freqs(idx_f);
                psd_plot = psd(idx_f, :);
                % convert to db and gaussfilt
                psd_plot = 10 * log10(psd_plot);
                psd_plot = imgaussfilt(psd_plot,'FilterSize',5);
                times_plot = times - T_idxevent{tri, coli_align2} / fs;
                
                
                % spectrogram subplot
                subp_bottom = subp_startTop - subp_height - (chi -1) * (subp_height + supb_deltaY);
                subplot('Position', [subp_left, subp_bottom, subp_width, subp_height])
                imagesc(times_plot, freqs_plot, psd_plot);
                set(gca,'YDir','normal')
                colormap(jet)
                colorbar
                if (chi == nchns)
                    xlabel('time/s')
                else
                    xticks([])
                end
                if(coli == 1)
                    ylabel('Frequency(Hz)')
                else
                    yticks([])
                end
                hold on
                % plot event lines
                for eventi = 1 : width(T_idxevent)
                    t_event = (T_idxevent{tri, eventi} - T_idxevent{tri, coli_align2}) / fs;
                    plot([t_event t_event], ylim, '--', 'LineWidth',1.5)
                    clear t_event
                end
                clear eventi
                
                
                % plot chn Names only in the first shown trial
                if(coli == 1)
                    annotation(gcf,'textbox',...
                        [areaName_left subp_bottom + subp_height / 2 0.03 0.03],...
                        'String',T_chnsarea{chi, 'brainarea'}, ...
                        'LineStyle', 'none', 'FontWeight', 'bold', 'FitBoxToText', 'off');
                end
                
                
                clear x freqs times psd idx_f freqs_plot psd_plot times_plot
                clear subp_bottom
            end
        end
        
    end %plot_spectrogram

end