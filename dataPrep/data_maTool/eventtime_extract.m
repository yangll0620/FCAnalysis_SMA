function eventtime_extract()
% extract the event time points of target onset, reach onset, touch, return
% onset and mouth.
%
animal = 'Pinky';
skbdatebk_file = ['F:\yang7003@umn\NMRC_umn\Projects\FCAnalysis\metainf\' animal '\pinky_skbinf.csv'];
inf_skbdatebk = readtable(skbdatebk_file);
coli_eventtimeextracted = find(cellfun(@(x) strcmp(x, 'eventtimeextracted'), inf_skbdatebk.Properties.VariableNames));
if isempty(coli_eventtimeextracted)
    eventtimeextracted = cell(size(inf_skbdatebk, 1),1);
    inf_skbdatebk = addvars(inf_skbdatebk,eventtimeextracted,'NewVariableNames','eventtimeextracted');
else
    inf_skbdatebk = readtable(skbdatebk_file);
end

addpath(fullfile('..','NHP_Pinky'))

NMRCdriver = 'Y:';
for i = 1:1%1: size(inf_skbdatebk, 1)
    
    %ma trc file folder and filename
    bkma = inf_skbdatebk.bkma(i);
    dateexp = datevec(num2str(inf_skbdatebk.dateexp(i)),'yymmdd');
    mafolder = fullfile(NMRCdriver, 'Animals2', 'Pinky', 'Recording', 'Raw', 'rawMA', ['MA' datestr(dateexp,'yyyymmdd')]);
    file = dir(fullfile(mafolder, ['pinky*' datestr(dateexp,'yyyymmdd') '_' num2str(bkma) '_cleaned.trc']));
    if length(file) ~= 1
        disp(['i = ' num2str(i) ' in folder ' mafolder ' find more than one or no matrc files!'])
        continue;
    end
    bktdt = inf_skbdatebk.bktdt(i);
    disp(['SKB time points for ' animal '_' datestr(dateexp,'mmddyy') '-mabk' num2str(bkma) '-tdtbk' num2str(bktdt)])
    file_matrc = fullfile(file.folder,file.name);
        
    % extract the x, y  and z coordinates for marker kluver and wrist, nan data were removed
    % ma_kluver/ ma_marker: time stamp and x, y, z coordinates of marker, ntimes * 4 (timestamp, x, y, z)
    [ma_kluver] = mamarkerdata_extract(file_matrc, 'Kluver');
    if ma_kluver == 0
        disp(['No ma data of kluver marker!!']);
        continue;
    end
    [ma_wrist, fs] = mamarkerdata_extract(file_matrc, 'Wrist');
    if ma_wrist == 0
        disp(['No ma data of kluver marker!!']);
        continue;
    end
    
    % extract the time points exist in both marker kluver and wrist
    % time stamps of marker wrist and marker kluver
    t_wrist = ma_wrist(:,1);
    t_kluver = ma_kluver(:,1);
    % time stamps for both marker wrist and kluver
    [~, idx_wrist, idx_kluver] = intersect(t_wrist, t_kluver);
    % synchronized ma_wrist and ma_kulver
    ma_wrist = ma_wrist(idx_wrist, :);
    ma_kluver = ma_kluver(idx_kluver, :);
    
    % openning and closing indices of kluver board extraction
    [idx_opennings, idx_closings, dis_rel] = idx_opencloseboard_extract(ma_kluver, fs);
    
    
    % reach onset, touch, return and mouth indices extraction
    [idx_target, idx_reachonset, idx_touch, idx_return, idx_mouth] ...
    = idx_reachonsettouchreturnmouth_extract(idx_opennings, idx_closings, ma_wrist, fs);
    
    t_target = idx_target' / fs;
    t_reachonset = idx_reachonset' / fs;
    t_touch = idx_touch' / fs;
    t_return = idx_return' / fs;
    t_mouth = idx_mouth' / fs;
    tbl_eventtime = table(t_target, t_reachonset, t_touch,t_return, t_mouth);
    
    savefolder = fullfile(NMRCdriver, 'Animals2', animal, 'Recording', 'Processed', 'DataDatabase', ...
        [animal '_' datestr(dateexp,'mmddyy')], ['Block-' num2str(bktdt)]);
    skbtime_filename = ['MASKBeventtime_' animal '_' datestr(dateexp,'yyyymmdd') '_mabk' num2str(bkma) '_tdtbk' num2str(bktdt) '.mat'];
    save(fullfile(savefolder, skbtime_filename), 'tbl_eventtime',...
        't_target', 't_reachonset', 't_touch', 't_return', 't_mouth');
    inf_skbdatebk.eventtimeextracted(i) = {'Yes'};
    inf_skbdatebk.PDCondition(i) = {parsePDCondtion_Pinky(datenum(dateexp))};
    
    clear bkma dateexp mafolder file file_matrc
    clear bktdt savefolder skbtime_filename 
    clear ma_kluver ma_wrist fs
    clear t_wrist t_kluver idx_wrist idx_kluver
    clear idx_opennings idx_closings dis_rel
    clear idx_target idx_reachonset idx_touch idx_return idx_mouth
    clear t_target t_reachonset t_touch t_return t_mouth tbl_eventtime
end
writetable(inf_skbdatebk, skbdatebk_file);
clear skbdatebk_file animal inf_skbdatebk NMRCdriver i


