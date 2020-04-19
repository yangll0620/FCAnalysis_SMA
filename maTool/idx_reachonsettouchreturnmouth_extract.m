function [idx_target, idx_reachonset, idx_touch, idx_return, idx_mouth] ...
    = idx_reachonsettouchreturnmouth_extract(idx_opennings, idx_closings, ma_wrist, fs)
%% extract all the indices of reach onset, touch, return and mouth for each trial
%   
%  Input:
%
%   openning: define as the time point the kluver board changing from close to open   
%   closing: define as the time point the kluver board changing from open to close
%
%
%   Input:
%       ma_kluver: timestamps and x, y, z coordinates of marker, ntimes * 4 (timestamp, x, y, z)
%
%   Output:
%       idx_target, idx_reachonset, idx_touch, idx_return, idx_mouth

triali = 0;
for openclosei = 1: length(idx_opennings)
    % temp_z: total interval used to detect reach onset, touch, return and mouth
    temp_z = ma_wrist(idx_opennings(openclosei): idx_closings(openclosei) + fs, 4);
    temp_y = ma_wrist(idx_opennings(openclosei): idx_closings(openclosei) + fs, 3);
    % real time relative distance in horizontal direction (dim 3)
    disrel_wrist_z = abs(temp_z - temp_z(1));
    disrel_wrist_y = abs(temp_y - temp_y(1));
    speed_z = smooth(abs(diff(temp_z)));
    if max(speed_z) - min(speed_z) < 0.1
        continue;
    end
    maxspeed_z = max(speed_z);
    % discrete speed_z into 0, and 1
    thre_speed = maxspeed_z * 0.1;
    speed_z01 = zeros(size(speed_z));
    speed_z01(speed_z>=thre_speed) = 1;
    
    %% reach onset index and touch index should be in interval [idx_open : idx_close]
    % using dim 3 data, dim 3 is the horizontal direction
    temp1_z = ma_wrist(idx_opennings(openclosei): idx_closings(openclosei), 4);
    % real time relative distance in horizontal direction (dim 3)
    disrel1_wrist_z = abs(temp1_z - temp1_z(1));
    [maxdis1, ~] = max(disrel1_wrist_z);
        
    % velocity in horizontal direction (dim 3)
    speed1_z = smooth(abs(diff(temp1_z)));
    maxspeed1_z = max(speed1_z);
    % discrete speed_z into 0, and 1
    thre1_speed = maxspeed1_z * 0.1;
    speed1_z01 = zeros(size(speed1_z));
    speed1_z01(speed1_z>=thre1_speed) = 1;
    
    % speed status of wrist marker
    status1_z = [0; diff(speed1_z01)]; % combine 0 as diff starts from the second point
    
    % status == 1: from still to movement, status == -1: from movement to still, status == 0: keep still or keep movement
    idx1_movestrs = find(status1_z == 1);
    idx1_movestops = find(status1_z == -1);
    
    [~, locs1_pkspeed] = findpeaks(speed1_z, 'MinPeakHeight', maxspeed1_z * 0.5);    
    
    %% detect the reachonset and touch time point
    % reachonset: identify as the one of the idx_movestrs (from still to movement)
    % touch: identify as the one of the idx_movestops (from movement to still)
    % max(disrel_wrist_z(idx_movestr :idx_movestop)) > maxdis * 0.8, meaning that in interval [idx_movestr:idx_movestop], and hand reaches around the food
    tag_reachonset_touch = 0;
    idi_pkspeed = 1;
    while(idi_pkspeed <= length(locs1_pkspeed))
        % find the idx for movestr and move stop with a peak speed in between
        idices_movestr = find(idx1_movestrs < locs1_pkspeed(idi_pkspeed));
        if isempty(idices_movestr)
            idi_pkspeed = idi_pkspeed + 1;
            continue
        end
        idx_movestr = idx1_movestrs(idices_movestr(end));
        idices_movestop = find(idx1_movestops >= locs1_pkspeed(idi_pkspeed));
        
        for idj = 1: length(idices_movestop)
            idx_movestop = idx1_movestops(idices_movestop(idj));
            
            % index in interval [idx_movestr:idx_movestop], and reach around the food
            if max(disrel1_wrist_z(idx_movestr :idx_movestop)) > maxdis1 * 0.8
                tag_reachonset_touch = 1;
                idxrel_reachonset = idx_movestr;
                idxrel_touch = idx_movestop;
                break;
            end
            clear idx_movestop
        end
        
        if tag_reachonset_touch == 1
            break;
        end
        
        idi_pkspeed = idi_pkspeed + 1;
        clear idices_movestr idices_movestop idj
        clear idx_movestr
    end
    if tag_reachonset_touch == 0
        disp(['idx = ' num2str(openclosei)  ', Can not find reach onset and touch at topen = ' num2str(idx_opennings(openclosei)/100) 's']);
        continue;
    end
        
    %% detect mouth index
    % mouth index should be after touch index
    temp2_z = ma_wrist(idx_opennings(openclosei) + idxrel_touch: idx_closings(openclosei) + fs, 4);
    temp2_y = ma_wrist(idx_opennings(openclosei) + idxrel_touch: idx_closings(openclosei) + fs, 3);
    temp1_y = ma_wrist(idx_opennings(openclosei): idx_closings(openclosei), 3);
    % real time relative distance respect to touch point in vertical direction (dim 2)
    disrel2_wrist_y = abs(temp2_y - temp1_y(idxrel_touch));
    % real time relative distance respect to touch point in horizontal direction (dim 3)
    disrel2_wrist_z = abs(temp2_z - temp1_z(idxrel_touch));
    [~, locs2_pkdis_y] = findpeaks(disrel2_wrist_y);
    tag_mouth = 0;
    thre2_dis = disrel1_wrist_z(idxrel_touch) * 0.6;    
    for loci = 1 : length(locs2_pkdis_y)
        loc = locs2_pkdis_y(loci);
        
        % in horizontal direction back to mouth
        if disrel2_wrist_z(loc) > thre2_dis
            idxrel_mouth = loc;
            tag_mouth = 1;
        end
        
        if tag_mouth == 1
            break;
        end
    end
    if tag_mouth == 0
        disp(['idx = ' num2str(openclosei)  ', Can not find mouth at topen = ' num2str(idx_opennings(openclosei)/100) 's']);
        continue;
    end
    idxrel_mouth = idxrel_mouth + idxrel_touch;
        
    %% Detect return
    % return: identify as the one of the idx_movestrs (from still to movement),
    % idxrel_return should be in interval [idxrel_touch idx_close]
    iopts_return = intersect(find(idx1_movestrs>idxrel_touch), find(idx1_movestrs<idxrel_mouth));
    
    if isempty(iopts_return)
        disp(['idx = ' num2str(openclosei)  ', Can not find return at topen = ' num2str(idx_opennings(openclosei)/100) 's']);
        continue;
    end
    iopt = length(iopts_return);
    tag_return = 0;
    while(iopt >0)
        idx_movestr = idx1_movestrs(iopts_return(iopt));
        
        % disrel1_wrist_z(idxrel_return) nearly equal disrel1_wrist_z(idxrel_touch)
        if disrel_wrist_z(idx_movestr) - disrel_wrist_z(idxrel_touch) < disrel1_wrist_z(idxrel_touch) * 0.1
            idxrel_return = idx_movestr;
            tag_return = 1;
            break;
        end
        iopt = iopt - 1;
        clear idx_movestr
    end
    if tag_return == 0
        disp(['idx = ' num2str(openclosei)  ', Can not find return at topen = ' num2str(idx_opennings(openclosei)/100) 's']);
        continue;
    end    
        
    
    triali = triali + 1;
    idx_target(triali) = idx_opennings(openclosei);
    idx_reachonset(triali) = idx_opennings(openclosei) + idxrel_reachonset -1;
    idx_touch(triali) = idx_opennings(openclosei) + idxrel_touch -1;
    idx_return(triali) = idx_opennings(openclosei) + idxrel_return -1;
    idx_mouth(triali) = idx_opennings(openclosei) + idxrel_mouth -1;
    
    
    if 1
        %% plot
        figure
        % relative distance respect to idx_open
        subplot(2,1,1);
        plot(disrel_wrist_z)
        hold on
        plot(disrel_wrist_y)
        title('distance relative to idx open')
        % reach onset and touch
        plot(idxrel_reachonset,disrel_wrist_z(idxrel_reachonset),'ro')
        plot(idxrel_touch,disrel_wrist_z(idxrel_touch),'bo')
        % mouth
        plot(idxrel_mouth,disrel_wrist_z(idxrel_mouth),'b*')
        % return
        plot(idxrel_return,disrel_wrist_z(idxrel_return),'r*')
        legend('distance z', 'distance y','reachonset', 'touch', 'mouth', 'return')
        xlabel('index')
        ylabel('distance')
        
        
        % speed
        subplot(2,1,2);
        plot(speed_z)
        hold on
        plot(speed_z01)
        % reach onset and touch
        plot(idxrel_reachonset,speed_z(idxrel_reachonset),'ro')
        plot(idxrel_touch,speed_z(idxrel_touch),'bo')
        % mouth
        plot(idxrel_mouth,speed_z(idxrel_mouth),'b*')
        % return
        plot(idxrel_return,speed_z(idxrel_return),'r*')
        legend('speed z','discrete speed z' ,'reachonset', 'touch', 'mouth', 'return')
        xlabel('index')
        ylabel('speed')
        title('speed relative to idx open')
    end
    
    clear temp_z temp_y disrel_wrist_z disrel_wrist_y 
    clear speed_z maxspeed_z thre_speed speed_z01
    clear temp1_z disrel1_wrist_z maxdis1 speed1_z maxspeed1_z thre1_speed speed1_z01 
    clear status1_z idx1_movestrs idx1_movestops locs1_pkspeed
    clear tag_reachonset_touch idi_pkspeed
    clear temp2_z temp2_y temp1_y
    clear disrel2_wrist_y disrel2_wrist_z thre2_dis locs2_pkdis_y 
    clear tag_mouth loci
    clear iopts_return 
    clear tag_return iopt
end
