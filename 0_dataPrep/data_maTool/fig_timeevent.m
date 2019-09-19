function fig_timeevent(idx_opennings, idx_closings, fs, idxrel_reachonset,idxrel_touch, idxrel_mouth, idxrel_return)
for openclosei = 1: length(idx_opennings)
    % temp_z: total interval used to detect reach onset, touch, return and mouth
    temp_z = ma_wrist(idx_opennings(openclosei): idx_closings(openclosei) + fs, 4);
    temp_y = ma_wrist(idx_opennings(openclosei): idx_closings(openclosei) + fs, 3);
    % real time relative distance in horizontal direction (dim 3)
    disrel_wrist_z = abs(temp_z - temp_z(1));
    disrel_wrist_y = abs(temp_y - temp_y(1));
    speed_z = smooth(abs(diff(temp_z)));
    
    maxspeed_z = max(speed_z);
    % discrete speed_z into 0, and 1
    thre_speed = maxspeed_z * 0.1;
    speed_z01 = zeros(size(speed_z));
    speed_z01(speed_z>=thre_speed) = 1;
    
    %% plot
    figure
    % relative distance respect to idx_open
    subplot(2,1,1);
    plot(disrel_wrist_z)
    hold on
    plot(disrel_wrist_y)
    title('relative distance repect to idx open')
    % reach onset and touch
    plot(idxrel_reachonset,disrel_wrist_z(idxrel_reachonset),'ro')
    plot(idxrel_touch,disrel_wrist_z(idxrel_touch),'bo')
    % mouth
    plot(idxrel_mouth,disrel_wrist_z(idxrel_mouth),'b*')
    % return
    plot(idxrel_return,disrel_wrist_z(idxrel_return),'r*')
    legend('distance z', 'distance y','reachonset', 'touch', 'mouth', 'return')
    
    
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
end