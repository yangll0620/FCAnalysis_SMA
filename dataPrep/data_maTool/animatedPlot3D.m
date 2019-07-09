function animatedPlot3D(dat3d)
% Functionality: 
%       show the trajectory of 3d data
%
%   Input:
%       dat3d:  3d data (nsamples * 3 (x, y, z))
%       joint_name:  the name of one joint, can be found in .trc file (i.e. 'Shoulder', 'Elbow', 'Wrist' (default), 'Kluver')

x = dat3d(:,1);
y = dat3d(:,2);
z = dat3d(:,3);

% decide the figure x, y, z lims
v_max = ceil(max(dat3d,[],1));
v_min = floor(min(dat3d,[],1));
x_minmax = [v_min(1), v_max(1)];
y_minmax = [v_min(2), v_max(2)];
z_minmax = [v_min(3), v_max(3)];

% Animated Plot in 3D
close all
figure
curve = animatedline('LineWidth', 2);
xv_lim = x_minmax;
yv_lim = y_minmax;
zv_lim = z_minmax;
set(gca, 'XLim', xv_lim, 'YLim', yv_lim, 'ZLim', zv_lim);
view(43, 24);
hold on;
for i = 1: length(z)
    addpoints(curve, x(i), y(i), z(i));
    head = scatter3(x(i), y(i), z(i),'filled');
    drawnow
    delete(head);
    pause(0.008)
end
end