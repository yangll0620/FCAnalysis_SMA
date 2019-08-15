function animatedPlotMA(file_matrc, marker_name)
% Functionality: 
%       show the trajectory of MA data for marker_name in file_matrc file
%
%   Example:
%       file_matrc = fullfile('F:','yang7003@umn','NMRC_umn', 'Projects', ...
%           'NWBStandardization','workingfolders','home','data_shared','raw',...
%           'bug','expdata', 'setupchair','bug-190111', 'ma','Bug_20190111_3_cleaned.trc');
% 
%       animatedPlotMA(file_matrc, 'Kluver')
%
%   Input:
%       file_matrc:  the full path of cleaned ma .trc file
%       marker_name:  the name of one marker, can be found in .trc file (i.e. 'Shoulder', 'Elbow', 'Wrist' (default), 'Kluver')
%
%
if nargin < 2
    marker_name = 'Wrist';
end

% extract ma_marker: the time stamp and x, y, z coordinates of marker, ntimes * 4 (timestamp, x, y, z)
ma_marker = mamarkerdata_extract(file_matrc, marker_name); 

% show the trajectory
x = ma_marker(:,2);
y = ma_marker(:,3);
z = ma_marker(:,4);

% decide the figure x, y, z lims
v_max = ceil(max(ma_marker(:,2:4),[],1));
v_min = floor(min(ma_marker(:,2:4),[],1));
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
t_startout = tic;
for i = 1: length(z)
    t_start = tic;
    addpoints(curve, x(i), y(i), z(i));
    head = scatter3(x(i), y(i), z(i),'filled');
    drawnow
    delete(head);
    pause(0.008)
    times(i) = toc(t_start);
end
disp(mean(times))
toc(t_startout)
end

