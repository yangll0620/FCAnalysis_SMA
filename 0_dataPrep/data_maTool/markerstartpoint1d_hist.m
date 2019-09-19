function stpoint_marker= markerstartpoint1d_hist(ma_marker1d)
% get the start point of marker in 1d using histogram 
% 
%   Method:
%       the mean of one histogram bin 
%
%   Input:
%       ma_marker: coordinates of marker, ntimes * 1 
%
%   Output:
%       stpoint_marker: the start point of marker
%

[N,~,BIN] = histcounts(ma_marker1d);
[maxbin] = find(N == max(N,[],'all'));
idx = find(BIN == maxbin);
stpoint_marker= mean(ma_marker1d(idx, :),1); % stpoint_marker is the start point of joint