function stpoint_marker = markerstartpoint3d_hist(ma_marker)
% get the start point of marker in 3d using histogram 
% 
%   Method:
%       the mean of one histogram bin which has the largest counts, the
%       dimension used is the two have the largest std among x, y, z dims
%
%   Input:
%       ma_marker: x, y, z coordinates of marker, ntimes * 3 (x, y, z)
%
%   Output:
%       stpoint_marker: the start point of marker in 3d
%

% find the dim that has the largest std among x, y, z dims
[~, idx_std] = sort([std(ma_marker(:,1)) std(ma_marker(:,2)) std(ma_marker(:,3))],'descend');

[N,~,~,BINX,BINY] = histcounts2(ma_marker(:,idx_std(1)), ma_marker(:,idx_std(2)));
[idx_row, idx_col] = find(N == max(N,[],'all'));
idx1 = find(BINX == idx_row);
idx2 = find(BINY == idx_col);
idx = intersect(idx1, idx2);
stpoint_marker = mean(ma_marker(idx, :),1); % stpoint_marker is the start point of marker