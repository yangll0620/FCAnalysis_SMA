function stpoint_joint = jointstartpoint_hist(ma_joint)
% get the start point of joint
% 
%   Method:
%       the mean of one histogram bin which has the largest counts
%
%   Input:
%       ma_joint: x, y, z coordinates of joint, ntimes * 3 (x, y, z)
%
%   Output:
%       stpoint_joint: the start point of joint
%   
[N,XEDGES,YEDGES,BINX,BINY] = histcounts2(ma_joint(:,1), ma_joint(:,2));
[idx_row, idx_col] = find(N == max(N,[],'all'));
idx1 = find(BINX == idx_row);
idx2 = find(BINY == idx_col);
idx = intersect(idx1, idx2);
stpoint_joint = mean(ma_joint(idx, :),1); % stpoint_joint is the start point of joint