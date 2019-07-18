function [idx_opennings, idx_closings, dis_rel] = idx_opencloseboard_extract(ma_kluver, fs)
%% extract the indix of openning and closing kluver board for each trial
%   openning: define as the time point the kluver board changing from close to open   
%   closing: define as the time point the kluver board changing from open to close
%
%  Method:
%       1. calculate the relative distance of marker kluver in each time stamp respect to
%         start point of the kluver board, start point is detected using 3d histogram 
%
%       2. use the average of  the peak relative distances as the maximal
%          relative distance of marker kluver, pkdis_kluver = mean(peaks_dis)
%
%       3. set shreshold for the status of kluver board openning and open success
%           thre_open = pkdis_kluver * 0.05;
%           thre_openSuccess = pkdis_kluver  * 0.8;
%
%       4. discrete the distance into 0 and 1 base on thre_open
%
%       5. extract the index where from close to open and from open to close 
%
%   Input:
%       ma_kluver: x, y, z coordinates of marker, ntimes * 3 (x, y, z)
%
%   Output:
%       idx_opennings, idx_closings: indices of openning and closing kluver
%       board for each trial (1 * ntrials)
%
%       dis_rel: relative distance of marker kluver respect to the stpoint_kluver


%% extract the indix of openning and closing kluver board for each trial
% real time distance of marker kluver relative to start point of kluver: dis_rel 
%  extract threshold for openning the kluver and discrete dis_rel into 0 and 1
% get the start point of kluver marker in 3d using histogram 
stpoint_kluver = markerstartpoint3d_hist(ma_kluver(:,2:4));
% real time distance
dis_rel = sqrt(sum((ma_kluver(:,2:4) - repmat(stpoint_kluver, [size(ma_kluver(:,2:4),1),1])).^2,2));
% smooth distance
dis_rel = smooth(dis_rel);

% find the peaks of relative distance dis_rel, criterial: MinPeakDistance = 3s, 
% MinPeakHeight = pkdis_kluver * 0.8
peaks_dis = findpeaks(dis_rel, 'MinPeakDistance', fs * 3);
peaks_dis = sort(peaks_dis,'descend');
if length(peaks_dis) >= 7
    pkdis_kluver = mean(peaks_dis(3:7));
else
    if length(peaks_dis) <=2
        disp('length(peaks_dis) <=2');
    else
        pkdis_kluver = mean(peaks_dis(3:end));
    end
end

% threshold for openning the kluver board is the 10% of peaks
thre_open = pkdis_kluver * 0.05;
thre_openSuccess = pkdis_kluver  * 0.8;

% discrete the distance into 0 and 1 base on thre_open
dis_discrete(dis_rel < thre_open) = 0;
dis_discrete(dis_rel >= thre_open) = 1;

%% extract the index where from close to open and from open to close 
% close and open status of kluver board
status = [0 diff(dis_discrete)]; % combine 0 as diff starts from the second point

% status == 1: from close to open, status == -1: from open to close 
idx_opennings = find(status == 1);
idx_closings = find(status == -1);

% refine the idx_opennings and idx_closings
idx_opennings1 = idx_opennings;
idx_closings1 = idx_closings;
clear idx_opennings idx_closings
triali = 0;
for i = 1: length(idx_opennings1)
    idx_open = idx_opennings1(i);
    idices = find(idx_closings1 > idx_open);
    if ~isempty(idices)
        idx_close = idx_closings1(idices(1));
        if max(dis_rel(idx_open:idx_close)) >= thre_openSuccess
            triali = triali + 1;
            idx_opennings(triali) = idx_open;
            idx_closings(triali) = idx_close;
        end
    end
end
