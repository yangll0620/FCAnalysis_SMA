function depth_filled = filled_dailyDepth(depth_1chn, iniDepth)
% filled the depth, if empty, use the previous one, 
%                   otherwise depth_1chn(i) + initDepth
% 
% Input
%       depth_1chn: the column Number ndays * 1
%       iniDepth: initial depth, scale
% output:
%       depth_filled: filled depth vector ndays * 1

depth_filled = zeros(size(depth_1chn));

% fill the first day
i = 1;
if isnan(depth_1chn(i))
    depth_filled(i) = iniDepth;
else
    depth_filled(i) = iniDepth + depth_1chn(i);
end

for i = 2: length(depth_1chn)
    
    if isnan(depth_1chn(i))
        depth_filled(i) = depth_filled(i-1);
        
    else
        depth_filled(i) = iniDepth + depth_1chn(i);
    end
end
