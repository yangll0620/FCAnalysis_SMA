function m0_restData_dateBlockYingUsed()
%   Manually marked the good and bad segments
% 
% 1. add variable segsRemain, for marking each segment with 1 (good) or 0 (not good)
%


%% folder generate
% the full path and the name of code file without suffix
codefilepath = mfilename('fullpath');


% find the codefolder
idx = strfind(codefilepath, 'code');
codefolder = codefilepath(1:idx + length('code')-1);
clear idx


% add util path
addpath(genpath(fullfile(codefolder,'util')));


% the corresponding pipeline folder for this code
[codecorresfolder, ~] = code_corresfolder(codefilepath, true, false);



%% save setup
savefolder = codecorresfolder; 


%% Start Here
load('/home/lingling/root2/Ying Yu/SKB_EventRelatedBeta/datasets_Allanimals/Resting/data_segment_Jo_Rest_NonMovwithMAbyCleandataNo60HzFilt (2).mat');

nfiles = length(data_segments);
for i = 1: nfiles
   
    dateBlock = getfield(data_segments, {i}, 'date');
   

   if i==1
        dateBlocks = dateBlock;
        predate = dateBlock;
   else
       if(~strcmp(predate, dateBlock))
           dateBlocks = cat(1, dateBlocks, dateBlock);
           predate = dateBlock;
       end      
   end    
   clear dateBlock
end

save(fullfile(savefolder, 'dateBlocksYingUsed_rest.mat'), 'dateBlocks');