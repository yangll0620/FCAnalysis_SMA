function m0_restData_uniqdatesYingUsed_extract()

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

% datafolder
[datafolder, ~, ~, ~] = exp_subfolders();



%% input setup
% dates used by Ying
filesYingUsed = fullfile(datafolder, 'DateYingUsed_Pinky_ForYLL.mat');

%% save setup
savefolder = codecorresfolder;


%% Start Here

load(filesYingUsed)
for i = 1: size(dateYingAnalyzed,1)
   onedate = string(dateYingAnalyzed(i, :));
   if i==1
        uniqDateBlocks = onedate;
        predate = onedate;
   else
       if(predate ~= onedate)
           uniqDateBlocks = cat(1, uniqDateBlocks, onedate);
           predate = onedate;
       end      
   end    
   clear onedate
end
save(fullfile(savefolder, 'uniqFileYingUsed_Pinky.mat'), 'uniqDateBlocks');