function [datafolder, codefolder, pipelinefolder, outputfolder] = exp_subfolders()
%% folder_extract() return absoluate folder path
% 
%     Usage:
%         
%         [datafolder, codefolder, pipelinefolder, outputfolder] = exp_subfolders()
% 
%     Returns:
%         datafolder: data folder 
%         
%         codefolder: code folder 
%         
%         pipelinefolder: pipeline folder
%         
%         outputfolder: output folder
% 

currfile = mfilename('fullpath');

% the exp folder path
expfolder = currfile(1 : strfind(currfile, 'exp')+length('exp') -1);

% data folder
datafolder = fullfile(expfolder, 'data');

% code folder
codefolder = fullfile(expfolder, 'code');

% pipeline folder
pipelinefolder = fullfile(expfolder, 'pipeline');

% output folder
outputfolder = fullfile(expfolder, 'output');

