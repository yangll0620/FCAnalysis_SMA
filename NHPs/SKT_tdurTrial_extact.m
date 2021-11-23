function [tdur_trial_normal, tdur_trial_mild, tdur_trial_moderate] = SKT_tdurTrial_extact(animal, varargin)
%   Inputs:
%       animal
%
%       Name-Value: 
%           'codesavefolder' - code saved folder


% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end


tdur_trial_normal = [];
tdur_trial_mild = [];
tdur_trial_moderate = [];

if strcmpi(animal, 'bug')
    tdur_trial_normal = [-1 1];
    tdur_trial_mild = [-1 1];
end
if strcmpi(animal, 'jo')
    tdur_trial_normal = [-1 1];
    tdur_trial_mild = [-1 1];
    tdur_trial_moderate = [-1 1];
    
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    tdur_trial_normal = [-1 1];
    tdur_trial_moderate = [-1 1];
end

if strcmpi(animal, 'pinky')
    tdur_trial_normal = [-1 1];
    tdur_trial_mild = [-1 1];
end