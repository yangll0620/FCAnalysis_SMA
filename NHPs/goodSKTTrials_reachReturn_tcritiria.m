function [t_minmax_reach_normal, t_minmax_return_normal, t_minmax_reach_mild, t_minmax_return_mild, t_minmax_reach_moderate, t_minmax_return_moderate] = goodSKTTrials_reachReturn_tcritiria(animal, varargin)
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


if strcmpi(animal, 'bug')
    t_minmax_reach_normal = [0.5, 2];
    t_minmax_return_normal = [0.5, 2];
    t_minmax_reach_mild = [0.5, 2];
    t_minmax_return_mild = [0.5, 2];
    
    t_minmax_reach_moderate = [0.5 2];
    t_minmax_return_moderate = [0.5 2];
    
end
if strcmpi(animal, 'jo') 
    t_minmax_reach_normal = [0.5, 1];
    t_minmax_return_normal = [0.4, 1];
    t_minmax_reach_mild = [0.6 1];
    t_minmax_return_mild = [0.8 1.3];
    t_minmax_reach_moderate = [0.6 1];
    t_minmax_return_moderate = [0.8 1.4];
    
end

if strcmpi(animal, 'kitty') % Kitty not have mild
    t_minmax_reach_normal = [0.5, 10];
    t_minmax_return_normal = [0.5, 10];
    t_minmax_reach_moderate = [0.5, 10];
    t_minmax_return_moderate = [0.5, 10];
    
    t_minmax_reach_mild = [];
    t_minmax_return_mild = [];
    
end

if strcmpi(animal, 'pinky')
    t_minmax_reach_normal = [0.5, 1];
    t_minmax_return_normal = [0.5, 1];
    t_minmax_reach_mild = [0.5 1];
    t_minmax_return_mild = [0.5 1];
    t_minmax_reach_moderate = [0.7 1.2];
    t_minmax_return_moderate = [0.8 1.2];
    
end