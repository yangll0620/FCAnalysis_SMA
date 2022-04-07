function [align2, t_AOI, align2name] = SKT_EventPhase_align2_tAOI_extract(event, animal, pdcond, varargin)
%
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

if strcmpi(event, 'preMove')
    align2 = SKTEvent.TargetOnset;
    t_AOI = [-1 -0.8];
end
if strcmpi(event, 'Anticipated')
    align2 = SKTEvent.TargetOnset;
    t_AOI = [-0.2 0];
end
if strcmpi(event, 'earlyReach')
    align2 = SKTEvent.ReachOnset;
    t_AOI = [0 0.2];
end
if strcmpi(event, 'Return')
    align2 = SKTEvent.ReturnOnset;
    t_AOI = [0 0.2];
end

if strcmpi(event, 'lateReach')
    align2 = SKTEvent.Reach;
    t_AOI = [-0.2 0];
end

if strcmpi(event, 'peakV')
    align2 = 'peakV';
    t_AOI = [-0.1 0.1];
end


align2name = char(align2);

% special case for Kitty
if strcmpi(animal, 'Kitty') && strcmpi(pdcond, 'normal') && strcmpi(event, 'preMove')
    align2 = SKTEvent.TargetOnset;
    t_AOI = [0.2 0.4];
    align2name = 'StartTrial';
end

if strcmpi(animal, 'Kitty') && strcmpi(pdcond, 'normal') && strcmpi(event, 'Anticipated')
    align2 = SKTEvent.ReachOnset;
    t_AOI = [-0.2 0];
    align2name = char(align2);
end

% align2name for Reach is Touch
if strcmpi(align2, 'Reach')
    align2name = 'Touch';
end

