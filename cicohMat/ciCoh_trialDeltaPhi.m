function [ciCoh, f_selected]= ciCoh_extract(lfptrials, fs, f_AOI, varargin)
%
%
%  Inputs:
%       lfptrials
%       fs:
%       f_AOI
%   
%       Name-Value: 
%           'codesavefolder' - code saved folder       
%   
%   Output:
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials
%       ciCoh:  nchns * nchns * nf
%       f_selected: nf * 1

% parse params
p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);
parse(p,varargin{:});

% copy code to savefolder if not empty
codesavefolder = p.Results.codesavefolder;
if ~isempty(codesavefolder) 
    copyfile2folder(mfilename('fullpath'), codesavefolder);
end



% extract ciCoh
[ciCoh, f_selected] = ciCohSKTAllchns_FFT_NoAmp(lfptrials, fs, f_AOI, 'codesavefolder', codesavefolder);