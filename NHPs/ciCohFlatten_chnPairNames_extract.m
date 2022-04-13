function [ciCoh_flatten, chnPairNames] = ciCohFlatten_chnPairNames_extract(ciCoh, T_chnsarea, varargin)
%
%   Input:
%       ciCoh: nchns * nchns * nf
%       deltaphis_allChnsTrials: nchns * nchns * nf * ntrials(nsegs)
%       T_chnsarea: nchns * 5
%
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


[nchns, ~, nf] = size(ciCoh);
chnPairNames = {};
ciCoh_flatten = zeros(nchns * (nchns -1)/2, nf);
ci = 0;
for chni = 1 : nchns -1
    for chnj = chni + 1  : nchns
        chnPairNames = [chnPairNames; {[T_chnsarea.brainarea{chni} '-'  T_chnsarea.brainarea{chnj}]}];
        
        ci = ci + 1;
        ciCoh_flatten(ci, :) = ciCoh(chni, chnj, :);
    end
end