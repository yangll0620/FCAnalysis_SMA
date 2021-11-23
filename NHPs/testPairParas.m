function testPairParas(varargin)

p = inputParser;
addParameter(p, 'codesavefolder', '', @isstr);

parse(p,varargin{:});

codesavefolder = p.Results.codesavefolder;