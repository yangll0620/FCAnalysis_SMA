function resampledX = resampleSeg(X, fs_old, fs_new)
%% downsample along the second dim X:  nchns * ntemp * ntrials
%
%   return:
%       resampledX: nchns * ntemp(resampled) * ntrials

resampledX = [];
for triali = 1: size(X, 3)
    x = X(:,:,triali); % x: nchns * ntemp 
    resamp = resample(x', round(fs_new), round(fs_old)); % resamp: ntemp(resampled) * nchns
    resamp = resamp'; % resamp: nchns * ntemp(resampled)
    
    resampledX = cat(3, resampledX, resamp); % resampledX: nchns * ntemp(resampled) * ntrials
end