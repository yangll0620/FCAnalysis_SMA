function tTrial_strEnd = ExtractTrial_tStrEnd(trialTimeXlsFile, dateofexp_str, bk)
opts = detectImportOptions(trialTimeXlsFile);
opts = setvartype(opts,{'date'},'string');
opts = setvartype(opts,{'bk', 'tStrInVideo', 'tEndInVideo'},'double');
t_trialtime = readtable(trialTimeXlsFile, opts);
clear opts


% extract one date-bk trials
rows = (t_trialtime.date== dateofexp_str & t_trialtime.bk == bk);
t_trialtime_1datebk = t_trialtime(rows, :);
clear rows

% frameTrial_strEnd
tTrial_strEnd = [t_trialtime_1datebk.tStrInVideo  t_trialtime_1datebk.tEndInVideo];
end