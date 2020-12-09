function filteredX = bp_filter(X, frebp, fs)
    %% band pass for each channel X: nchns *  ntemp

    [nchns, ntemp] = size(X);

    if ntemp == 1
        X = X';
        [nchns, ntemp] = size(X);
    end

    filteredX = zeros(nchns, ntemp);

    for chi = 1:nchns
        filteredX(chi, :) = filter_bpbutter(X(chi, :), frebp, fs);
    end

end