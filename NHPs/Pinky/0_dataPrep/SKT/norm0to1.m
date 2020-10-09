function normX = norm0to1(X)
% 
% if X is matrix, normalized to [0 1] along the first dimension

[m, n]= size(X);

if n == 1
    X = reshape(X, n, m);
    [m, n]= size(X);
end


minX = min(X,[],2);
maxX = max(X,[],2);

normX = (X - repmat(minX, 1, n))./repmat(maxX, 1, n);