function eucvec = eucliser(u, euc_weights, FixedPoint)

% can selectively calculate Euclidean distance from different combinations
% of dimensions

    diffs = u - FixedPoint;
    weighted_diffs = diffs.*euc_weights;
    squares = weighted_diffs.^2;
    sums = sum(squares,2);
    eucvec = sqrt(sums);

end