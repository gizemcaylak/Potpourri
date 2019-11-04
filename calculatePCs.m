function [ PCs, EigVals ] = calculatePCs( X, NumPCs, varargin )
    p = inputParser;
    validMatrix = @(x) validateattributes(x, {'numeric'}, ...
        {'2d', 'nonempty', 'nonsparse', 'real', 'nonnan'});
    validScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'nonempty', 'real', 'nonnan', 'finite', 'nonnegative'});
    validPositiveScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'nonempty', 'real', 'nonnan', 'finite', 'positive'});
    addRequired(p, 'X', validMatrix);
    addRequired(p, 'NumPCs', validScalar);
    addParameter(p, 'BatchSize', 10000, validPositiveScalar);
    parse(p, X, NumPCs, varargin{:});
    [nSample, nVariant] = size(X);
    assert(p.Results.NumPCs < nVariant, ['Expected ''NumPCs'' ', ...
        'to be lesser than the number of columns in ''X''.']);
    if(NumPCs == 0)
       PCs = zeros(size(X, 1), 0);
       EigVals = zeros(1, 0);
    else
        if(p.Results.BatchSize >= nVariant)
            Xn = double(X);
            Xn = (Xn - mean(Xn, 1)) ./ std(Xn, 1);
            A = (Xn * Xn');
        else
            index = 0;
            A = zeros(nSample, nSample);
            while(index < nVariant)
                nItem = min(p.Results.BatchSize, nVariant - index);
                indices = index + (1:nItem);
                Xn = double(X(:, indices));
                Xn = (Xn - mean(Xn, 1)) ./ std(Xn, 1);
                A = A + Xn * Xn';
                index = index + nItem;
            end
        end
        [V, EigVals] = eig(A);
        [EigVals, ind] = sort(diag(EigVals), 'descend');
        PCs = V(:, ind(1:NumPCs));
    end
end

