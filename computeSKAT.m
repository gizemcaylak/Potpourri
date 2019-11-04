function [ C, PCs ] = computeSKAT( X, Y, varargin )
    p = inputParser;
    p.CaseSensitive = false;
    validMatrix = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonempty', 'real', 'nonsparse', 'nonnan', 'finite'});
    validColumn = @(x) validateattributes(x, ...
        {'numeric','logical', 'string'}, ...
        {'column', 'nonempty', 'nonsparse', 'real', 'nonnan'});
    validPCs = @(x) validateattributes(x, {'numeric'}, ...
        {'2d', 'real', 'nonsparse', 'nonnan', 'finite'});
    validIntegerScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'nonempty', 'real', 'nonnan', 'nonnegative', 'integer'});
    validPositiveScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', 'nonempty', 'real', 'nonnan', 'positive', 'integer'});
    validResponseVarType = @(x) any(validatestring(x, ...
        {'continuous', 'dichotomous', 'ordinal'}));
    addRequired(p, 'X', validMatrix);
    addRequired(p, 'Y', validColumn);
    addParameter(p, 'PCs', NaN, validPCs);
	addParameter(p, 'NumPCs', 0, validIntegerScalar);
    addParameter(p, 'BatchSize', 10000, validPositiveScalar);
    addParameter(p, 'ResponseVarType', 'continuous', validResponseVarType);
    parse(p, X, Y, varargin{:});
    [nSample, nVariant] = size(X);
    invalidateMismatch(X, Y, 'X', 'Y', 'row');
    assert(p.Results.NumPCs < nVariant, ['Expected ''NumPCs'' ', ...
        'to be lesser than the number of columns in ''X''.']);
    param = p.Results;
    checkUsingDefaults = @(p,varname) any(strcmp(p.UsingDefaults,varname));
    if(checkUsingDefaults(p, 'PCs'))
        PCs = calculatePCs(X, param.NumPCs, 'BatchSize', param.BatchSize);
    else
        PCs = param.PCs;
        invalidateMismatch(X, PCs, 'X', 'PCs', 'row');
    end
    cov = [ones(nSample, 1), PCs];     %% Covariates  - n x (NumPCs + 1) matrix
    param.Ordinal = false;
    switch(lower(param.ResponseVarType))
        case 'continuous'
            y = cov * (cov \ Y);       %% Estimated Y - n x 1 column vector
            r = Y - y;                 %% Residuals   - n x 1 column vector
            param.Categorical = false;
        case 'ordinal'
            param.Coding = 'ordinal';
            param.Categorical = true;
            param.Ordinal = true;
%         case 'nominal'
%             param.Coding = 'onevsone';
%             param.Categorical = true;
        case 'dichotomous'
            param.Coding = 'onevsone';
            param.Categorical = true;
    end
    if(param.Categorical)
        [classNames, ~, Yp] = unique(Y);
        nClass = length(classNames);
        if(nClass > 2 && ~param.Ordinal)
            error(['The value of ''ResponseType'' is invalid. ', ...
                'Expected ', 'number of classes in ''Y'' to be ', ...
                '2 for the ''dichotomous'' option.']);
        end
        if(~param.Ordinal && nClass == 2)
            param.Coding = 'onevsone'; 
        end
        if(nClass ~= 1) 
            t = templateLinear('Lambda', 0, 'Learner', 'logistic');
            [Mdl] = fitcecoc(cov, Yp, 'Coding', param.Coding, ...
                'Learners', t, 'FitPosterior', true);
            [~, ~, ~, posterior] =  predict(Mdl, cov);
            r = (Yp - 1) - sum(repmat(0:nClass-1, nSample, 1) .* posterior, 2);
        else
            r = zeros(size(Yp));
        end
    end
    Xmax = max(abs(min(min(X))), max(max(X)));
    Norm = Xmax * max(sum(r(r>0)), sum(abs(r(r<0))));
    if(param.BatchSize >= nVariant)
        C = ((X' * r) / Norm) .^ 2; %% SKAT Scores - m x 1 column vector
    else
        C = zeros(nVariant, 1);
        index = 0;
        while(index < nVariant)
            nItem = min(param.BatchSize, nVariant - index);
            indices = index + (1:nItem);
            C(indices) = ((double(X(:, indices)') * r) / Norm) .^ 2;
            index = index + nItem;
        end
    end
end

