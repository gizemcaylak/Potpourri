function [ I, Info ] = spadis( C, W, k, varargin )
    p = inputParser;
    p.CaseSensitive = false;
    validScoring = @(x) validateattributes(x, {'numeric'}, ...
        {'column','nonempty','nonsparse','real','nonnan','nonnegative'});
    validNetwork = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'square', 'nonempty', 'real', 'nonnan'});
    validIntegerScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','nonempty','real','nonnan','positive','integer'});
    validDelta = @(x) validateattributes(x, {'numeric'}, ...
        {'nonempty','nonsparse','vector','real','nonnan','nonnegative'});
    validBeta = @(x) validateattributes(x, {'numeric'}, ...
        {'2d','nonempty','nonsparse','real','nonnan','nonnegative'});
    addRequired(p, 'C', validScoring);
    addRequired(p, 'W', validNetwork);
    addRequired(p, 'k', validIntegerScalar);
    addParameter(p, 'Delta', [], validDelta);
    addParameter(p, 'Beta', [], validBeta);
    addParameter(p, 'NumDelta', 10, validIntegerScalar);
    addParameter(p, 'NumBeta', 20, validIntegerScalar);    
    addParameter(p, 'MaxIter', 10, validIntegerScalar);
    parse(p, C, W, k, varargin{:});
    param = p.Results;
    nVariant = size(C, 1);
    invalidateMismatch(C, W, 'C', 'W', 'row');
    assert(k < nVariant, ['Expected cardinality constaint ''k'' to be', ...
        ' lesser than the number of rows in ''C''.']);
    if(~issparse(W)); W = sparse(W); end
    checkUsingDefaults = @(p,varname) any(strcmp(p.UsingDefaults,varname));
    Info = struct();
    if(checkUsingDefaults(p, 'Delta'))
        [DeltaMax, DeltaMin] = spadis_drange(C, W, k, 'Partition', 2^(param.MaxIter));
        param.Delta = logspace(log10(DeltaMin), log10(DeltaMax), param.NumDelta + 1);
        param.Delta = param.Delta(2:end);
        Info.DeltaMin = DeltaMin;
        Info.DeltaMax = DeltaMax;
    end
    nDelta = length(param.Delta);
    if(checkUsingDefaults(p, 'Beta'))
        param.Beta = zeros(nDelta, param.NumBeta);
        for iDelta = 1:nDelta
            D = param.Delta(iDelta);
            [BetaMin, BetaMax] = spadis_mex(C, W, k, D, 'BetaRange');
            BetaMin = BetaMin * 0.999;
            BetaMax = max(BetaMin, BetaMax) * 1.001;
            betaValues = logspace(log10(BetaMin), log10(BetaMax), param.NumBeta);
            param.Beta(iDelta, :) = [betaValues];
%             param.Beta(iDelta, :) = [0, betaValues, Inf];
        end
        Info.BetaMin = param.Beta(:, 1);
        Info.BetaMax = param.Beta(:, end);
    else % not using default beta values
        if(size(param.Beta, 1) > size(param.Beta, 2))
            param.Beta = param.Beta';
        end
        param.Beta = repmat(param.Beta, nDelta, 1);
    end
    nBeta = size(param.Beta, 2);
    I = zeros(nVariant, nDelta, nBeta, 'logical');
    Info.FunVal = zeros(nDelta, nBeta);
    for iDelta = 1:nDelta
        D = param.Delta(iDelta);
        Beta = param.Beta(iDelta, :);
        [indicators, funVal] = spadis_mex(C, W, k, D, Beta);
        I(:,iDelta,:) = indicators;
        Info.FunVal(iDelta, :) = funVal;
    end
    Info.Delta = param.Delta';
    Info.Beta = param.Beta;
    I1 = I(:, 1:end-1, :);
    I2 = I(:, 2:end, :);
    JaccardCof = sum(I1 & I2, 1) ./ sum(I1 | I2, 1);
    Info.JCofDelta = reshape(JaccardCof, nDelta - 1, nBeta);
    I1 = I(:, :, 1:end-1);
    I2 = I(:, :, 2:end);
    JaccardCof = sum(I1 & I2, 1) ./ sum(I1 | I2, 1);
    Info.JCofBeta = reshape(JaccardCof, nDelta, nBeta - 1);
    [~, si] = sort(C, 'descend');
    Info.IndicatorsBetaMin = zeros(nVariant, 1, 'logical');
    Info.IndicatorsBetaMin(si(1:k)) = true;
    I2 = repmat(Info.IndicatorsBetaMin, 1, nDelta, nBeta);
    JaccardCof = sum(I & I2, 1) ./ sum(I | I2, 1);
    Info.JCofBetaMin = reshape(JaccardCof, nDelta, nBeta);
    I2 = repmat(I(:, :, end), 1, 1, nBeta);
    JaccardCof = sum(I & I2, 1) ./ sum(I | I2, 1);
    Info.JCofBetaMax = reshape(JaccardCof, nDelta, nBeta);
end











