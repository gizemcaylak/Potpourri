function [ I, Info] = spadis_logistic( X, Y, W, k, R, omega, varargin )
    p = inputParser;
    p.CaseSensitive = false;
    validMatrix = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'nonempty', 'real', 'nonsparse', 'nonnan', 'finite'});
    validPhenotype = @(x) validateattributes(x, ...
        {'numeric','logical', 'char'}, ...
        {'vector', 'nonempty', 'real', 'nonnan', 'nonsparse'});
    validNetwork = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'2d', 'square', 'nonempty', 'real', 'nonnan'});
    validIntegerScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','nonempty','real','nonnan','positive','integer'});
    validFloat = @(x) validateattributes(x, {'numeric'}, ...
        {'nonempty','real','nonnan','nonnegative'});
    validVector = @(x) validateattributes(x, {'numeric'}, ...
        {'nonempty','nonsparse','vector','real','nonnan','nonnegative'});
    validScoring = @(x) any(validatestring(x, {'skat', 'chi2'}));
    validPCs = @(x) validateattributes(x, {'numeric'}, ...
        {'2d', 'real', 'nonsparse', 'nonnan', 'finite'});
    validNumPCs = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','nonempty','real','nonnan','nonnegative','integer'});
    validCoding = @(x) any(validatestring(x, ...
        {'allpairs', 'onevsone', 'binarycomplete', 'denserandom', ...
        'onevsall', 'ordinal', 'sparserandom', 'ternarycomplete'}));
    validCriterion = @(x) any(validatestring(x, ...
        {'accuracy', 'f1score', 'mcc', 'precision', 'recall'}));
    validVerbose = @(x) validateattributes(x, {'numeric', 'logical'}, ...
        {'scalar','nonempty','real','nonnan','nonnegative','integer'});
    validConditionPositive = @(x) validateattributes(x, ...
        {'numeric','logical', 'char'}, ...
        {'scalar', 'nonempty', 'real', 'nonnan', 'nonsparse'});
    addRequired(p, 'X', validMatrix);
    addRequired(p, 'Y', validPhenotype);
    addRequired(p, 'W', validNetwork);
    addRequired(p, 'k', validIntegerScalar);
    addRequired(p, 'R', validVector);
    addRequired(p, 'omega', validFloat);
    addParameter(p, 'Delta', [], validVector);
    addParameter(p, 'Beta', [], validVector);
    addParameter(p, 'NumDelta', 10, validIntegerScalar);
    addParameter(p, 'NumBeta', 20, validIntegerScalar);   
    addParameter(p, 'MaxIter', 10, validIntegerScalar);
    addParameter(p, 'Scoring', 'skat', validScoring);
    addParameter(p, 'PCs', NaN, validPCs);
    addParameter(p, 'NumPCs', 0, validNumPCs);
    addParameter(p, 'BatchSize', 10000, validIntegerScalar);
    addParameter(p, 'Coding', 'onevsone', validCoding);
    addParameter(p, 'CVPartition', [], @(x) isa(c, 'cvpartition'));
    addParameter(p, 'KFold', 10, validIntegerScalar);
    addParameter(p, 'Criterion', 'mcc', validCriterion);
    addParameter(p, 'Verbose', 0, validVerbose);
    addParameter(p, 'ConditionPositive', 0, validConditionPositive);
    parse(p, X, Y, W, k, R, omega, varargin{:});
    param = p.Results;
    [nSample, nVariant] = size(X);
    invalidateMismatch(X, Y, 'X', 'Y', 'row');
    invalidateMismatch(X, W, 'X', 'W', 'column');
    assert(k < nVariant, ['Expected cardinality constaint ''k'' to be', ...
        ' lesser than the number of columns in ''X''.']);
    reporter = timereporter(p.Results.Verbose);
    obj1 = reporter.printRunning('SPADIS_Logistic', 1);
    Info = struct();
    [ClassNames, ~, Yind] = unique(Y);
    if(length(ClassNames) > 2)
        error('Only dichotomous phenotypes are supported.');
    end
    checkUsingDefaults = @(p,varname) any(strcmp(p.UsingDefaults,varname));
    if(checkUsingDefaults(p, 'ConditionPositive') && ...
            (isstring(Y) || ~ismember(0, ClassNames)))
        error(['The value of condition positive is ambiguous. Please ', ...
            'specify it using the ''ConditionPositive'' argument.']);
    end
    [condPosFound, condPosIndex] = ismember(...
        param.ConditionPositive, ClassNames);
    if(~condPosFound)
        error('The value of condition positive is invalid.');
    end
    Yp = logical((Yind == condPosIndex) - 1);
    SKAToptions = structsubset(param, {'NumPCs', 'BatchSize'});
    SKAToptions.ResponseVarType = 'dichotomous';
    if(~checkUsingDefaults(p, 'PCs')); SKAToptions.PCs = param.PCs; end
    [C, Info.PCs] = computeSKAT(X, Yp, SKAToptions);
    C = C + C.*(omega*R);
    SPADISfields = {'NumDelta', 'NumBeta', 'MaxIter'};
    SPADISoptions = structsubset(p, param, {'Delta', 'Beta'});
    SPADISoptions = structsubset(param, SPADISfields, SPADISoptions);
    [I, SPADISinfo] = spadis(C, W, k, SPADISoptions);
    [Info] = structconcat(Info, SPADISinfo);
    NumDelta = size(Info.Delta, 1);
    NumBeta = size(Info.Beta, 2);
    if(checkUsingDefaults(p, 'CVPartition'))
        param.CVPartition = cvpartition(nSample, 'KFold', param.KFold);
    end
    Info.CVPartition = param.CVPartition;
    SPADISoptions = structsubset(Info, {'Delta', 'Beta'});
    nFold = param.CVPartition.NumTestSets;
    Yhat = zeros(nSample, NumDelta, NumBeta, 'logical');
    indicatorsList = {};
    skatAll = C;

    for iFold = 1:nFold
        obj = reporter.printRunning(['Cross-validation Set ', ...
            num2str(iFold)], 1, NumDelta * NumBeta, 2);
        tr_set = training(param.CVPartition, iFold);
        te_set = test(param.CVPartition, iFold);
        [C] = computeSKAT(X(tr_set, :), Yp(tr_set), SKAToptions);
        C = C + C.*(omega*R);
        [indicators] = spadis(C, W, k, SPADISoptions);
        indicatorsList{iFold} = indicators;
        for iDelta = 1:NumDelta
            for iBeta = 1:NumBeta
                Xselected = double(X(:, indicators(:, iDelta, iBeta)));
                Xtrain = Xselected(tr_set, :);
                Xtest = Xselected(te_set, :);
                Mdl = fitclinear(Xtrain, Yp(tr_set), 'Learner', 'logistic','Regularization', 'ridge');
                [labels] = predict(Mdl, Xtest);
                Yhat(te_set, iDelta, iBeta) = labels;
                obj = obj.printProgress();
            end
        end
        obj.printDone();
    end
    Info.ClassPerf = evaluateclass(repmat(Yp, 1, NumDelta, NumBeta), ...
                                    Yhat, param.Criterion);
    [ClassPerf, deltaInd] = max(Info.ClassPerf, [], 1);
    [Info.MaxClassPerf, Info.BetaIndex] = max(ClassPerf, [], 2);
    Info.DeltaIndex = deltaInd(Info.BetaIndex);
    [~, Info.StatsMaxPerf] = evaluateclass(Yp, ...
        Yhat(:, Info.DeltaIndex, Info.BetaIndex));
    Info.DeltaSelected = Info.Delta(Info.DeltaIndex);
    Info.BetaSelected = Info.Beta(Info.DeltaIndex, Info.BetaIndex);
    I = I(:, Info.DeltaIndex, Info.BetaIndex);
    Info.CardinalityConstraint = k;
    Info.ClassNames = ClassNames;
    Info.NumDelta = NumDelta;
    Info.NumBeta = NumBeta;
    Info.Criterion = param.Criterion;
    obj1.printDone();
end













