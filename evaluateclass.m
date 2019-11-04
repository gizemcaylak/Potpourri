function [ P, stats ] = evaluateclass( G, Ghat, varargin )
    p = inputParser;
    validGroup = @(x) validateattributes(x, {'logical', 'numeric'}, ...
        {'nonempty'});
    allStats = {'Accuracy', 'Prevalence', 'Precision', 'NPV', 'Recall', ...
         'Specificity', 'F1score', 'MCC', 'Informedness', 'Markedness'}';
    validCriterion = @(x) any(validatestring(x, allStats));
    addRequired(p, 'G', validGroup);
    addRequired(p, 'Ghat', validGroup);
    addOptional(p, 'Criterion', 'mcc', validCriterion);
    parse(p, G, Ghat);
    invalidateMismatch(G, Ghat, 'G', 'Ghat', 'row');
    invalidateMismatch(G, Ghat, 'G', 'Ghat', 'dimension');
    nObs = size(G, 1);
    dim = size(G);
    if(numel(dim) == 2); dim = [dim, 1]; end
    if(numel(dim) == 1); dim = [dim, 1, 1]; end
    criterionIndex = strcmpi(allStats, p.Results.Criterion);
    P = zeros(dim(2:end));
    for i = (numel(G) / nObs):-1:1
        ind = (1:nObs) + (i - 1) * nObs;
        [conf] = confusionmat(G(ind), Ghat(ind));
        if(nargout >= 2)
            S = computeStats(conf, allStats);
            stats(i) = cell2struct(S, allStats);
            P(i) = S{criterionIndex};
        else
            S = computeStats(conf, {p.Results.Criterion});
            P(i) = S{1}; 
        end
    end
end

function [ S ] = computeStats(conf, statList)
    S = cell(size(statList));
    TP = conf(2, 2);
    TN = conf(1, 1);
    FP = conf(1, 2);
    FN = conf(2, 1);
    nObs = sum(sum(conf));
    for i = 1:length(statList)
        switch(lower(statList{i}))
            case 'accuracy'
                S{i} = (TP + TN) / nObs;
            case 'prevalence'
                S{i} = (TP + FN) / nObs;
            case 'precision'
                S{i} = TP / (TP + FP);
            case 'npv'
                S{i} = TN / (TN + FN);
            case 'recall'
                S{i} = TP / (TP + FN);
            case 'specificity'
                S{i} = TN / (TN + FP);
            case 'f1score'
                S{i} = (2 * TP) / (2 * TP + FP + FN);
            case 'mcc'
                S{i} = (TP*TN - FP*FN) / sqrt((TP + FP) * (TP + FN) * ...
                                              (TN + FP) * (TN + FN));
            case 'informedness'
                S{i} = (TP / (TP + FN)) + (TN / (TN + FP)) - 1;
            case 'markedness'
                S{i} = (TP / (TP + FP)) + (TN / (TN + FN)) - 1;
            otherwise
                error(['Invalid statistic requested:', statList{i}]);
        end
    end
end

function value = getfieldi(S,field)
    names   = fieldnames(S);
    isField = strcmpi(field,names);  
    if any(isField)
      value = S.(names{isField});
    else
      value = [];
    end
end
