function [] = func_spadis(NSF)
    combinedPath = [ 'autism_data/'];
    resultPath = [ 'autism_results/'];
    load([combinedPath, 'combined.mat']);
    load([combinedPath, 'evaluation_sets.mat']);

    nFold = 10;
    nPC = 4;

    Y = double(combined.subjects.AFFECTION == 2);
    Ytest = Y(test_set);
    Y = Y(trval_set);
    disp('[Running] Loading variants...');
    load([combinedPath, 'variants.mat']);
    disp('[Done] Loading variants.');
    Xtest = variants(test_set, :);
    variants = variants(trval_set, :);
    load([combinedPath, 'skat_scores_all.mat']);

    data.X = variants;
    clear variants
    disp('[Running] Running SPADIS...');
    [indicators, info] = spadis_all(data, tr_sets, NSF);
    disp('[Done] Running SPADIS...');
    disp('[Running] Saving results of SPADIS...');
    save([resultPath, strcat('spadis_', num2str(NSF),'_results.mat')], 'indicators', 'info');
    disp('[Done] Saving results of SPADIS...');
    
    disp('[Running] Evaluating results of SPADIS...');
    [Mdl, fitInfo] = fitclinear(double(data.X(:, indicators)) ...
        ,double(data.Y), 'ObservationsIn','rows' ...
        ,'Learner','logistic', 'Regularization','ridge');
    [Yhat] = predict(Mdl, double(Xtest(:, indicators)));
    stats = struct();
    [table] = confusionmat(Ytest, Yhat); % Actual vs Predicted
    stats.precision = table(2,2) / sum(table(:, 2));
    stats.recall = table(2,2) / sum(table(2, :));
    stats.accuracy = (table(1,1) + table(2,2)) / sum(sum(table));
    stats.fscore = 2 / ((1/stats.precision) + (1/stats.recall));
    stats.mcc = findMCC(table);
    disp('[Done] Evaluating results of SPADIS.');
    disp('[Running] Saving evaluation results of SPADIS...');
    save([resultPath, strcat('spadis_eval_',num2str(NSF),'_results.mat')], 'table', 'stats', 'Mdl', 'fitInfo', 'Ytest', 'Yhat');
    disp('[Done] Saving evaluation results of SPADIS...');
    
    

    S = struct();
    
    load([resultPath, strcat('spadis_eval_',num2str(NSF),'_results.mat')]);
    disp('[Running] Evaluating results of SPADIS...');
    [Mdl, fitInfo] = fitclinear(double(data.X(:, indicators)) ...
        ,double(data.Y), 'ObservationsIn','rows' ...
        ,'Learner','logistic', 'Regularization','ridge');
    [Yhat] = predict(Mdl, double(Xtest(:, indicators)));
    stats = struct();
    [table] = confusionmat(Ytest, Yhat); % Actual vs Predicted
    stats.precision = table(2,2) / sum(table(:, 2));
    stats.recall = table(2,2) / sum(table(2, :));
    stats.accuracy = (table(1,1) + table(2,2)) / sum(sum(table));
    stats.fscore = 2 / ((1/stats.precision) + (1/stats.recall));
    stats.mcc = findMCC(table);
    
    S.SPADIS = struct();
    S.SPADIS.stats = stats;
    S.SPADIS.indicators = indicators;
    S.SPADIS.table = table;
    disp('[Done] Evaluating results of SPADIS.');
    
    disp('[Running] Saving evaluation results of baseline methods...');
    save([resultPath, strcat('baseline_eval_', num2str(NSF),'_results.mat')], 'S');
    disp('[Done] Saving evaluation results of baseline methods...');
    
    
end














