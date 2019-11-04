function potpourri( X, Y, W, k, R, SNP_info, b, omega, maxMarginalSignificance, outputFileName)
    p = inputParser;
    p.CaseSensitive = false;
    validInfo = @(x) validateattributes(x, {'string'}, ...
        {'2d', 'cell'});
    validIntegerScalar = @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','nonempty','real','nonnan','positive','integer'});
    validFloat = @(x) validateattributes(x, {'numeric'}, ...
        {'nonempty','real','nonnan','nonnegative'});
    addRequired(p, 'SNP_info', validInfo);
    addRequired(p, 'b', validIntegerScalar);
    addRequired(p, 'maxMarginalSignificance', validIntegerScalar);
    addRequired(p, 'omega', validFloat);

    [indicators, Info] = spadis_logistic(double(X), Y, W, k, R, omega);
    group_ids = zeros(size(indicators,1),1);
    chr_info =  str2double(SNP_info(:,2));
    group_ind = 1;
    for chr = 1:max(chr_info)
        ind_chr = indicators(chr_info == chr);
        chr_size = size(ind_chr,1);
        group_ids_chr = zeros(chr_size,1);
        indices = find(ind_chr);
        for i = 1:size(indices,1)
            ind_chr(max(1,indices(i)-b):min(chr_size,indices(i)+b)) = 1;
            group_ids_chr(max(1,indices(i)-b):min(chr_size,indices(i)+b)) = group_ind;
            group_ind = group_ind + 1;
        end
        indicators(chr_info == chr) = ind_chr;
        group_ids(chr_info == chr) = group_ids_chr;
    end
    Y = logical(Y);
    epistasis_test_mex(transpose(X(Y,indicators)), transpose(X(~Y,indicators)), SNP_info(indicators,:), group_ids(indicators), R(indicators), cellstr(outputFileName), maxMarginalSignificance);
end











