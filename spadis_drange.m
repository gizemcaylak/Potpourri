function [ Dmax, Dmin ] = spadis_drange(C, W, k, varargin )
    if(nargin < 3); error('Not enough input arguments.'); end
    if(nnz(W) == 0); error('W cannot be empty!'); end
    p = inputParser;
    f_pos_numeric = @(x) (isnumeric(x) && x>0);
    addParameter(p, 'Partition', 1000, f_pos_numeric);
    addParameter(p, 'Verbose', false, @islogical);
    parse(p, varargin{:});
%     [Dmin] = spadis_mexa(C, W, k, 'DeltaRange');
    Dmin = full(double(min(min(W(W > 0)))));
    Dmax = full(double(sum(sum(W))));
    Beta = Inf;
    func_run = @(param) run_spadis(C, W, k, param, Beta);
    func_dist = @(param, F) ((F == k) / param) + (F ~= k) * (-1 / Dmin);
    optB.Verbose = p.Results.Verbose;
    optB.Direction = 'decreasing';
    optB.Partition = p.Results.Partition;
    optB.Tolerance = 0;
    [~, info] = do_binary_search(func_run, func_dist, Dmin, Dmax, optB);
    Dmax = info.param;
end

function [F] = run_spadis(C, W, k, D, Beta)
    [~, F] = spadis_mex(C, W, k, D, Beta);
    F = imag(F);
end
