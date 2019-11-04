function [selected_out, info] = do_binary_search(func_run, func_dist, min, max, varargin)
    if(nargin < 4); error('Not enough input arguments.'); end
    p = inputParser;
    p.CaseSensitive = false;
    check_nn_numeric = @(x) (isnumeric(x) && x>=0);
	addParameter(p, 'Partition', 10000, @isnumeric);
    addParameter(p, 'Verbose', false, @islogical);
    addParameter(p, 'Direction', 'increasing', @ischar);
    addParameter(p, 'Tolerance', 0, check_nn_numeric);
    parse(p, varargin{:});
	opt = p.Results;
    switch(opt.Direction)
        case 'increasing'
            opt.DirectionValue = 1;
        case 'decreasing'
            opt.DirectionValue = -1;
        otherwise
            error(['Direction must be either increasing or decreasing.']);
    end
    info = struct();
    info.Partition = p.Results.Partition;
    info.Direction = p.Results.Direction;
    info.Tolerance = p.Results.Tolerance;
    opt.pmin = 0;
    opt.pmax = 100;
	opt.calc_value = @(pv) min * ((max / min) ^ (pv / 100));
    qr.xmin = 0;
    qr.xmax = opt.Partition - 1;
    xmid = calcmid(qr);
	[out, param, info.dist, info.noi] = run(func_run, func_dist, xmid, opt, 0);
    selected_out = out;
    info.param = param;
    [stop, qr] = adjust(qr, info.dist, opt);
	if(stop);  return; end
    while(true)
       xmid = calcmid(qr);
       range = qr.xmax - qr.xmin;
       if(range <= 1)
            if(qr.xmin == 0)
				[out, param, dist, info.noi] = run(func_run, func_dist, qr.xmin, opt, info.noi);
                if(abs(info.dist) > abs(dist)); selected_out = out; info.dist = dist; info.param = param; end
            end
            if(qr.xmax == opt.Partition - 1)
                [out, param, dist, info.noi] = run(func_run, func_dist, qr.xmax, opt, info.noi);
            end
			if(abs(info.dist) > abs(dist)); selected_out = out; info.dist = dist; info.param = param; end
			return;
       end
       [out, param, dist, info.noi] = run(func_run, func_dist, xmid, opt, info.noi);
	   [stop, qr] = adjust(qr, dist, opt);
	   if(stop); selected_out = out; info.dist = dist; info.param = param; return; end
	   if(abs(info.dist) > abs(dist)); selected_out = out; info.dist = dist; info.param = param; end
    end
end

function [out, param, dist, noi] = run(func_run, func_dir, x, opt, noi)
	pv = opt.pmin + (opt.pmax - opt.pmin) * (1 + x) / opt.Partition;
	param = opt.calc_value(pv);
	[out] = func_run(param);
	[dist] = func_dir(param, out);
	noi = noi + 1;
    if(opt.Verbose)
		disp(['Iteration ', num2str(noi, 4), ' - ', 'percentage : ', num2str(pv, 4), '%', ', param : ', num2str(param, 4), ', dist : ', num2str(dist, 4)]);
    end
end

function [stop, qr] = adjust(qr, dist, opt)
    stop = false;
    if(abs(dist) <= opt.Tolerance)
        stop = true;
        return;
    end 
    if((dist * opt.DirectionValue) > 0)  % Decrease d
        qr.xmax = calcmid(qr);
    else
        qr.xmin = calcmid(qr);
    end
end

function [mid] = calcmid(qr)
    mid = ((qr.xmin + qr.xmax) / 2);
end



