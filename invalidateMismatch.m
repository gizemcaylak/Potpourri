function [ ] = invalidateMismatch( var1, var2, varname1,  ...
                                    varname2, match_type, varargin)
    p = inputParser;
    validString = @(x) validateattributes(x, {'char', 'string'}, ...
        {'scalartext', 'nonempty'});
    validMatchType = @(x) any(validatestring(x, ...
        {'row', 'column', 'dimension', 'length'}));
    addRequired(p, 'var1', @(x) true);
    addRequired(p, 'var2', @(x) true);
    addRequired(p, 'varname1', validString);
    addRequired(p, 'varname2', validString);
    addRequired(p, 'match_type', validMatchType);
    parse(p, var1, var2, varname1, varname2, match_type, varargin{:});
    checkRow = false;
    checkColumn = false;
    checkLength = false;
    checkDimension = false;
    switch(p.Results.match_type)
        case 'row'
            checkRow = true;
        case 'column'
            checkColumn = true;
        case 'dimension'
%             checkRow = true;
%             checkColumn = true;
            checkDimension = true;
        case 'length'
            checkLength = true;
    end
    dimensionMismatchId = 'MATLAB:invalidateMismatch:dimensionMismatch';
    if(checkRow && size(var1, 1) ~= size(var2, 1))
        errorMessage = ['Expected number of rows in ''', ...
            varname1, ''' and ''', varname2, ''' to be equal.'];
        throwAsCaller(MException(dimensionMismatchId, errorMessage));
    end
    if(checkColumn && size(var1, 2) ~= size(var2, 2))
        errorMessage = ['Expected number of columns in ''', ...
            varname1, ''' and ''', varname2, ''' to be equal.'];
        throwAsCaller(MException(dimensionMismatchId, errorMessage));
    end
    if(checkLength && length(var1) ~= length(var2))
        errorMessage = ['Expected lengths of ''', ...
            varname1, ''' and ''', varname2, ''' to be equal.'];
        throwAsCaller(MException(dimensionMismatchId, errorMessage));
    end
    if(checkDimension && ~isequal(size(var1), size(var2)))
        errorMessage = ['Expected dimensions of ''', ...
            varname1, ''' and ''', varname2, ''' to be equal.'];
        throwAsCaller(MException(dimensionMismatchId, errorMessage));
    end
    
end

