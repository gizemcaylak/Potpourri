function [ target ] = structsubset( varargin )
    p = inputParser;
    parserEnabled = nargin > 1 && isa(varargin{1}, 'inputParser');
    if(parserEnabled)
        addRequired(p, 'Parser', @(x) isa(x, 'inputParser'));
    end
    addRequired(p, 'Source', @isstruct);
    addRequired(p, 'FieldNames', @iscell);
    addOptional(p, 'Target', struct(), @isstruct);
    parse(p, varargin{:});
    checkUsingDefaults = @(p,varname) any(strcmp(p.UsingDefaults,varname));
    source = p.Results.Source;
    fieldnames = p.Results.FieldNames;
    target = p.Results.Target;
    D = cellfun(@(f) getfield(source,f),fieldnames, 'UniformOutput',false);
    if(parserEnabled)
        mask = arrayfun(@(f) ~checkUsingDefaults(p.Results.Parser, f), ...
            fieldnames, 'UniformOutput', true);
    else
        mask = ones(size(D), 'logical');
    end
    target = structconcat(target, cell2struct(D(mask), fieldnames(mask), 2));
end

