function [ target ] = structconcat( target, source )
    p = inputParser;
    addRequired(p, 'target', @isstruct);
    addRequired(p, 'source', @isstruct);
    parse(p, target, source);
    f = fieldnames(source);
    for i = 1:length(f)
        target.(f{i}) = source.(f{i}); 
    end
end


