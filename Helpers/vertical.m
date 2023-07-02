function V = vertical(V);

if length(size(V))>2 | size(V,1)*size(V,2)>numel(V) error('Input is not a vector!'); end
if size(V,2) > 1 V = V'; end 