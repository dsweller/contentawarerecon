function xout = my_zpad_crop(x,sz)
% One-stop function for padding or cropping an image or numeric/logical signal.
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

sz_x = size(x); sz_x(end+1:length(sz)) = 1;
sz(end+1:length(sz_x)) = sz_x(length(sz)+1:end);

idxs_out = arrayfun(@(minlen) 1:minlen,min(sz,sz_x),'UniformOutput',false);
idxs_in = idxs_out;

idxs_center = arrayfun(@(minlen,maxlen) floor(maxlen/2)-floor(minlen/2)+(1:minlen),min(sz,sz_x),max(sz,sz_x),'UniformOutput',false);
idxs_out(sz_x < sz) = idxs_center(sz_x < sz); % output is centered for zero-padding
idxs_in(sz_x > sz) = idxs_center(sz_x > sz); % input is centered for cropping

if isnumeric(x)
    if issparse(x)
        xout = sparse([],[],[],sz(1),sz(2),nnz(x));
    else
        xout = zeros(sz,'like',x);
    end
elseif islogical(x)
    xout = false(sz);
elseif iscell(x)
    xout = cell(sz);
else
    error('Unsupported type');
end
xout(idxs_out{:}) = x(idxs_in{:});

end
