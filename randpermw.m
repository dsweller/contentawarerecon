function out = randpermw(weights,K)
% Compute weighted random permutation of K values from N things (where N is
% the length of weights). This extends randperm() to do variable density
% sampling. But, this code is kind of slow.
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

cweights = cumsum(weights(:));
out = rand(1,K);
if K > length(cweights), error('K is too large relative to number of weights.'); end

used = false(size(cweights));
for ii = 1:K
    outval = find(~used & (out(ii).*cweights(end) <= cweights),1);
    if isempty(outval), outval = find(~used,1,'last'); end
    out(ii) = outval;
    used(out(ii)) = true;
    cweights(out(ii):end) = cweights(out(ii):end) - weights(out(ii));
end

out = sort(out);

end
