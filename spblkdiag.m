function blk = spblkdiag(varargin)
% Construct sparse block diagonal matrix from sparse submatrices, avoiding
% memory issue trying to sparsify output of MATLAB's blkdiag().
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

sz1s = cellfun(@(X) size(X,1),varargin);
sz2s = cellfun(@(X) size(X,2),varargin);
off1s = cumsum([0,sz1s(1:end-1)]);
off2s = cumsum([0,sz2s(1:end-1)]);

% ensure inputs are sparse
sps = cellfun(@issparse,varargin);
if ~all(sps)
    varargin(~sps) = cellfun(@(X) sparse(X),varargin(~sps),'UniformOutput',false);
end

% find indices and values of nonzero elements
[ind1s,ind2s,vals] = cellfun(@(X) find(X),varargin,'UniformOutput',false);

% offset indices into block matrix
ind1s = arrayfun(@(ii) ind1s{ii}(:) + off1s(ii),1:length(varargin),'UniformOutput',false);
ind1s = cat(1,ind1s{:});
ind2s = arrayfun(@(ii) ind2s{ii}(:) + off2s(ii),1:length(varargin),'UniformOutput',false);
ind2s = cat(1,ind2s{:});
vals = cellfun(@(val) val(:),vals,'UniformOutput',false);
vals = cat(1,vals{:});

% construct sparse block matrix
blk = sparse(ind1s,ind2s,vals,sum(sz1s),sum(sz2s));

end
