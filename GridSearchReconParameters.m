function [params,optval,info,varargout] = GridSearchReconParameters(params_range,reconfun,optfun,search_spec)
% Use iterative coarse-to-fine grid search algorithm to iterate over
% parameter space to find "best" reconstruction according to specified
% quality metric (e.g., NRMSE). The optfun() function takes the log-scale
% parameters passed to reconfun() as a second argument.
%
% Note: parameters are computed in the log scale so be sure reconfun() and
% optfun() do the appropriate conversion.
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

Nparams = size(params_range,1);

if ~exist('search_spec','var') || isempty(search_spec), search_spec = struct(); end
if ~isfield(search_spec,'verbose') || isempty(search_spec.verbose), search_spec.verbose = false; end
if ~isfield(search_spec,'ngrid') || isempty(search_spec.ngrid), search_spec.ngrid = 5; end % >= 2, >= 4 for niters > 1
if ~isfield(search_spec,'niters') || isempty(search_spec.ngrid), search_spec.niters = 1; end

params = [];
optval = Inf;
strparams = ['%g',repmat(', %g',1,Nparams-1)];

if any(search_spec.ngrid < 2), error('Does not support single-point grid.'); end
if search_spec.niters > 1 && any(search_spec.ngrid < 4), error('Grid must be >= 4 for multiple iterations.'); end

if isscalar(search_spec.ngrid), search_spec.ngrid = search_spec.ngrid.*ones(1,Nparams); end

% setup initial search space
info = struct('optvals',NaN([prod(search_spec.ngrid),search_spec.niters]),'params',{cell(Nparams,search_spec.niters)});
logparams_step = (log(params_range(:,2)) - log(params_range(:,1))).'./(search_spec.ngrid-1);
logparams_center = mean(log(params_range),2).';

Noutputs = max(1,abs(nargout)-3);
outputs = cell(1,Noutputs); 
outputs_test = cell(1,Noutputs);

for ii = 1:search_spec.niters
    % get search space
    logparams_search = arrayfun(@(c,n,s) c + ((0:n-1)-(n-1)/2).*s,logparams_center,search_spec.ngrid,logparams_step,'UniformOutput',false);
    info.params(:,ii) = logparams_search;
    [logparams_search{:}] = ndgrid(logparams_search{:});
    logparams_search = reshape(cat(Nparams+1,logparams_search{:}),[],Nparams);
    
    if ii > 1
        subsbest = cell(1,Nparams);
        if Nparams > 1
            [subsbest{:}] = ind2sub(search_spec.ngrid,ibest);
        else
            subsbest{1} = ibest;
        end
        subsin = cell(1,Nparams);
        subsin(mod(search_spec.ngrid,2) == 1) = cellfun(@(s) s+(-1:1),subsbest(mod(search_spec.ngrid,2) == 1),'UniformOutput',false);
        subsin(mod(search_spec.ngrid,2) ~= 1) = cellfun(@(s) s+[-1,1],subsbest(mod(search_spec.ngrid,2) ~= 1),'UniformOutput',false);
        [subsin{:}] = ndgrid(subsin{:});
        subsout = cell(1,Nparams);
        subsout(mod(search_spec.ngrid,2) == 1) = arrayfun(@(n) [1,(n+1)/2,n],search_spec.ngrid(mod(search_spec.ngrid,2) == 1),'UniformOutput',false);
        subsout(mod(search_spec.ngrid,2) ~= 1) = arrayfun(@(n) [1,n],search_spec.ngrid(mod(search_spec.ngrid,2) ~= 1),'UniformOutput',false);
        [subsout{:}] = ndgrid(subsout{:});
        subsvalid = arrayfun(@(ii) subsin{ii} >= 1 & subsin{ii} <= search_spec.ngrid(ii),1:Nparams,'UniformOutput',false);
        subsvalid = all(cat(Nparams+1,subsvalid{:}),Nparams+1);
        subsin = cellfun(@(s) s(subsvalid),subsin,'UniformOutput',false);
        subsout = cellfun(@(s) s(subsvalid),subsout,'UniformOutput',false);
        if Nparams > 1
            indsin = sub2ind(search_spec.ngrid,subsin{:});
            indsout = sub2ind(search_spec.ngrid,subsout{:});
        else
            indsin = subsin{1};
            indsout = subsout{1};
        end
        info.optvals(indsout,ii) = info.optvals(indsin,ii-1);
    end
    
    % compute reconstruction for each parameter choice and get best one
    ibest = []; optval_best = Inf;
    for ip = 1:size(logparams_search,1)
        if isfinite(info.optvals(ip,ii))
            optval_test = info.optvals(ip,ii);
        else
            [outputs_test{:}] = reconfun(logparams_search(ip,:));
            optval_test = optfun(outputs_test{1},logparams_search(ip,:));
            if isempty(params) || optval_test < optval
                params = logparams_search(ip,:);
                optval = optval_test;
                outputs = outputs_test;
                if search_spec.verbose
                    fprintf(1,['[%s] Found new optimum: (', strparams, ') => %g.\n'],mfilename,params,optval);
                end
            end
            info.optvals(ip,ii) = optval_test;
        end
        if isempty(ibest) || optval_test < optval_best
            ibest = ip;
            optval_best = optval;
        end
        drawnow;
    end
    if search_spec.verbose
        fprintf(1,['[%s] Finished level %d. Optimum: (', strparams, ') => %g.\n'],mfilename,ii,params,optval);
        drawnow;
    end
    
    % construct next level
    logparams_center = logparams_search(ibest,:);
    logparams_step = (2.*logparams_step)./(search_spec.ngrid-1);
end

varargout = outputs;

end
