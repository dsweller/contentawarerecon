function [params,optval,info,varargout] = OptimizeReconParameters(initialsimplex,lbs,ubs,reconfun,optfun,opt_spec)
% Use Nelder-Mead algorithm to iterate over parameter space to find "best"
% reconstruction according to specified quality metric (e.g., NRMSE). This
% function has similar capability to MATLAB's fminsearch(), but isn't
% encumbered by MATLAB's inexplicable inability to easily specify an
% initial simplex. Default parameters adapt based on number of parameters
% (see F Gao and L Han, Comput. Optim. Appl. 51(1), pp. 259-277, 2010).
%
% Note: parameters are assumed to extend from -Infty to Infty (if not
% bounded). For regularization parameters that are usually nonnegative, it
% may be advisable to optimize over their log() values instead. Ensure
% reconfun() does the appropriate conversions.
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

if isscalar(initialsimplex)
    Nparams = initialsimplex;
    initialsimplex = [];
else
    Nparams = size(initialsimplex,1);
end

if ~exist('opt_spec','var') || isempty(opt_spec), opt_spec = struct(); end
if ~isfield(opt_spec,'niters') || isempty(opt_spec.niters), opt_spec.niters = 10*Nparams; end
if ~isfield(opt_spec,'verbose') || isempty(opt_spec.verbose), opt_spec.verbose = false; end
if ~isfield(opt_spec,'tolfun') || isempty(opt_spec.tolfun), opt_spec.tolfun = 1e-4; end
if ~isfield(opt_spec,'tolparams') || isempty(opt_spec.tolparams), opt_spec.tolparams = 1e-4; end
if ~isfield(opt_spec,'r') || isempty(opt_spec.r), opt_spec.r = 1; end % 0 < r <= 1
if ~isfield(opt_spec,'gamma') || isempty(opt_spec.gamma), opt_spec.gamma = min(2,1+2/Nparams); end % gamma > 1
if ~isfield(opt_spec,'beta') || isempty(opt_spec.beta), opt_spec.beta = max(0.5,0.75-0.5/Nparams); end % 0 < beta < 1
if ~isfield(opt_spec,'delta') || isempty(opt_spec.delta), opt_spec.delta = max(0.5,1-1/Nparams); end % 0 < delta < 1

if isempty(lbs), lbs = -Inf(Nparams,1); end
if isempty(ubs), ubs = Inf(Nparams,1); end

bounded = isfinite(lbs) & isfinite(ubs);
bounded_below = isfinite(lbs) & ~isfinite(ubs);
bounded_above = isfinite(ubs) & ~isfinite(lbs);
unbounded = ~isfinite(lbs) & ~isfinite(ubs);

% initialization
simplex = initialsimplex;
while isempty(simplex) || isdegenerate(simplex)
    % randomize initial simplex
    simplex = rand(Nparams,Nparams+1);
    if any(bounded) % uniform distribution
        simplex(bounded,:) = bsxfun(@plus,bsxfun(@times,ubs(bounded)-lbs(bounded),simplex(bounded,:)),lbs(bounded));
    end
    if any(bounded_below) % exponential distribution
        simplex(bounded_below,:) = bsxfun(@minus,lbs(bounded_below),log(simplex(bounded_below,:)));
    end
    if any(bounded_above) % reversed exponential distribution
        simplex(bounded_above,:) = bsxfun(@plus,ubs(bounded_above),log(simplex(bounded_above,:)));
    end
    if any(unbounded) % Laplacian (two-sided exponential) distribution
        simplex(unbounded,:) = (-1).^(simplex(unbounded,:) < 0.5).*log(abs(2.*simplex(unbounded,:)-1));
    end
end

t0 = tic();

% evaluate initial simplex points and keep output for best initial params
funcCount = 0;
Noutputs = abs(nargout)-3;
outputs = cell(1,Noutputs); outputs_test = cell(1,Noutputs);
[outputs{:}] = reconfun(simplex(:,1)); funcCount = funcCount + 1;
optval = optfun(outputs{1});
optvals = [optval,NaN(1,Nparams)];
params = simplex(:,1);
for n = 2:Nparams+1
    [outputs_test{:}] = reconfun(simplex(:,n)); funcCount = funcCount + 1;
    optvals(n) = optfun(outputs_test{1});
    if optvals(n) < optval
        % new minimum
        optval = optvals(n);
        outputs = outputs_test;
        params = simplex(:,n);
    end
end

% sort initial simplex
[optvals,indssort] = sort(optvals);
simplex = simplex(:,indssort);

% storage for tracking progress
flats = NaN(opt_spec.niters,1);
dims = NaN(opt_spec.niters,1);
optvalhist = [optval;NaN(opt_spec.niters,1)];
steps = char(zeros(opt_spec.niters,1));

if opt_spec.verbose, fprintf(1,'[%s] Initial simplex (optval = %g)\n',mfilename,optval); drawnow; end

% Nelder-Mead algorithm
for ii = 1:opt_spec.niters
    % check for convergence (flat and small)
    flats(ii) = optvals(end) - optvals(1);
    dims(ii) = getmaxdimension(simplex,lbs,ubs);
    if flats(ii) <= opt_spec.tolfun*max(abs(optvals)) && dims(ii) <= opt_spec.tolparams
        break;
    end

    % midpoint on best segment of simplex
    midpoint = mean(simplex(:,1:end-1),2);
    
    % reflect worst point about midpoint
    reflection = (1+opt_spec.r).*midpoint - opt_spec.r.*simplex(:,end);
    reflection = project_bounds(reflection,lbs,ubs);
    if isdegenerate([simplex(:,1:end-1),reflection])
        optval_r = [];
    else
        [outputs_test{:}] = reconfun(reflection); funcCount = funcCount + 1;
        optval_r = optfun(outputs_test{1});
    end

    % test reflection
    if ~isempty(optval_r) && optval_r < optvals(1) 
        % new minimum
        optval = optval_r;
        outputs = outputs_test;
        params = reflection;
        
        % expand along reflection to try to get even better point
        expansion = opt_spec.gamma.*reflection - (opt_spec.gamma-1).*midpoint;
        expansion = project_bounds(expansion,lbs,ubs);
        if isdegenerate([simplex(:,1:end-1),expansion])
            optval_e = [];
        else
            [outputs_test{:}] = reconfun(expansion); funcCount = funcCount + 1;
            optval_e = optfun(outputs_test{1});
        end
        
        if ~isempty(optval_e) && optval_e < optval_r
            optvals = [optval_e,optvals(1:end-1)];
            simplex = [expansion,simplex(:,1:end-1)];
            
            % new minimum
            optval = optval_e;
            outputs = outputs_test;
            params = expansion;
            
            steps(ii) = 'E';
            if opt_spec.verbose, fprintf(1,'[%s] Iter %d: Expand (optval = %g)\n',mfilename,ii,optval); drawnow; end
        else % just use reflected point
            optvals = [optval_r,optvals(1:end-1)];
            simplex = [reflection,simplex(:,1:end-1)];
            
            steps(ii) = 'R';
            if opt_spec.verbose, fprintf(1,'[%s] Iter %d: Reflect (optval = %g)\n',mfilename,ii,optval); drawnow; end
        end
        optvalhist(ii+1) = optval;
        continue;
    end
    
    if ~isempty(optval_r) && optvals(1) <= optval_r && optval_r <= optvals(end-1)
        % replace worst with reflection
        indinsert = find(optval_r >= optvals(1:end-1),1,'last');
        optvals = [optvals(1:indinsert),optval_r,optvals(indinsert+1:end-1)];
        simplex = [simplex(:,1:indinsert),reflection,simplex(:,indinsert+1:end-1)];

        steps(ii) = 'r';
        if opt_spec.verbose, fprintf(1,'[%s] Iter %d: Reflect (optval = %g)\n',mfilename,ii,optval); drawnow; end
        optvalhist(ii+1) = optval;
        continue;
    end
    
    % contract instead
    if isempty(optval_r) || optval_r >= optvals(end) % contract inside
        contraction = (1-opt_spec.beta).*midpoint + opt_spec.beta.*simplex(:,end);
        contract_type = 'inside';
    else % contract outside
        contraction = (1-opt_spec.beta).*midpoint + opt_spec.beta.*reflection;
        contract_type = 'outside';
    end
    if isdegenerate([simplex(:,1:end-1),contraction])
        optval_c = [];
    else
        [outputs_test{:}] = reconfun(contraction); funcCount = funcCount + 1;
        optval_c = optfun(outputs_test{1});
    end

    if ~isempty(optval_c) && optval_c < optvals(end) % replace with contracted point
        indinsert = find(optval_c >= optvals(1:end-1),1,'last');
        if isempty(indinsert) % new minimum
            indinsert = 0;
            optval = optval_c;
            outputs = outputs_test;
            params = contraction;
            steps(ii) = upper(contract_type(1));
        else
            steps(ii) = lower(contract_type(1));
        end
        optvals = [optvals(1:indinsert),optval_c,optvals(indinsert+1:end-1)];
        simplex = [simplex(:,1:indinsert),contraction,simplex(:,indinsert+1:end-1)];

        if opt_spec.verbose, fprintf(1,'[%s] Iter %d: Contract %s (optval = %g)\n',mfilename,ii,contract_type,optval); drawnow; end
        optvalhist(ii+1) = optval;
        continue;
    end
    
    % shrink
    simplex(:,2:end) = bsxfun(@plus,(1-opt_spec.delta).*simplex(:,1),opt_spec.delta.*simplex(:,2:end));
    
    for n = 2:Nparams+1
        [outputs_test{:}] = reconfun(simplex(:,n)); funcCount = funcCount + 1;
        optvals(n) = optfun(outputs_test{1});
        if optvals(n) < optval
            % new minimum
            optval = optvals(n);
            outputs = outputs_test;
            params = simplex(:,n);
        end
    end
    
    % sort
    [optvals,indssort] = sort(optvals);
    simplex = simplex(:,indssort);

    if indssort(1) > 1
        steps(ii) = 'S';
    else
        steps(ii) = 's';
    end
    if opt_spec.verbose, fprintf(1,'[%s] Iter %d: Shrink (optval = %g)\n',mfilename,ii,optval); drawnow; end
    optvalhist(ii+1) = optval;
end

te = toc(t0);
info = struct('te',te,'flats',flats,'dims',dims,'optvals',optvalhist,'steps',steps,'funcCount',funcCount);

varargout = outputs;

end

function point = project_bounds(point,lbs,ubs)

below = isfinite(lbs) & point < lbs;
above = isfinite(ubs) & point > ubs;
point(below) = lbs(below);
point(above) = ubs(above);

end

function out = isdegenerate(simplex)

Nparams = size(simplex,1);

% construct edge matrix
edges = bsxfun(@minus,reshape(simplex.',Nparams+1,1,Nparams),reshape(simplex.',1,Nparams+1,Nparams));
edge_lengths = sqrt(sum(edges.^2,3));
max_edge_length = max(edge_lengths(:));
edge_lengths = edge_lengths + diag(Inf(Nparams+1,1));

if any(edge_lengths <= eps(max_edge_length))
    out = true; % overlapping points
    return;
end

% get determinant of coordinate matrix
coordmatrix = bsxfun(@minus,simplex(:,2:end),simplex(:,1));
volume = abs(det(coordmatrix))/factorial(Nparams);

if volume <= eps(prod(sqrt(sum(coordmatrix.^2,1))))
    out = true; % linear dependence among remaining vertices
else
    out = false;
end

end

function dim = getmaxdimension(simplex,lbs,ubs)

Nparams = size(simplex,1);

% construct edge matrix
edges = bsxfun(@minus,reshape(simplex.',Nparams+1,1,Nparams),reshape(simplex.',1,Nparams+1,Nparams));

% rescale bounded dimensions (others are assumed already normalized)
bounded = isfinite(lbs) & isfinite(ubs);
if any(bounded)
    edges(bounded) = bsxfun(@rdivide,edges(bounded),reshape(ubs(bounded)-lbs(bounded),1,1,sum(bounded)));
end

% return max of edge norms
dim = max(max(sqrt(sum(edges.^2,3))));

end
