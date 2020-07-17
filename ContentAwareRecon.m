function [x,info] = ContentAwareRecon(x0,data,DSFTspec,Smat,Qs,params,transform_specs,recon_spec)
% CS/SENSE reconstruction with (optional) content-aware weighted
% regularization (this is the core algorithm):
%
% x = argmin 0.5*||sqrt(data_weights).*(data-DSFTop*Smat*x)||^2
%             + beta*sum_p w(Q(p))*|[Psi*x]_p|
%             + lambda*sum_p w(Q(p))*||[TV(x)]_{p,:}||
%
% Approximate |[Psi*x]_p| with sqrt(|[Psi*x]_p|^2+epsilon)-sqrt(epsilon) 
% and ||[TV(x)]_{p,:}|| with sqrt(||[TV(x)]_{p,:}||^2+epsilon)-sqrt(epsilon)
%
% Then, the gradient is
% Smat'*DSFTop'*(data_weights.*(DSFTop*Smat*x-data))
% + beta*Psi'*(w(Q(p))./sqrt(|Psi*x|.^2+epsilon).*(Psi*x))
% + lambda*TV'(w(Q(p))./sqrt(|Dx*x|.^2+|Dy*x|.^2+|Dt*x|.^2+epsilon).*(TV(x)))
% where TV(x) = [Dx*x,Dy*x,Dt*x] are finite differences, and TV'(dx,dy,dt) = [Dx'*dx,Dy'*dy,Dt'*dt]
% and inner multiplication w(Q(p))./sqrt(|Dx*x|.^2+|Dy*x|.^2+|Dt*x|.^2+epsilon).*(TV(x))
% operates elementwise on all of dx = Dx*x, dy = Dy*x, and dt = Dt*x.
%
% We use nonlinear conjugate gradient with backtracking linesearch in this.
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

Reweight_epsilon = 1; % according to Candes' paper, select on the order of the nonzero coefficient values

if any(cellfun(@(spec) any(~isfield(spec,{'op','op_tr','op_Q','epsilon'})),transform_specs))
    error('Transform specifications not properly initialized.');
end

if ~isfield(recon_spec,'useReweighted') || isempty(recon_spec.useReweighted), recon_spec.useReweighted = false; end
if ~isfield(recon_spec,'useQ') || isempty(recon_spec.useQ), recon_spec.useQ = false; end
if recon_spec.useQ && recon_spec.useReweighted, error('Cannot use both Q and reweighted l1 simultaneously.'); end

if ~isfield(recon_spec,'niters') || isempty(recon_spec.niters), recon_spec.niters = 10; end
if ~isfield(recon_spec,'alpha0') || isempty(recon_spec.alpha0), recon_spec.alpha0 = 1; end
if ~isfield(recon_spec,'c0') || isempty(recon_spec.c0), recon_spec.c0 = 0.1; end
if ~isfield(recon_spec,'alpha_tau') || isempty(recon_spec.alpha_tau), recon_spec.alpha_tau = 0.5; end
if ~isfield(recon_spec,'tolx') || isempty(recon_spec.tolx), recon_spec.tolx = 0; end
if ~isfield(recon_spec,'tolfun') || isempty(recon_spec.tolfun), recon_spec.tolfun = 0; end
if ~isempty(params) && (~isfield(recon_spec,'penalties') || isempty(recon_spec.penalties)), recon_spec.penalties = 'l1'; end
if recon_spec.useQ && (~isfield(recon_spec,'updateQ') || isempty(recon_spec.updateQ)), recon_spec.updateQ = Inf; end
if recon_spec.useQ && (~isfield(recon_spec,'lastupdateQ') || isempty(recon_spec.lastupdateQ)), recon_spec.lastupdateQ = Inf; end
if recon_spec.useQ && (~isfield(recon_spec,'Qweightfun') || isempty(recon_spec.Qweightfun)), recon_spec.Qweightfun = @(Q) 1-Q; end
if recon_spec.useQ && (~isfield(recon_spec,'Qwindow') || isempty(recon_spec.Qwindow)), recon_spec.Qwindow = []; end
if recon_spec.useReweighted && (~isfield(recon_spec,'updateReweighted') || isempty(recon_spec.updateReweighted)), recon_spec.updateReweighted = Inf; end
if recon_spec.useReweighted && (~isfield(recon_spec,'lastupdateReweighted') || isempty(recon_spec.lastupdateReweighted)), recon_spec.lastupdateReweighted = Inf; end
if recon_spec.useReweighted && (~isfield(recon_spec,'Reweightfun') || isempty(recon_spec.Reweightfun)), recon_spec.Reweightfun = @(w) 1./(w+Reweight_epsilon); end

t0 = tic();

normsqdy = sum(reshape(bsxfun(@times,DSFTspec.weights,abs(data).^2),[],1));
AWy = sum(bsxfun(@times,conj(Smat),DSFTspec.op_tr(bsxfun(@times,DSFTspec.weights,data))),4);
if isempty(x0)
    Sgrammat = sum(abs(Smat).^2,4);
    x0 = bsxfun(@rdivide,AWy,Sgrammat); % normalize combined images for each time point
end

% initialize
x = x0;

% compute transforms
if isfield(DSFTspec,'op_gram')
    AWAx = sum(bsxfun(@times,conj(Smat),DSFTspec.op_gram(bsxfun(@times,Smat,x))),4);
else
    AWAx = sum(bsxfun(@times,conj(Smat),DSFTspec.op_tr(bsxfun(@times,DSFTspec.weights,DSFTspec.op(bsxfun(@times,Smat,x))))),4);
end

penalties = cell(size(params));
phis = cell(size(params));
if ~isempty(params)
    if ~iscell(recon_spec.penalties), recon_spec.penalties = repmat({recon_spec.penalties},1,length(params)); end
    for iparam = 1:length(params)
        switch lower(recon_spec.penalties{iparam})
            case 'l1'
                penalties{iparam} = @(w) sqrt(sum(abs(w).^2,4)+transform_specs{iparam}.epsilon)-sqrt(transform_specs{iparam}.epsilon);
                phis{iparam} = @(w) (sum(abs(w).^2,4)+transform_specs{iparam}.epsilon).^-0.5;
            otherwise
                error('Unsupported penalty function "%s".',recon_spec.penalties{iparam});
        end
    end
    Txs = cellfun(@(spec) spec.op(x),transform_specs,'UniformOutput',false);
else
    Txs = {};
end

% compute weights
nextupdate = [];
if recon_spec.useQ
    recalculate = cellfun(@isempty,Qs);
    if any(recalculate)
        Q = compute_Q(x,recon_spec.Qwindow);
        Qs(recalculate) = cellfun(@(spec) spec.op_Q(Q),transform_specs(recalculate),'UniformOutput',false);
    end
    weights = cellfun(@(Q) recon_spec.Qweightfun(Q),Qs,'UniformOutput',false);
    mean_weights = cellfun(@(w) mean(w(:)),weights);
    weights(mean_weights > 0) = arrayfun(@(ind) weights{ind}./mean_weights(ind),find(mean_weights > 0),'UniformOutput',false);
    
    if isfinite(recon_spec.updateQ) && recon_spec.updateQ >= 1
        nextupdate = recon_spec.updateQ;
        if nextupdate > recon_spec.lastupdateQ, nextupdate = []; end
    end
elseif recon_spec.useReweighted
    weights = cellfun(@(Tx) recon_spec.Reweightfun(sqrt(sum(abs(Tx).^2,4))),Txs,'UniformOutput',false);
    
    if isfinite(recon_spec.updateReweighted) && recon_spec.updateReweighted >= 1
        nextupdate = recon_spec.updateReweighted;
        if nextupdate > recon_spec.lastupdateReweighted, nextupdate = []; end
    end
else
    weights = num2cell(ones(size(params))); % no weighting
end

% get initial function value
fval = compute_objective(x,AWAx,Txs,weights,AWy,normsqdy,params,penalties);

% output storage
fvals = [fval;NaN(recon_spec.niters,1)];
alphas = NaN(recon_spec.niters,1);
gammas = NaN(recon_spec.niters,1);
xdiffs = NaN(recon_spec.niters,1);
fdiffs = NaN(recon_spec.niters,1);
tes = NaN(recon_spec.niters,1);
if isfield(recon_spec,'outfun') && isa(recon_spec.outfun,'function_handle')
    outvals = [recon_spec.outfun(x);NaN(recon_spec.niters,1)];
end

% get initial gradient
xgrad = compute_gradient(AWAx,Txs,weights,AWy,params,transform_specs,phis);
xgradprod = sum(abs(xgrad(:)).^2);
xdir = -xgrad; % initial descent direction
mval = -xgradprod;

gradreset = false;

% main loop
ii = 0;
while ii < recon_spec.niters
    ii = ii + 1;
    
    % gradient transforms
    if isfield(DSFTspec,'op_gram')
        AWAxdir = sum(bsxfun(@times,conj(Smat),DSFTspec.op_gram(bsxfun(@times,Smat,xdir))),4);
    else
        AWAxdir = sum(bsxfun(@times,conj(Smat),DSFTspec.op_tr(bsxfun(@times,DSFTspec.weights,DSFTspec.op(bsxfun(@times,Smat,xdir))))),4);
    end
    if ~isempty(params)
        Txdirs = cellfun(@(spec) spec.op(xdir),transform_specs,'UniformOutput',false);
    else
        Txdirs = {};
    end

    % line search
    alpha = recon_spec.alpha0;
    ftest = compute_objective(x+alpha.*xdir,AWAx+alpha.*AWAxdir,cellfun(@(Tx,Txdir) Tx+alpha.*Txdir,Txs,Txdirs,'UniformOutput',false),weights,AWy,normsqdy,params,penalties);
    normxdir = norm(xdir(:)); normx = norm(x(:));
    while ftest > fval + recon_spec.c0 * alpha * mval
        % shrink alpha and try again
        alpha = recon_spec.alpha_tau * alpha;
        ftest = compute_objective(x+alpha.*xdir,AWAx+alpha.*AWAxdir,cellfun(@(Tx,Txdir) Tx+alpha.*Txdir,Txs,Txdirs,'UniformOutput',false),weights,AWy,normsqdy,params,penalties);
        if alpha * normxdir <= eps(normx) % stagnation
            alpha = 0;
            ftest = fval;
            break; 
        end 
    end

    alphas(ii) = alpha;
    xdiffs(ii) = alpha*normxdir;
    fvals(ii+1) = ftest;

    % check for stagnation
    if alpha * normxdir <= recon_spec.tolx*normx && abs(fval - ftest) <= recon_spec.tolfun*abs(fval)
        if ~isempty(nextupdate)
            nextupdate = ii; % converged early, update Q again immediately
        else
            gammas(ii) = 0;
            if isfield(recon_spec,'outfun') && isa(recon_spec.outfun,'function_handle')
                outvals(ii+1) = outvals(ii);
            end
            tes(ii) = toc(t0);
            break;
        end
    end
    
    % update x, objective, and transforms (fast)
    x = x + alpha.*xdir;
    fval = ftest;
    AWAx = AWAx + alpha.*AWAxdir;
    Txs = cellfun(@(Tx,Txdir) Tx+alpha.*Txdir,Txs,Txdirs,'UniformOutput',false);
    if isfield(recon_spec,'outfun') && isa(recon_spec.outfun,'function_handle')
        outvals(ii+1) = recon_spec.outfun(x);
    end
    
    % update weights (optional)
    if ~isempty(nextupdate) && ii == nextupdate
        if recon_spec.useQ
            Q = compute_Q(x,recon_spec.Qwindow);
            Qs = cellfun(@(spec) spec.op_Q(Q),transform_specs,'UniformOutput',false);
            weights = cellfun(@(Q) recon_spec.Qweightfun(Q),Qs,'UniformOutput',false);
            mean_weights = cellfun(@(w) mean(w(:)),weights);
            weights(mean_weights > 0) = arrayfun(@(ind) weights{ind}./mean_weights(ind),find(mean_weights > 0),'UniformOutput',false);
        elseif recon_spec.useReweighted
            weights = cellfun(@(Tx) recon_spec.Reweightfun(sqrt(sum(abs(Tx).^2,4))),Txs,'UniformOutput',false);
        end
        
        % recompute function value and reset gradient
        fval = compute_objective(x,AWAx,Txs,weights,AWy,normsqdy,params,penalties);
        gradreset = true;
        
        fvals(ii+1) = fval;
        fdiffs(ii) = fval - ftest;
        if recon_spec.useQ
            nextupdate = nextupdate + recon_spec.updateQ;
            if nextupdate > recon_spec.lastupdateQ, nextupdate = []; end
        elseif recon_spec.useReweighted
            nextupdate = nextupdate + recon_spec.updateReweighted;
            if nextupdate > recon_spec.lastupdateReweighted, nextupdate = []; end
        end
    end
    
    % update gradient and descent direction
    xgradprev = xgrad;
    xgradprodprev = xgradprod;
    xgrad = compute_gradient(AWAx,Txs,weights,AWy,params,transform_specs,phis);
    xgradprod = sum(abs(xgrad(:)).^2);

    if gradreset
        xdir = -xgrad;
        gamma = 0;
        gradreset = false;
    else
        gamma = max(0,(xgradprod-sum(reshape(real(conj(xgrad).*xgradprev),[],1)))/xgradprodprev); % modified Polak-Ribiere
        xdir = -xgrad + gamma.*xdir;
    end

    mval = sum(reshape(real(conj(xdir).*xgrad),[],1));
    if mval > 0 % verify gradient direction
        xdir = -xgrad; % CG is in wrong direction
        gamma = 0;
        mval = -xgradprod;
    end
    
    gammas(ii) = gamma;
    tes(ii) = toc(t0);
    drawnow;
end

info = struct('niters',ii,'fvals',fvals,'tes',tes,'alphas',alphas,'gammas',gammas,'xdiffs',xdiffs,'fdiffs',fdiffs);
if isfield(recon_spec,'outfun') && isa(recon_spec.outfun,'function_handle')
    info.outvals = outvals;
end

end

function fval = compute_objective(x,AWAx,Txs,weights,AWy,normsqdy,params,penalties)

% data fit term = 0.5*(x'(AWAx) + normsqdy - 2Re{x'AWy})
fval = 0.5*(sum(reshape(real(conj(x).*(AWAx-2.*AWy)),[],1))+normsqdy);

% transforms
for iparam = 1:length(params)
    fval = fval + params(iparam)*sum(reshape(bsxfun(@times,weights{iparam},penalties{iparam}(Txs{iparam})),[],1));
end

end

function xgrad = compute_gradient(AWAx,Txs,weights,AWy,params,transform_specs,phis)

% data fit gradient A'W(Ax-y)
xgrad = AWAx-AWy;

% transforms
for iparam = 1:length(params)
    xgrad = xgrad + params(iparam).*transform_specs{iparam}.op_tr(bsxfun(@times,bsxfun(@times,weights{iparam},phis{iparam}(Txs{iparam})),Txs{iparam}));
end

end
