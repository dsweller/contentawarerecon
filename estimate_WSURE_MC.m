function WSURE_est = estimate_WSURE_MC(x_est,x0,data,DSFTspec,Smat,Lambda,noisecov,epsilon,b,c,Qs,params,transform_specs,recon_spec)
% Perform estimation of weighted predicted squared error via Monte Carlo
% WSURE (see Ramani et al., IEEE TMI, 2013).
%
% This code assumes you have already generated the random direction vector
% b as well as the vector c = S'*F'*W*Sigma*Lambda^{-1}*b. The Lambda
% matrix should be provided just by its diagonal. The epsilon value should
% be small to approximate the derivative, but not too small as to cause
% numerical issues.
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

% compute u(y), u(y+epsilon*Lambda*b)
difference = data + epsilon.*bsxfun(@times,Lambda,b);
difference = ContentAwareRecon(x0,difference,DSFTspec,Smat,Qs,params,transform_specs,recon_spec) - x_est;
difference = real(c(:)'*difference(:))./epsilon;

% get rest of WSURE estimate
WSURE_est = norm(reshape(bsxfun(@times,sqrt(DSFTspec.weights),DSFTspec.op(bsxfun(@times,Smat,x_est)) - data),[],1))^2 - trace(real(noisecov))*sum(DSFTspec.weights) + 2*difference;

end
