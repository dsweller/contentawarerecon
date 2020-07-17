function TTVspec = make_TTV_spec(epsilon)
% Construct temporal total variation operator for use with Content Aware
% Reconstruction. This operator is returned as a struct with the following
% fields:
%
% .op() - function that computes forward finite difference transform
% .op_tr() - function that computes adjoint finite difference transform
% .epsilon - value added to make TV differentiable at x = 0 (default=1e-6)
% .op_Q() - function for computing Q value
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

if ~exist('epsilon','var') || isempty(epsilon), epsilon = 1e-6; end

TTVspec = struct('epsilon',epsilon);
TTVspec.op = @opfun;
TTVspec.op_tr = @optrfun;
TTVspec.op_Q = @Qfun;

end

function Q = Qfun(Q) % average Q-values in adjacent voxels in time dimension

navg = ones(size(Q));
window = ones(1,1,2);
Q = Q(:,:,end:-1:1);
Q = convn(Q,window,'same');
navg = convn(navg,window,'same');
Q = Q./navg;
Q = Q(:,:,end:-1:1);

end

function out = opfun(in) % zero-boundary conditions

out = in - in(:,:,[1,1:end-1]);

end

function out = optrfun(in)

out = cat(3,-in(:,:,2),in(:,:,2:end-1)-in(:,:,3:end),in(:,:,end));
    
end
