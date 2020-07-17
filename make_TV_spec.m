function TVspec = make_TV_spec(epsilon)
% Construct 2D total variation operator for use with Content Aware
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

TVspec = struct('epsilon',epsilon);
TVspec.op = @opfun;
TVspec.op_tr = @optrfun;
TVspec.op_Q = @Qfun;

end

function Q = Qfun(Q) % average adjacent voxels in 2D

navg = ones(size(Q));
window = ones(2,2);
Q = Q(end:-1:1,end:-1:1,:);
Q = convn(Q,window,'same');
navg = convn(navg,window,'same');
Q = Q./navg;
Q = Q(end:-1:1,end:-1:1,:);

end

function out = opfun(in) % zero-boundary conditions

out1 = in - in([1,1:end-1],:,:);
out2 = in - in(:,[1,1:end-1],:);

out = cat(4,out1,out2);
    
end

function out = optrfun(in)

out = cat(1,-in(2,:,:,1),in(2:end-1,:,:,1)-in(3:end,:,:,1),in(end,:,:,1));
out = out + cat(2,-in(:,2,:,2),in(:,2:end-1,:,2)-in(:,3:end,:,2),in(:,end,:,2));
    
end
