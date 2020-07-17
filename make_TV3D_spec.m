function TV3Dspec = make_TV3D_spec(epsilon)
% Construct 3D total variation operator for use with Content Aware
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

TV3Dspec = struct('epsilon',epsilon);
TV3Dspec.op = @opfun;
TV3Dspec.op_tr = @optrfun;
TV3Dspec.op_Q = @Qfun;

end

function Q = Qfun(Q) % average Q-values for adjacent voxels in 3D (or space and time)

navg = ones(size(Q));
window = ones(2,2,2);
Q = Q(end:-1:1,end:-1:1,end:-1:1);
Q = convn(Q,window,'same');
navg = convn(navg,window,'same');
Q = Q./navg;
Q = Q(end:-1:1,end:-1:1,end:-1:1);

end

function out = opfun(in) % zero-boundary conditions

out1 = in - in([1,1:end-1],:,:);
out2 = in - in(:,[1,1:end-1],:);
out3 = in - in(:,:,[1,1:end-1]);

out = cat(4,out1,out2,out3);
    
end

function out = optrfun(in)

out = cat(1,-in(2,:,:,1),in(2:end-1,:,:,1)-in(3:end,:,:,1),in(end,:,:,1));
out = out + cat(2,-in(:,2,:,2),in(:,2:end-1,:,2)-in(:,3:end,:,2),in(:,end,:,2));
out = out + cat(3,-in(:,:,2,3),in(:,:,2:end-1,3)-in(:,:,3:end,3),in(:,:,end,3));
    
end
