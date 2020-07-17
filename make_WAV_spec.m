function WAVspec = make_WAV_spec(sz_x,FilterType,FilterPar,FilterScale,epsilon)
% Construct wavelet operator for use with Content Aware Reconstruction.
% This operator is returned as a struct with the following fields:
%
% .op() - forward wavelet transform
% .op_tr() - adjoint wavelet transform
% .epsilon - value added to make 1-norm differentiable (default=1e-6)
% .qmf - quadrature mirror filter used by orthogonal wavelet
% .scale - number of decompositions (scales) for wavelet transform
% .op_Q() - function for computing Q value
%
% Uses MakeONFilter, FWT2_PO, IWT2_PO from David Donoho's Wavelab
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

if ~exist('epsilon','var') || isempty(epsilon), epsilon = 1e-6; end

qmf = MakeONFilter(FilterType,FilterPar);

WAVspec = struct('epsilon',epsilon,'scale',FilterScale,'qmf',qmf);

% get appropriate size
sz_pad = 2^nextpow2(max(sz_x(1:2)));
sz_pad = [sz_pad,sz_pad];

WAVspec.op = @(x) opfun(x,sz_pad,FilterScale,qmf);
WAVspec.op_tr = @(w) optrfun(w,sz_x(1:2),FilterScale,qmf);
WAVspec.op_Q = @(Q) Qfun(Q,sz_pad,FilterScale);

end

function Q = Qfun(Q,sz_pad,scale)

% get appropriate size
Q = my_zpad_crop(Q,sz_pad);
Qpart = Q;

for j = 1:scale
    % average and downscale
    Qpart = 0.5.*(Qpart(1:2:end,:,:) + Qpart(2:2:end,:,:));
    Qpart = 0.5.*(Qpart(:,1:2:end,:) + Qpart(:,2:2:end,:));
    
    Q(1:2*size(Qpart,1),1:2*size(Qpart,2),:) = repmat(Qpart,[2,2,1]);
end

end

function w = opfun(x,sz_pad,scale,qmf)

[~,~,sz] = size(x);
x = my_zpad_crop(x,[sz_pad,sz]);
if isreal(x)
    w = cast(FWT2_PO(double(x(:,:,1)),scale,qmf),class(x)); %#ok<ZEROLIKE>
else
    w = complex(cast(FWT2_PO(double(real(x(:,:,1))),scale,qmf),class(x)),cast(FWT2_PO(double(imag(x(:,:,1))),scale,qmf),class(x))); %#ok<ZEROLIKE>
end
sz_w = [size(w),sz];
if prod(sz) > 1
    w(:,:,2:prod(sz)) = 0;
    for n = 2:prod(sz)
        if isreal(x)
            w(:,:,n) = cast(FWT2_PO(double(x(:,:,n)),scale,qmf),class(x)); %#ok<ZEROLIKE>
        else
            w(:,:,n) = complex(cast(FWT2_PO(double(real(x(:,:,n))),scale,qmf),class(x)),cast(FWT2_PO(double(imag(x(:,:,n))),scale,qmf),class(x))); %#ok<ZEROLIKE>
        end
    end
    w = reshape(w,sz_w);
end

end

function x = optrfun(w,sz_x,scale,qmf)

[~,~,sz] = size(w);
if isreal(w)
    x = cast(IWT2_PO(double(w(:,:,1)),scale,qmf),class(w)); %#ok<ZEROLIKE>
else
    x = complex(cast(IWT2_PO(double(real(w(:,:,1))),scale,qmf),class(w)),cast(IWT2_PO(double(imag(w(:,:,1))),scale,qmf),class(w))); %#ok<ZEROLIKE>
end
x = my_zpad_crop(x,sz_x);
sz_x = [sz_x,sz];
if prod(sz) > 1
    x(:,:,2:prod(sz)) = 0;
    for n = 2:prod(sz)
        if isreal(w)
            x(:,:,n) = my_zpad_crop(cast(IWT2_PO(double(w(:,:,n)),scale,qmf),class(w)),sz_x(1:2)); %#ok<ZEROLIKE>
        else
            x(:,:,n) = my_zpad_crop(complex(cast(IWT2_PO(double(real(w(:,:,n))),scale,qmf),class(w)),cast(IWT2_PO(double(imag(w(:,:,n))),scale,qmf),class(w))),sz_x(1:2)); %#ok<ZEROLIKE>
        end
    end
    x = reshape(x,sz_x);
end

end

