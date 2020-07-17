function DSFTspec = make_DSFT_spec(Nx,Ny,samp_omegas,weights,varargin)
% Construct discrete space Fourier transform operator
% The input samp_omegas is either a logical matrix of 2D DFT samples or a
% Mx2 numeric list of NUFFT frequencies in radians.
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

if islogical(samp_omegas)
    % custom DFT operator (like Gdft without mask and with does_many == 1)
    DSFTspec = struct('type','DFT','M',sum(samp_omegas(:)),'samp',samp_omegas);
    if isempty(weights), weights = 1; end
    DSFTspec.weights = weights;
	DSFTspec.op = @(x) fft2_forward(x,DSFTspec);
    DSFTspec.op_tr = @(y) fft2_adjoint(y,DSFTspec);
	DSFTspec.op_gram = @(x) fft2_gram(x,DSFTspec);
else
    % DSFT gram matrix [AWA]_{no,ni} = sum_m w_m*exp{-j(om_m)(ni-no)}
    % is block Toeplitz, so double the size to make it block circulant:
    % t_{n1,n2} = sum_m w_m*exp{j(om_{m,1}n1+om_{m,2}n2)}, and
    %     [  t_{0,0}  ...   t_{0,N2-1}  0   t_{0,1-N2}  ...   t_{0,-1} ]
    %     [   ...              ...     ...     ...              ...    ]
    %     [t_{N1-1,0} ... t_{N1-1,N2-1} 0 t_{N1-1,1-N2} ... t_{N1-1,-1}]
    % C = [    0      ...       0       0       0       ...      0     ]
    %     [t_{1-N1,0} ... t_{1-N1,N2-1} 0 t_{1-N1,1-N2} ... t_{1-N1,-1}]
    %     [   ...              ...     ...     ...              ...    ]
    %     [ t_{-1,0}  ...   t{-1,N2-1}  0  t_{-1,1-N2}  ...  t_{-1,-1} ]
    if ~iscell(samp_omegas)
        samp_omegas = {samp_omegas};
        if ~isempty(weights)
            weights = {weights};
        end
    end
    Cs = complex(zeros((2*Nx-1)*Ny,length(samp_omegas)));
    [n1s,n2s] = ndgrid([0:Nx-1,1-Nx:-1],0:Ny-1);
    n12s = [n1s(:),n2s(:)];
    for ic = 1:length(samp_omegas)
        M = size(samp_omegas{ic},1);
        Nblock = min(size(n12s,1),ceil(2^28/M));
        for n=1:Nblock:size(n12s,1)
            fprintf(1,'Generating circulant extension for frame #%d/%d: %d/%d.\n',ic,length(samp_omegas),ceil(n/Nblock),ceil(size(n12s,1)/Nblock)); drawnow;
            nrange = n:min(size(n12s,1),n+Nblock-1);
            if isempty(weights)
                Cs(nrange,ic) = sum(exp(complex(0,n12s(nrange,:)*samp_omegas{ic}.')),2);
            else
                Cs(nrange,ic) = exp(complex(0,n12s(nrange,:)*samp_omegas{ic}.'))*weights{ic};
            end
        end
    end
    Cs = reshape(Cs,2*Nx-1,Ny,length(samp_omegas));
    Cs = cat(2,Cs(:,1:Ny,:),zeros(2*Nx-1,1,length(samp_omegas),'like',Cs),conj(Cs([1,end:-1:2],Ny:-1:2,:)));
    Cs = cat(1,Cs(1:Nx,:,:),zeros(1,2*Ny,length(samp_omegas),'like',Cs),Cs(Nx+1:end,:,:));
    Cs = real(fft2(Cs)); % enforce symmetry
    fprintf(1,'Done generating circulant extension(s).\n'); drawnow;
    
    % custom NUFFT operators (like Gnufft with fully parallel fft/ifft)
    if nargin >= 5
        Js = varargin{1};
    else
        Js = [6,6];
    end
    if nargin >= 6
        overgridding = varargin{2};
    else
        overgridding = 2;
    end
    if nargin >= 7
        NUFFTkernel = varargin{3};
    else
        NUFFTkernel = 'kaiser';
    end
    sts = cellfun(@(omegas) nufft_init(omegas,[Nx,Ny],Js,round(overgridding*[Nx,Ny]),[Nx,Ny]./2,NUFFTkernel),samp_omegas,'UniformOutput',false);
    DSFTspec = struct('type','NUFFT','N',sts{1}.Nd,'K',sts{1}.Kd);
    if length(samp_omegas) > 1
        DSFTspec.omegas = samp_omegas;
        DSFTspec.M = cellfun(@(st) st.M,sts);
        sns = cellfun(@(st) st.sn,sts,'UniformOutput',false);
        DSFTspec.sn = cat(3,sns{:});
        ps = cellfun(@(st) st.p,sts,'UniformOutput',false);
        DSFTspec.p = spblkdiag(ps{:});
        if ~isempty(weights)
            weights = cat(1,weights{:});
        else
            weights = 1;
        end
    else
        DSFTspec.omegas = samp_omegas{1};
        DSFTspec.M = sts{1}.M;
        DSFTspec.sn = sts{1}.sn;
        DSFTspec.p = sts{1}.p;
        if ~isempty(weights)
            weights = weights{1};
        else
            weights = 1;
        end
    end
    DSFTspec.weights = weights;
    DSFTspec.C = Cs;
    DSFTspec.op = @(x) nufft_forward(x,DSFTspec);
    DSFTspec.op_tr = @(y) nufft_adjoint(y,DSFTspec);
    DSFTspec.op_gram = @(x) dsft_gram(x,DSFTspec);
end

end

function y = fft2_forward(x,spec)

[Nx,Ny,Nt,NCha] = size(x);

y = fftshift(fftshift(fft2(ifftshift(ifftshift(x,2),1)),2),1)./sqrt(Nx*Ny);
if size(spec.samp,3) > 1
    y = reshape(y,Nx*Ny*Nt,NCha);
    y = y(spec.samp(:),:);
else
    y = reshape(y,Nx*Ny,Nt,NCha);
    y = y(spec.samp(:),:,:);
end
% y = bsxfun(@times,sqrt(spec.weights),y);

end

function x = fft2_adjoint(y,spec)

if size(spec.samp,3) > 1
    [Nx,Ny,Nt] = size(spec.samp);
    [~,NCha] = size(y);
else
    [Nx,Ny] = size(spec.samp);
    [~,Nt,NCha] = size(y);
end

% y = bsxfun(@times,sqrt(spec.weights),y);
if size(spec.samp,3) > 1
    x = zeros([Nx*Ny*Nt,NCha],'like',y);
    x(spec.samp(:),:) = y;
else
    x = zeros([Nx*Ny,Nt,NCha],'like',y);
    x(spec.samp(:),:,:) = y;
end

x = reshape(x,Nx,Ny,Nt,NCha);

x = sqrt(Nx*Ny).*fftshift(fftshift(ifft2(ifftshift(ifftshift(x,2),1)),2),1);

end

function x = fft2_gram(x,spec)

samp_weights = zeros(size(spec.samp),'like',spec.weights);
samp_weights(spec.samp) = spec.weights;

x = fft2(ifftshift(ifftshift(x,2),1));
x = bsxfun(@times,ifftshift(ifftshift(samp_weights,2),1),x);
x = fftshift(fftshift(ifft2(x),2),1);

end

function y = nufft_forward(x,spec)

[~,~,Nt,NCha] = size(x);

y = fft2(bsxfun(@times,spec.sn,x),spec.K(1),spec.K(2));
if iscell(spec.omegas)
    y = reshape(y,prod(spec.K)*Nt,NCha);
    y = cast(spec.p*double(y),'like',x)./sqrt(prod(spec.N));
else
    y = reshape(y,prod(spec.K),Nt*NCha);
    y = reshape(cast(spec.p*double(y),'like',x),[],Nt,NCha)./sqrt(prod(spec.N));
end

% if ~isempty(spec.weights)
%     y = bsxfun(@times,sqrt(spec.weights),y);
% end

end

function x = nufft_adjoint(y,spec)

if iscell(spec.omegas)
    Nt = length(spec.omegas);
    NCha = size(y,2);
else
    [~,Nt,NCha] = size(y);
end

% if ~isempty(spec.weights)
%     y = bsxfun(@times,sqrt(spec.weights),y);
% end

if iscell(spec.omegas)
    x = reshape(cast(spec.p'*double(y),'like',y),spec.K(1),spec.K(2),Nt,NCha);
else
    x = reshape(cast(spec.p'*double(reshape(y,[],Nt*NCha)),'like',y),spec.K(1),spec.K(2),Nt,NCha);
end
x = prod(spec.K).*ifft2(x);
x = x(1:spec.N(1),1:spec.N(2),:,:);
x = bsxfun(@times,conj(spec.sn),x)./sqrt(prod(spec.N));

end

function x = dsft_gram(x,spec)

[Nx,Ny,Nt,NCha] = size(x);

% recenter and zeropad
xpad = zeros(2*Nx,2*Ny,Nt,NCha,'like',x);
xpad(1:Nx,1:Ny,:,:) = x;
xpad = circshift(xpad,[-floor(Nx/2),-floor(Ny/2),0,0]);

% 2D fft-based convolution
xpad = bsxfun(@times,spec.C,fft2(xpad,2*Nx,2*Ny));
xpad = ifft2(xpad);

% recenter and keep just original [Nx,Ny] part of each image
xpad = circshift(xpad,[floor(Nx/2),floor(Ny/2),0,0]);
x = xpad(1:Nx,1:Ny,:,:)./(Nx*Ny);

end
