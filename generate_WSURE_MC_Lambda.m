function Lambda = generate_WSURE_MC_Lambda(data,DSFTspec)
% Constructs a suitable Lambda matrix diagonal for acquired k-space data.
% This process works by taking the data and finding a suitable "low-pass"
% approximation in 2D k-space. This works by gridding the data onto a
% low-resolution Cartesian grid via the NUFFT operator and resampling the
% low-resolution grid onto the acquired k-space locations. The Lambda
% matrix diagonal is normalized so the maximum value is 1, and the minimum
% value is 1e-4. This mitigates numerical instability.
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

Nx_new = 32; Ny_new = 32; % low FOV grid to smooth out k-space
data = bsxfun(@times,DSFTspec.weights,data);

switch upper(DSFTspec.type)
    case 'DFT'
        % get k-space coordinates for Cartesian subgridding
        if size(DSFTspec.samp,3) > 1
            [Nx,Ny,Nt] = size(DSFTspec.samp);
            [~,NCha] = size(data);
        else
            [Nx,Ny] = size(DSFTspec.samp);
            [~,Nt,NCha] = size(data);
        end
        [omegas_x,omegas_y] = ndgrid(((0:Nx-1) - floor(Nx/2)).*(2*pi/Nx),((0:Ny-1) - floor(Ny/2)).*(2*pi/Ny));
        
        % make smoothed image from k-space via NUFFT and then go back again
        if size(DSFTspec.samp,3) > 1
            % construct NUFFT operators
            samps = arrayfun(@(t) DSFTspec.samp(:,:,t),1:size(DSFTspec.samp,3),'UniformOutput',false);
            omegas = cellfun(@(samp) [omegas_x(samp(:)),omegas_y(samp(:))],samps,'UniformOutput',false);
            sts = cellfun(@(omega) nufft_init(omega,[Nx_new,Ny_new],[6,6],2.*[Nx_new,Ny_new],[Nx_new,Ny_new]./2,'kaiser'),omegas,'UniformOutput',false);
            sns = cellfun(@(st) st.sn,sts,'UniformOutput',false);
            sns = cat(3,sns{:});
            ps = cellfun(@(st) st.p,sts,'UniformOutput',false);
            ps = spblkdiag(ps{:});
            
            % recover subgrid images
            img = cast(reshape(ps'*double(data),[sts{1}.Kd(1),sts{1}.Kd(2),Nt,NCha]),'like',data);
            img = prod(sts{1}.Kd).*ifft2(img);
            img = img(1:sts{1}.Nd(1),1:sts{1}.Nd(2),:,:);
            img = bsxfun(@times,conj(sns),img)./sqrt(prod(sts{1}.Nd));

            % back to k-space
            img = fft2(bsxfun(@times,sns,img),sts{1}.Kd(1),sts{1}.Kd(2));
            img = reshape(img,prod(sts{1}.Kd)*Nt,NCha);
            Lambda = cast(ps*double(img),'like',data)./sqrt(prod(sts{1}.Nd));
        else
            omega = [omegas_x(DSFTspec.samp(:)),omegas_y(DSFTspec.samp(:))];
            st = nufft_init(omega,[Nx_new,Ny_new],[6,6],2.*[Nx_new,Ny_new],[Nx_new,Ny_new]./2,'kaiser');
            img = cast(reshape(st.p'*double(reshape(data,[],Nt*NCha)),[st.Kd(1),st.Kd(2),Nt,NCha]),'like',data);
            img = prod(st.Kd).*ifft2(img);
            img = img(1:st.Nd(1),1:st.Nd(2),:,:);
            img = bsxfun(@times,conj(st.sn),img)./sqrt(prod(st.Nd));
            
            % back to k-space
            img = fft2(bsxfun(@times,st.sn,img),st.Kd(1),st.Kd(2));
            img = reshape(img,prod(st.Kd),Nt*NCha);
            Lambda = reshape(cast(st.p*double(img),'like',data),[],Nt,NCha)./sqrt(prod(st.Nd));
        end
    case 'NUFFT'
        if iscell(DSFTspec.omegas)
            Nt = length(DSFTspec.omegas);
            NCha = size(data,2);
        else
            [~,Nt,NCha] = size(data);
        end

        % make smoothed image from k-space via NUFFT and then go back again
        if iscell(DSFTspec.omegas)
            % construct NUFFT operators
            sts = cellfun(@(omega) nufft_init(omega,[Nx_new,Ny_new],[6,6],2.*[Nx_new,Ny_new],[Nx_new,Ny_new]./2,'kaiser'),DSFTspec.omegas,'UniformOutput',false);
            sns = cellfun(@(st) st.sn,sts,'UniformOutput',false);
            sns = cat(3,sns{:});
            ps = cellfun(@(st) st.p,sts,'UniformOutput',false);
            ps = spblkdiag(ps{:});
            
            % recover subgrid images
            img = cast(reshape(ps'*double(data),[sts{1}.Kd(1),sts{1}.Kd(2),Nt,NCha]),'like',data);
            img = prod(sts{1}.Kd).*ifft2(img);
            img = img(1:sts{1}.Nd(1),1:sts{1}.Nd(2),:,:);
            img = bsxfun(@times,conj(sns),img)./sqrt(prod(sts{1}.Nd));

            % back to k-space
            img = fft2(bsxfun(@times,sns,img),sts{1}.Kd(1),sts{1}.Kd(2));
            img = reshape(img,prod(sts{1}.Kd)*Nt,NCha);
            Lambda = cast(ps*double(img),'like',data)./sqrt(prod(sts{1}.Nd));
        else
            st = nufft_init(DSFTspec.omegas,[Nx_new,Ny_new],[6,6],2.*[Nx_new,Ny_new],[Nx_new,Ny_new]./2,'kaiser');
            img = cast(reshape(st.p'*double(reshape(data,[],Nt*NCha)),[st.Kd(1),st.Kd(2),Nt,NCha]),'like',data);
            img = prod(st.Kd).*ifft2(img);
            img = img(1:st.Nd(1),1:st.Nd(2),:,:);
            img = bsxfun(@times,conj(st.sn),img)./sqrt(prod(st.Nd));
            
            % back to k-space
            img = fft2(bsxfun(@times,st.sn,img),st.Kd(1),st.Kd(2));
            img = reshape(img,prod(st.Kd),Nt*NCha);
            Lambda = reshape(cast(st.p*double(img),'like',data),[],Nt,NCha)./sqrt(prod(st.Nd));
        end
    otherwise
        Lambda = 1;
end

Lambda = abs(Lambda);
Lambda = Lambda./max(Lambda(:)); % max Lambda = 1
Lambda = max(Lambda,1e-4); % limit scaling

end
