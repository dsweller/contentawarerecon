% Test content aware ESPIRiT reconstruction on Cartesian undersampled brain
% from Weller et al., "Denoising Sparse Images from GRAPPA Using the
% Nullspace Method." Magn. Reson. Med., vol. 68, no. 4, pp. 1176-1189, Oct.
% 2012. This test automatically optimizes the regularization parameters for
% wavelet and total variation regularization. The tuned reconstruction
% PSERs and MSSIMs are shown in Figure 5, and reconstructed images are in
% Figure 6 and in supplementary Figures S1-S4. The ground truth is in
% Figure 2(a).
% 
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

% path setup
setup_ESPIRiT;

setup_Wavelab;

rng('default');

%% get file name for brain image (make sure calibration area appropriate for # of channels later)
filename = 'data/t1mprage_vol2_slc181_244_13ch.mat'; % contains 64/256 slices; download separately

%% simulated acquisition parameters
accel = 3; % reasonable range between 3 and 6
calibSize = [32,32]; % reasonable for 13 coil channels
vd_power = 1; % variable density
SNRs = [20,16,13,10,6,3]; % dB (reasonable range between 3 and 20 dB)

%% reconstruction parameters
kernelSize = [7,7];
eigThresh_k = 0.02;
Nmaps = 1;

beta_range = [1e-4,1e-1];
lambda_range = [1e-4,1e-1];
WAVtype = 'Daubechies'; % see MakeONFilter() for options
WAVpar = 4; % see MakeONFilter() for options
WAVscale = 4;
epsilon = 1e-6;

niters = 500;
alpha0 = 1;
c0 = 0.1;
alpha_tau = 0.5;
tolx = 1e-6;
tolfun = 1e-6;
updateQ = 25;
lastupdateQ = 250;
penalty = 'l1';
Qweightfun = @(Q) 1-Q;

updateReweighted = 25;
lastupdateReweighted = 250;
Reweightfun = @(w) 1./(sqrt(sum(abs(w).^2,4))+1);

niters_params = 50;
tolfun_params = 1e-4;
tolparams_params = 1e-4;

save_complex = false;

%% read data file
load(filename);

img = double(img);
noisecov = double(noisecov);

%% variable density subsampling (includes full center of k-space)
samp = my_zpad_crop(true(calibSize),[Nx,Ny]);
samp_outside = false(Nx*Ny-calibSize(1)*calibSize(2),1);
psamp = sqrt(bsxfun(@plus,(linspace(-1,1,Nx).').^2,linspace(-1,1,Ny).^2));
psamp = (1-psamp./max(psamp(:))).^vd_power;
psamp = psamp(~samp(:));
samp_outside(randpermw(psamp,round(Nx*Ny/accel-calibSize(1)*calibSize(2)))) = true;
samp(~samp) = samp_outside;

% readout direction (z) is fully sampled, so we don't bother FFT'ing in that direction
data = fftshift(fftshift(fft2(ifftshift(ifftshift(img,2),1)),2),1)./sqrt(Nx*Ny);
calibdata = my_zpad_crop(data,calibSize);
data = reshape(data,Nx*Ny,Nslices,NCha);
data = data(samp(:),:,:);
M = size(data,1);

%% Voronoi-based density correction factors
subs = zeros(M,2);
[subs(:,1),subs(:,2)] = ind2sub([Nx,Ny],find(samp(:)));
w = make_dcf(bsxfun(@minus,subs,1+floor([Nx,Ny]./2)),'max');

%% set up DSFT (DFT operator for Cartesian data)
DSFTspec = make_DSFT_spec(Nx,Ny,samp,w);

%% set up wavelet, total variation transforms
WAVspec = make_WAV_spec([Nx,Ny],WAVtype,WAVpar,WAVscale,epsilon);
TVspec = make_TV_spec(epsilon);

%% storage
optval_Qs = NaN(size(SNRs));
beta_Qs = NaN(size(SNRs));
lambda_Qs = NaN(size(SNRs));
x_Qs = arrayfun(@(SNR) zeros(Nx,Ny,Nslices),SNRs,'UniformOutput',false);
info_Qs = cell(size(SNRs));
info_sel_Qs = cell(size(SNRs));

optval_RWs = NaN(size(SNRs));
beta_RWs = NaN(size(SNRs));
lambda_RWs = NaN(size(SNRs));
x_RWs = arrayfun(@(SNR) zeros(Nx,Ny,Nslices),SNRs,'UniformOutput',false);
info_RWs = cell(size(SNRs));
info_sel_RWs = cell(size(SNRs));

optval_noQs = NaN(size(SNRs));
beta_noQs = NaN(size(SNRs));
lambda_noQs = NaN(size(SNRs));
x_noQs = arrayfun(@(SNR) zeros(Nx,Ny,Nslices),SNRs,'UniformOutput',false);
info_noQs = cell(size(SNRs));
info_sel_noQs = cell(size(SNRs));

%% single-slice parameter selection
iselect = floor(Nslices/2) + 1; % use middle by default

% do ESPIRiT calbiration for slice
Smat = reshape(ESPIRiT_kernels_calibrate(reshape(calibdata(:,:,iselect,:),[calibSize,NCha]),Nx,Ny,struct('kernelSize',kernelSize,'eigThresh_k',eigThresh_k,'Nmaps',Nmaps)),[Nx,Ny,1,NCha]);

% generate noise for slice
noisescales = (norm(img(:))/sqrt(Nx*Ny*Nslices*trace(noisecov))).*10.^(-SNRs./20);
noise = sqrtm(noisecov/2)*complex(randn([NCha,M],class(data)),randn([NCha,M],class(data)));
noise = reshape(noise.',M,1,NCha);

% ground truth
img_sos = zeros(Nx,Ny,Nslices);
if save_complex
    img_sos(:,:,iselect) = sum(conj(Smat).*reshape(img(:,:,iselect,:),Nx,Ny,1,NCha),4)./sum(abs(Smat).^2,4);
else
    img_sos(:,:,iselect) = abs(sum(conj(Smat).*reshape(img(:,:,iselect,:),Nx,Ny,1,NCha),4)./sum(abs(Smat).^2,4));
end

% error metric to use
optfun = @(x) norm(reshape(dispfun(abs(x)-abs(img_sos(:,:,iselect))),[],1))/norm(reshape(dispfun(img_sos(:,:,iselect)),[],1)); % NRMSE

% run tuning reconstructions for each SNR level
for iSNR = 1:length(SNRs)
    %% tune reconstruction without Q
    Qs = {0,0};
    data_noisy = data(:,iselect,:)+noisescales(iSNR).*noise;
    initsimplex = log([beta_range([1,1,2]);lambda_range([1,2,1])]);

    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
    [optparams,optval_noQs(iSNR),info_sel_noQs{iSNR},x,info_noQs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x = abs(x); end
    x_noQs{iSNR}(:,:,iselect) = x;
    beta_noQs(iSNR) = exp(optparams(1)); 
    lambda_noQs(iSNR) = exp(optparams(2));
    fprintf(1,'[SNR = %d] Done tuning unweighted reconstruction.\n',SNRs(iSNR)); drawnow;

    %% tune reweighted reconstruction
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useReweighted',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateReweighted',updateReweighted,'penalties',penalty,'Reweightfun',Reweightfun,'lastupdateReweighted',lastupdateReweighted));
    [optparams,optval_RWs(iSNR),info_sel_RWs{iSNR},x,info_RWs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x = abs(x); end
    x_RWs{iSNR}(:,:,iselect) = x;
    beta_RWs(iSNR) = exp(optparams(1)); 
    lambda_RWs(iSNR) = exp(optparams(2));
    fprintf(1,'[SNR = %d] Done tuning reweighted reconstruction.\n',SNRs(iSNR)); drawnow;

    %% tune reconstruction with Q
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
    [optparams,optval_Qs(iSNR),info_sel_Qs{iSNR},x,info_Qs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x = abs(x); end
    x_Qs{iSNR}(:,:,iselect) = x;
    beta_Qs(iSNR) = exp(optparams(1)); 
    lambda_Qs(iSNR) = exp(optparams(2));
    fprintf(1,'[SNR = %d] Done tuning Q-weighted reconstruction with w(Q) = %s.\n',SNRs(iSNR),func2str(Qweightfun)); drawnow;

end

%% recon other slices with tuned parameters
for islice = [1:iselect-1,iselect+1:Nslices]
    %% do ESPIRiT calbiration
    Smat = reshape(ESPIRiT_kernels_calibrate(reshape(calibdata(:,:,islice,:),[calibSize,NCha]),Nx,Ny,struct('kernelSize',kernelSize,'eigThresh_k',eigThresh_k,'Nmaps',Nmaps)),[Nx,Ny,1,NCha]);
    
    %% generate noise
    noise = sqrtm(noisecov/2)*complex(randn([NCha,M],class(data)),randn([NCha,M],class(data)));
    noise = reshape(noise.',M,1,NCha);
    
    %% ground truth
    if save_complex
        img_sos(:,:,islice) = sum(conj(Smat).*reshape(img(:,:,islice,:),Nx,Ny,1,NCha),4)./sum(abs(Smat).^2,4);
    else
        img_sos(:,:,islice) = abs(sum(conj(Smat).*reshape(img(:,:,islice,:),Nx,Ny,1,NCha),4)./sum(abs(Smat).^2,4));
    end
    
    %% error metric to use
    errfun = @(x) norm(reshape(dispfun(abs(x)-abs(img_sos(:,:,islice))),[],1))/norm(reshape(dispfun(img_sos(:,:,islice)),[],1)); % NRMSE
    
    %% run reconstructions for each SNR level
    for iSNR = 1:length(SNRs)
        %% prep
        Qs = {0,0};
        data_noisy = data(:,islice,:)+noisescales(iSNR).*noise;
        
        %% perform reconstruction without Q
        x = ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,[beta_noQs(iSNR),lambda_noQs(iSNR)],{WAVspec,TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
        if ~save_complex, x = abs(x); end
        x_noQs{iSNR}(:,:,islice) = x;
        fprintf(1,'[Slice %d, SNR = %d] Done with unweighted reconstruction (err = %g).\n',islice,SNRs(iSNR),errfun(x)); drawnow;
        
        %% perform reweighted reconstruction
        x = ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,[beta_RWs(iSNR),lambda_RWs(iSNR)],{WAVspec,TVspec},struct('useReweighted',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateReweighted',updateReweighted,'penalties',penalty,'Reweightfun',Reweightfun,'lastupdateReweighted',lastupdateReweighted));
        if ~save_complex, x = abs(x); end
        x_RWs{iSNR}(:,:,islice) = x;
        fprintf(1,'[Slice %d, SNR = %d] Done tuning reweighted reconstruction (err = %g).\n',islice,SNRs(iSNR),errfun(x)); drawnow;
        
        %% perform reconstruction with Q
        x = ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,[beta_Qs(iSNR),lambda_Qs(iSNR)],{WAVspec,TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
        if ~save_complex, x = abs(x); end
        x_Qs{iSNR}(:,:,islice) = x;
        fprintf(1,'[Slice %d, SNR = %d] Done with Q-weighted reconstruction with w(Q) = %s (err = %g).\n',islice,SNRs(iSNR),func2str(Qweightfun),errfun(x)); drawnow;
    end
end

%% cleanup
img = single(img_sos);
x_noQs = cellfun(@single,x_noQs,'UniformOutput',false);
x_RWs = cellfun(@single,x_RWs,'UniformOutput',false);
x_Qs = cellfun(@single,x_Qs,'UniformOutput',false);
clear samp_outside subs psamp Smat x data_noisy noise reconfun optfun DSFTspec WAVspec TVspec ans img_sos data calibdata;
