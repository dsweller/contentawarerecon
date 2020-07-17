% Test content aware reconstruction on undersampled spiral T1-weighted
% brain data previously described by Weller et al., IEEE TMI 33(2), pp.
% 351-361, Feb. 2014. The image is shown in Figure 2(b), with Q and weight
% maps in Figure 1, and PSER and MSSIM values in supplementary Figure S5.
% 
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

% path setup
setup_ESPIRiT;

setup_Wavelab;

rng('default');

filename = 'data/weller_spgr_spiral.mat';

%% simulated acquisition parameters
accel = 3; % range of 3 through 6
calibSize = [28,28]; % enough for 8 coil channels
SNRs = [3,6,10,13,16,20]; % amplify noise

%% reconstruction parameters
Js = [6,6];
overgridding = 2;
NUFFTkernel = 'kaiser';

kernelSize = [7,7];
eigThresh_k = 0.02;
Nmaps = 1;

beta_range = [1e-3,1e-1]; % to be consistent with range for WSURE-based tuning
lambda_range = [1e-3,1e-1]; % to be consistent with range for WSURE-based tuning
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

ks = double(ks);
dcf = double(dcf);
data = double(data);
noisecov = double(noisecov);

%% reconstruct ground truth image 
img = nufft_adj(reshape(bsxfun(@times,dcf,data),Nsamp*Nleaves,1,NCha),...
    nufft_init([real(ks(:)),imag(ks(:))].*(2*pi/Nx),[Nx,Ny],Js,overgridding*[Nx,Ny],[Nx,Ny]./2,NUFFTkernel))./sqrt(Nx*Ny);

%% obtain calibration image data (real scans might use actual calibration scan)
calibdata = reshape(my_zpad_crop(fftshift(fftshift(fft2(ifftshift(ifftshift(img,2),1)),2),1),calibSize),[calibSize,NCha])./sqrt(Nx*Ny);

%% undersample k-space
indleaves = sort(randperm(Nleaves,round(Nleaves/accel)));
omegas = ks(:,indleaves);
omegas = [real(omegas(:)),imag(omegas(:))];
M = size(omegas,1);

%% make DCF for each group
w = make_dcf(omegas,'prev');
omegas = omegas.*(2*pi/Nx);
fprintf(1,'Finished calculating density correction factors.\n'); drawnow;

%% subsample
data = reshape(data(:,indleaves,:),M,1,NCha);
fprintf(1,'Done generating data.\n'); drawnow;

%% create noise
noisecc = sqrt(diag(real(noisecov)));
noisecc = noisecov./(noisecc*noisecc');
noisescales = (norm(img(:))/sqrt(numel(img))).*10.^(-SNRs./20);

noise = reshape((sqrtm(noisecc/2)*complex(randn(NCha,M),randn(NCha,M))).',M,1,NCha);

%% do ESPIRiT calbiration
Smat = reshape(ESPIRiT_kernels_calibrate(calibdata,Nx,Ny,struct('kernelSize',kernelSize,'eigThresh_k',eigThresh_k,'Nmaps',Nmaps)),[Nx,Ny,1,NCha]);

%% combine image channels using SENSE maps and noise covariance
img = bsxfun(@rdivide,sum(bsxfun(@times,conj(Smat),reshape((noisecov\(reshape(img,Nx*Ny,NCha).')).',Nx,Ny,1,NCha)),4),real(sum(conj(Smat).*reshape((noisecov\(reshape(Smat,Nx*Ny,NCha).')).',Nx,Ny,1,NCha),4)));
if ~save_complex, img = abs(img); end

optfun = @(x) norm(reshape(dispfun(abs(x)-abs(img)),[],1))/norm(reshape(dispfun(img),[],1));

%% set up DSFT (DFT operator for Cartesian data)
DSFTspec = make_DSFT_spec(Nx,Ny,omegas,w,Js,overgridding,NUFFTkernel);

%% set up wavelet, TV transforms (for individual slices)
WAVspec = make_WAV_spec([Nx,Ny],WAVtype,WAVpar,WAVscale,epsilon);
TVspec = make_TV_spec(epsilon);

%% output storage
optval_noQs = NaN(size(SNRs));
beta_noQs = NaN(size(SNRs));
lambda_noQs = NaN(size(SNRs));
x_noQs = cell(size(SNRs));
info_sel_noQs = cell(size(SNRs));
info_noQs = cell(size(SNRs));

optval_RWs = NaN(size(SNRs));
beta_RWs = NaN(size(SNRs));
lambda_RWs = NaN(size(SNRs));
x_RWs = cell(size(SNRs));
info_sel_RWs = cell(size(SNRs));
info_RWs = cell(size(SNRs));

optval_Qs = NaN(size(SNRs));
beta_Qs = NaN(size(SNRs));
lambda_Qs = NaN(size(SNRs));
x_Qs = cell(size(SNRs));
info_sel_Qs = cell(size(SNRs));
info_Qs = cell(size(SNRs));

%% reconstruct for different noise levels
for ind = 1:numel(SNRs)
    %% prep
    Qs = {0,0};
    data_noisy = data+noisescales(ind).*noise;
    initsimplex = log([beta_range([1,1,2]);lambda_range([1,2,1])]);

    %% tune reconstruction without Q
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
    [optparams,optval_noQs(ind),info_sel_noQs{ind},x_noQs{ind},info_noQs{ind}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_noQs{ind} = abs(x_noQs{ind}); end
    beta_noQs(ind) = exp(optparams(1)); 
    lambda_noQs(ind) = exp(optparams(2));
    fprintf(1,'[SNR=%g] Done with unweighted reconstruction.\n',SNRs(ind)); drawnow;
    
    %% perform reweighted reconstruction
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useReweighted',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateReweighted',updateReweighted,'penalties',penalty,'Reweightfun',Reweightfun,'lastupdateReweighted',lastupdateReweighted));
    [optparams,optval_RWs(ind),info_sel_RWs{ind},x_RWs{ind},info_RWs{ind}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_RWs{ind} = abs(x_RWs{ind}); end
    beta_RWs(ind) = exp(optparams(1)); 
    lambda_RWs(ind) = exp(optparams(2));
    fprintf(1,'[SNR=%g] Done with reweighted reconstruction.\n',SNRs(ind)); drawnow;

    %% perform reconstruction with Q
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
    [optparams,optval_Qs(ind),info_sel_Qs{ind},x_Qs{ind},info_Qs{ind}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_Qs{ind} = abs(x_Qs{ind}); end
    beta_Qs(ind) = exp(optparams(1)); 
    lambda_Qs(ind) = exp(optparams(2));
    fprintf(1,'[SNR=%g] Done with Q-weighted reconstruction with w(Q) = %s.\n',SNRs(ind),func2str(Qweightfun)); drawnow;
end

%% cleanup
clear dcf ks Smat data_noisy noise reconfun optfun DSFTspec WAVspec TVspec ans;
