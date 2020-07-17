% Test automatic parameter selection with content aware reconstruction on
% undersampled spiral T1-weighted brain data previously described by Weller
% et al., IEEE TMI 33(2), pp. 351-361, Feb. 2014. The image is shown in
% Figure 2(b), with Q and weight maps in Figure 1, PSER and MSSIM values in
% Figure 7, and reconstructions and difference images in Figure 8 (with
% accel = 3 and SNR = 10 dB).
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

beta_range = [1e-3,1e-1];
lambda_range = [1e-3,1e-1];
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

SURE_epsilon = 1e-3;
ngrid_params = 5;
niters_params = 3;

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

%% obtain calibration image data
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

errfun = @(x) norm(reshape(dispfun(abs(x)-abs(img)),[],1))/norm(reshape(dispfun(img),[],1));

%% set up DSFT (DFT operator for Cartesian data)
DSFTspec = make_DSFT_spec(Nx,Ny,omegas,w,Js,overgridding,NUFFTkernel);

%% set up wavelet, TV transforms (for individual slices)
WAVspec = make_WAV_spec([Nx,Ny],WAVtype,WAVpar,WAVscale,epsilon);
TVspec = make_TV_spec(epsilon);

%% output storage
errval_noQs = NaN(size(SNRs));
optval_noQs = NaN(size(SNRs));
beta_noQs = NaN(size(SNRs));
lambda_noQs = NaN(size(SNRs));
x_noQs = cell(size(SNRs));
info_sel_noQs = cell(size(SNRs));
info_noQs = cell(size(SNRs));

errval_RWs = NaN(size(SNRs));
optval_RWs = NaN(size(SNRs)); % unused
beta_RWs = NaN(size(SNRs));
lambda_RWs = NaN(size(SNRs));
x_RWs = cell(size(SNRs));
info_sel_RWs = cell(size(SNRs)); % unused
info_RWs = cell(size(SNRs));

errval_Qs = NaN(size(SNRs));
optval_Qs = NaN(size(SNRs)); % unused
beta_Qs = NaN(size(SNRs));
lambda_Qs = NaN(size(SNRs));
x_Qs = cell(size(SNRs));
info_sel_Qs = cell(size(SNRs)); % unused
info_Qs = cell(size(SNRs));

%% reconstruct for different noise levels
for ind = 1:numel(SNRs)
    %% prep
    Qs = {0,0};
    data_noisy = data+noisescales(ind).*noise;

    %% prepare WSURE estimation
    noisecov_scaled = noisescales(ind)^2.*noisecc;
    Lambda = generate_WSURE_MC_Lambda(data_noisy,DSFTspec);
    b = (1/sqrt(2)).*complex(2*(rand(size(data_noisy),class(data_noisy)) >= 0.5)-1, 2*(rand(size(data_noisy),class(data_noisy)) >= 0.5)-1); % binary rvs
    c = b./Lambda;
    c = reshape(reshape(c,[],NCha)*(noisecov_scaled.'),size(c));
    c = sum(bsxfun(@times,conj(Smat),DSFTspec.op_tr(bsxfun(@times,DSFTspec.weights,c))),4);
    
    %% tune reconstruction without Q
    wsurefun = @(x,params) estimate_WSURE_MC(x,[],data_noisy,DSFTspec,Smat,Lambda,noisecov_scaled,SURE_epsilon,b,c,Qs,exp(params),{WAVspec,TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
    [optparams,optval_noQs(ind),info_sel_noQs{ind},x_noQs{ind},info_noQs{ind}] = GridSearchReconParameters([beta_range;lambda_range],reconfun,wsurefun,struct('ngrid',ngrid_params,'niters',niters_params,'verbose',true));
    errval_noQs(ind) = errfun(x_noQs{ind});
    if ~save_complex, x_noQs{ind} = abs(x_noQs{ind}); end
    beta_noQs(ind) = exp(optparams(1)); 
    lambda_noQs(ind) = exp(optparams(2));
    fprintf(1,'[SNR=%g] Done with unweighted reconstruction (err = %g).\n',SNRs(ind),errval_noQs(ind)); drawnow;
    
    %% perform reweighted reconstruction
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useReweighted',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateReweighted',updateReweighted,'penalties',penalty,'Reweightfun',Reweightfun,'lastupdateReweighted',lastupdateReweighted));
    [x_RWs{ind},info_RWs{ind}] = reconfun(optparams);
    errval_RWs(ind) = errfun(x_RWs{ind});
    if ~save_complex, x_RWs{ind} = abs(x_RWs{ind}); end
    beta_RWs(ind) = exp(optparams(1)); % same as for unweighted
    lambda_RWs(ind) = exp(optparams(2)); % same as for unweighted
    fprintf(1,'[SNR=%g] Done with reweighted reconstruction (err = %g).\n',SNRs(ind),errval_RWs(ind)); drawnow;

    %% perform reconstruction with Q
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
    [x_Qs{ind},info_Qs{ind}] = reconfun(optparams);
    errval_Qs(ind) = errfun(x_Qs{ind});
    if ~save_complex, x_Qs{ind} = abs(x_Qs{ind}); end
    beta_Qs(ind) = exp(optparams(1)); % same as for unweighted
    lambda_Qs(ind) = exp(optparams(2)); % same as for unweighted
    fprintf(1,'[SNR=%g] Done with Q-weighted reconstruction with w(Q) = %s (err = %g).\n',SNRs(ind),func2str(Qweightfun),errval_Qs(ind)); drawnow;
end

%% cleanup
clear dcf ks Smat data_noisy noise b c Lambda wsurefun reconfun errfun DSFTspec WAVspec TVspec ans;
