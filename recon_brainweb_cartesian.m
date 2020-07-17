% Test content aware ESPIRiT reconstruction on Cartesian undersampled brain
% from Brainweb database with simulated 8-channel array coil created using
% Fessler's image reconstruction toolbox. This test automatically optimizes
% the regularization parameters for wavelet and total variation
% regularization. The other file (recon_brainweb_cartesian_notune.m) is
% used just to test convergence and uses the optimized parameters found
% from running this reconstruction. Real data are used instead of this data
% in the revised manuscript.
% 
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

% path setup
setup_ESPIRiT;

setup_Wavelab;

rng('default');

filename = 'data/brainweb_8ch.mat';

%% simulated acquisition parameters
accel = 6; % reasonable range between 6 and 10
calibSize = [26,26];
vd_power = 1; % variable density
SNRs = [20,16,13,10]; % dB (reasonable range between 10 and 20 dB)

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

data = fftshift(fftshift(fft2(ifftshift(ifftshift(img,2),1)),2),1)./sqrt(Nx*Ny);
calibdata = my_zpad_crop(data,calibSize);
data = reshape(data,Nx*Ny,1,NCha);
data = data(samp(:),:,:);
M = size(data,1);

%% Voronoi-based density correction factors
subs = zeros(M,2);
[subs(:,1),subs(:,2)] = ind2sub([Nx,Ny],find(samp(:)));
w = make_dcf(bsxfun(@minus,subs,1+floor([Nx,Ny]./2)),'max');

%% do ESPIRiT calbiration
Smat = ESPIRiT_kernels_calibrate(calibdata,Nx,Ny,struct('kernelSize',kernelSize,'eigThresh_k',eigThresh_k,'Nmaps',Nmaps));
Smat = reshape(Smat,[Nx,Ny,1,NCha]);

%% generate noise
noisescales = (norm(img(:))/sqrt(Nx*Ny*trace(noisecov))).*10.^(-SNRs./20);
noise = sqrtm(noisecov/2)*complex(randn([NCha,M],class(data)),randn([NCha,M],class(data)));
noise = reshape(noise.',M,1,NCha);

%% other stuff
img = sum(conj(Smat).*reshape(img,Nx,Ny,1,NCha),4)./sum(abs(Smat).^2,4);
if ~save_complex, img = abs(img); end

%% set up DSFT (DFT operator for Cartesian data)
DSFTspec = make_DSFT_spec(Nx,Ny,samp,w);

%% set up wavelet, total variation transforms
WAVspec = make_WAV_spec([Nx,Ny],WAVtype,WAVpar,WAVscale,epsilon);
TVspec = make_TV_spec(epsilon);

%% error metric to use
optfun = @(x) norm(reshape(dispfun(abs(x)-abs(img)),[],1))/norm(reshape(dispfun(img),[],1)); % NRMSE

%% storage
optval_Qs = NaN(size(SNRs));
beta_Qs = NaN(size(SNRs));
lambda_Qs = NaN(size(SNRs));
x_Qs = cell(size(SNRs));
info_Qs = cell(size(SNRs));
info_sel_Qs = cell(size(SNRs));

optval_noQs = NaN(size(SNRs));
beta_noQs = NaN(size(SNRs));
lambda_noQs = NaN(size(SNRs));
x_noQs = cell(size(SNRs));
info_noQs = cell(size(SNRs));
info_sel_noQs = cell(size(SNRs));

%% run reconstructions for each SNR level
for iSNR = 1:length(SNRs)
    %% tune reconstruction without Q
    Qs = {0,0};
    data_noisy = data+noisescales(iSNR).*noise;
    initsimplex = log([beta_range([1,1,2]);lambda_range([1,2,1])]);

    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
    [optparams,optval_noQs(iSNR),info_sel_noQs{iSNR},x_noQs{iSNR},info_noQs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_noQs{iSNR} = abs(x_noQs{iSNR}); end
    beta_noQs(iSNR) = exp(optparams(1)); 
    lambda_noQs(iSNR) = exp(optparams(2));
    fprintf(1,'[SNR = %d] Done with unweighted reconstruction.\n',SNRs(iSNR)); drawnow;

    %% perform reconstruction with Q
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{WAVspec,TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
    [optparams,optval_Qs(iSNR),info_sel_Qs{iSNR},x_Qs{iSNR},info_Qs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_Qs{iSNR} = abs(x_Qs{iSNR}); end
    beta_Qs(iSNR) = exp(optparams(1)); 
    lambda_Qs(iSNR) = exp(optparams(2));
    fprintf(1,'[SNR = %d] Done with Q-weighted reconstruction with w(Q) = %s.\n',SNRs(iSNR),func2str(Qweightfun)); drawnow;

end

%% cleanup
clear samp_outside subs psamp Smat data_noisy noise reconfun optfun DSFTspec WAVspec TVspec ans;
