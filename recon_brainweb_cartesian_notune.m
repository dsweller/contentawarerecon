% Test content aware ESPIRiT reconstruction on Cartesian undersampled brain
% from Brainweb simulation with 8-channel simulated circular array coil
% from Fessler's image reconstruction toolbox. This code uses preselected
% lambda and beta regularization parameters (e.g., those tuned with the
% recon_brainweb_cartesian.m script) and does a full length reconstruction
% (no early termination) for convergence analysis. These are used to
% generate Figures 3 and 4 (with accel = 6 and accel = 8).
% 
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

% path setup
setup_ESPIRiT;

setup_Wavelab;

rng('default');

filename = 'data/brainweb_8ch.mat';

%% simulated acquisition parameters - these will be overriden by values loaded below
accel = 6; % try either 6 or 8
calibSize = [26,26];
vd_power = 1; % variable density
SNRs = [20,13]; % dB

%% reconstruction parameters
kernelSize = [7,7];
eigThresh_k = 0.02;
Nmaps = 1;

WAVtype = 'Daubechies'; % see MakeONFilter() for options
WAVpar = 4; % see MakeONFilter() for options
WAVscale = 4;
epsilon = 1e-6;

niters = 500;
alpha0 = 1;
c0 = 0.1;
alpha_tau = 0.5;
tolx = 0;
tolfun = 0;
updateQ = 25;
lastupdateQ = 250;
penalty = 'l1';
Qweightfun = @(Q) 1-Q;

save_complex = false;

%% read data file
load(filename);

img = double(img);
noisecov = double(noisecov);

%% get parameters
if ~exist('beta_Qs','var') || ~exist('beta_noQs','var') || ~exist('lambda_Qs','var') || ~exist('lambda_noQs','var') || any(cellfun(@isempty,{beta_Qs,beta_noQs,lambda_Qs,lambda_noQs}))
    [filename,pathname] = uigetfile('*.mat','Select results file to load parameters...','.');
    if isequal(filename,0) || isequal(pathname,0), return; end
    load(fullfile(pathname,filename),'beta_Qs','beta_noQs','lambda_Qs','lambda_noQs','accel','SNRs','calibSize','vd_power'); % match simulation parameters
end

%% subsample
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

%% set up wavelet, TV transforms
WAVspec = make_WAV_spec([Nx,Ny],WAVtype,WAVpar,WAVscale,epsilon);
TVspec = make_TV_spec(epsilon);

%% error metric to use
optfun = @(x) norm(reshape(dispfun(abs(x)-abs(img)),[],1))/norm(reshape(dispfun(img),[],1)); % NRMSE

%% storage
optval_Qs = NaN(size(SNRs));
x_Qs = cell(size(SNRs));
info_Qs = cell(size(SNRs));

optval_noQs = NaN(size(SNRs));
x_noQs = cell(size(SNRs));
info_noQs = cell(size(SNRs));

%% run reconstructions for each SNR level
for iSNR = 1:length(SNRs)
    Qs = {0,0};
    data_noisy = data+noisescales(iSNR).*noise;

    %% perform reconstruction with Q
    [x_Qs{iSNR},info_Qs{iSNR}] = ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,[beta_Qs(iSNR),lambda_Qs(iSNR)],{WAVspec,TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
    optval_Qs(iSNR) = optfun(x_Qs{iSNR});
    if ~save_complex, x_Qs{iSNR} = abs(x_Qs{iSNR}); end
    fprintf(1,'[SNR = %d] Done with Q-weighted reconstruction with w(Q) = %s (optval = %g).\n',SNRs(iSNR),func2str(Qweightfun),optval_Qs(iSNR)); drawnow;

    %% tune reconstruction without Q
    [x_noQs{iSNR},info_noQs{iSNR}] = ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,[beta_noQs(iSNR),lambda_noQs(iSNR)],{WAVspec,TVspec},struct('useQ',false,'niters',info_Qs{iSNR}.niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',0,'tolfun',0,'penalties',penalty));
    optval_noQs(iSNR) = optfun(x_noQs{iSNR});
    if ~save_complex, x_noQs{iSNR} = abs(x_noQs{iSNR}); end
    fprintf(1,'[SNR = %d] Done with unweighted reconstruction (optval = %g).\n',SNRs(iSNR),optval_noQs(iSNR)); drawnow;

end

%% cleanup
clear samp_outside subs psamp Smat data_noisy noise reconfun optfun DSFTspec WAVspec TVspec ans;
