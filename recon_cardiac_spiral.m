% Test content aware reconstruction on undersampled continuous acquisition 
% described in Xue Feng et al, Mag Reson Med 75(4), pp. 1546-1555, 2016.
% (DOI: 10.1002/mrm.25738)
% Note this is a continuous acquisition across cardiac phases and is not a
% conventional cine sequence that shares data across heartbeats. This code
% uses retrospective spiral undersampling of reconstructed image, adding
% noise as desired. The image is shown in Figure 2(d), the PSER and MSSIM
% values are plotted in Figure 11, and an example reconstruction is shown
% in Figure 12 (accel = 3, SNR = 16 dB). 
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

% path setup
setup_ESPIRiT;

rng('default');

filename = 'data/xf4j_cardiac_spiral.mat';

%% simulated acquisition parameters
accel = 3; % reasonable range between 3 and 6
calibSize = [20,20]; % enough for 3 coil channels
SNRs = [20,16,13,10]; % in dB

%% reconstruction parameters
Js = [6,6];
overgridding = 2;
NUFFTkernel = 'kaiser';

kernelSize = [7,7];
eigThresh_k = 0.02;
Nmaps = 1;

lambdaT_range = [1e-3,1];
epsilon = 1e-6;

niters = 50;
alpha0 = 1;
c0 = 0.1;
alpha_tau = 0.5;
tolx = 0;
tolfun = 0;
updateQ = 25;
lastupdateQ = 50;
penalty = 'l1';
Qweightfun = @(Q) 1-Q;

updateReweighted = 25;
lastupdateReweighted = 50;
Reweightfun = @(w) 1./(sqrt(sum(abs(w).^2,4))+1);

niters_params = 50;
tolfun_params = 1e-4;
tolparams_params = 1e-4;

save_complex = false;

%% read data file
load(filename);

data = double(data); % square image so Nx = Ny
ks = double(ks);
dcf = double(dcf);
noisecov = double(noisecov);

%% reconstruct unaccelerated images (full set)
img = nufft_adj(bsxfun(@times,dcf(:),reshape(data,Nsamp*Nleaves,Nt,NCha)),...
    nufft_init([real(ks(:)),imag(ks(:))].*(2*pi/Nx),[Nx,Ny],Js,round(overgridding.*[Nx,Ny]),[Nx,Ny]./2,NUFFTkernel))./sqrt(Nx*Ny);

%% obtain calibration image data from first full group of interleaves
calibdata = reshape(my_zpad_crop(fftshift(fftshift(fft2(ifftshift(ifftshift(img(:,:,1,:),2),1)),2),1),calibSize),[calibSize,NCha])./sqrt(Nx*Ny);

%% rotate interleaves
rots = (0:round(Nleaves*Nt/accel)-1).*(2*pi/(1+sqrt(5))); % golden angle sampling across times
rotframes = floor((0:length(rots)-1).*(accel/Nleaves))+1;
indframes = arrayfun(@(ii) find(rotframes == ii),1:Nt,'UniformOutput',false);
if any(cellfun(@isempty,indframes)), error('Too undersampled: some frames are empty.'); end
omegas = [real(ks(:,1)),imag(ks(:,1))].*(2*pi/Nx);
omegas = cellfun(@(ii) reshape(cat(3,omegas*[cos(rots(ii));-sin(rots(ii))],omegas*[sin(rots(ii));cos(rots(ii))]),size(omegas,1)*length(ii),2),indframes,'UniformOutput',false);

%% make DCF for each group
w = cell(size(omegas));
for ind = 1:numel(omegas)
    fprintf(1,'Calculating density correction factor for frame %d.\n',ind); drawnow;
    w{ind} = make_dcf(omegas{ind}.*(Nx/(2*pi)),'prev');
end
fprintf(1,'Finished calculating density correction factors.\n'); drawnow;

%% sample using DSFT
[n1s,n2s] = ndgrid((0:Nx-1)-floor(Nx/2),(0:Ny-1)-floor(Ny/2));
n12s = [n1s(:),n2s(:)].';
data = cell(1,1,Nt);
for ic = 1:length(omegas)
    M = size(omegas{ic},1);
    Mblock = min(M,ceil(2^28/size(n12s,2)));
    data{ic} = zeros(M,NCha,'like',img);
    for m=1:Mblock:M
        fprintf(1,'Generating data for frame #%d/%d: %d/%d.\n',ic,length(omegas),ceil(m/Mblock),ceil(M/Mblock)); drawnow;
        mrange = m:min(M,m+Mblock-1);
        data{ic}(mrange,:) = exp(complex(0,-omegas{ic}(mrange,:)*n12s))*reshape(img(:,:,ic,:),Nx*Ny,NCha);
    end
end
data = cat(1,data{:}); % stack samples for different times
data = data./sqrt(Nx*Ny); % for unitary NUFFT scaling
fprintf(1,'Done generating data.\n'); drawnow;

%% set up DSFTs for preliminary period and all frames
DSFTspec = make_DSFT_spec(Nx,Ny,omegas,w,Js,overgridding,NUFFTkernel);
w = cat(1,w{:});

%% set up temporal TV transform
TTVspec = make_TTV_spec(epsilon);

%% generate noise
noisecc = sqrt(diag(real(noisecov)));
noisecc = noisecov./(noisecc*noisecc.');

noise = reshape((sqrtm(noisecc./2)*(reshape(complex(randn(size(data)),randn(size(data))),[],NCha).')).',size(data));
noisescales = (norm(img(:))/sqrt(numel(img))).*10.^(-SNRs./20);

%% do ESPIRiT calbiration
Smat = ESPIRiT_kernels_calibrate(calibdata,Nx,Ny,struct('kernelSize',kernelSize,'eigThresh_k',eigThresh_k,'Nmaps',Nmaps));
Smat = reshape(Smat,[Nx,Ny,1,NCha]);

%% combine image channels using SENSE maps and noise covariance
img = bsxfun(@rdivide,sum(bsxfun(@times,conj(Smat),reshape((noisecov\(reshape(img,Nx*Ny*Nt,NCha).')).',Nx,Ny,Nt,NCha)),4),real(sum(conj(Smat).*reshape((noisecov\(reshape(Smat,Nx*Ny,NCha).')).',Nx,Ny,1,NCha),4)));
if ~save_complex, img = abs(img); end

optfun = @(x) norm(reshape(dispfun(abs(x)-abs(img)),[],1))/norm(reshape(dispfun(img),[],1));

%% storage
lambdaT_noQs = NaN(size(SNRs));
lambdaT_RWs = NaN(size(SNRs));
lambdaT_Qs = NaN(size(SNRs));
optval_noQs = NaN(size(SNRs));
optval_RWs = NaN(size(SNRs));
optval_Qs = NaN(size(SNRs));
x_noQs = cell(size(SNRs));
x_RWs = cell(size(SNRs));
x_Qs = cell(size(SNRs));
info_noQs = cell(size(SNRs));
info_RWs = cell(size(SNRs));
info_Qs = cell(size(SNRs));
info_sel_noQs = cell(size(SNRs));
info_sel_RWs = cell(size(SNRs));
info_sel_Qs = cell(size(SNRs));

%% loop over noise amplifications
for ind = 1:length(SNRs)
    %% amplify noise
    data_noisy = data+noisescales(ind).*noise;
    
    %% tune reconstruction without Q
    Qs = {0};
    initsimplex = log(lambdaT_range);
    
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{TTVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
    [optparams,optval_noQs(ind),info_sel_noQs{ind},x_noQs{ind},info_noQs{ind}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_noQs{ind} = abs(x_noQs{ind}); end
    lambdaT_noQs(ind) = exp(optparams(1));
    fprintf(1,'[SNR=%g] Done tuning unweighted reconstruction.\n',SNRs(ind)); drawnow;
    
    %% tune reweighted reconstruction
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{TTVspec},struct('useReweighted',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateReweighted',updateReweighted,'penalties',penalty,'Reweightfun',Reweightfun,'lastupdateReweighted',lastupdateReweighted));
    [optparams,optval_RWs(ind),info_sel_RWs{ind},x_RWs{ind},info_RWs{ind}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_RWs{ind} = abs(x_RWs{ind}); end
    lambdaT_RWs(ind) = exp(optparams(1));
    fprintf(1,'[SNR=%g] Done with reweighted reconstruction.\n',SNRs(ind)); drawnow;

    %% tune reconstruction with Q
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,Smat,Qs,exp(params),{TTVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
    [optparams,optval_Qs(ind),info_sel_Qs{ind},x_Qs{ind},info_Qs{ind}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x_Qs{ind} = abs(x_Qs{ind}); end
    lambdaT_Qs(ind) = exp(optparams(1));
    fprintf(1,'[SNR=%g] Done tuning Q-weighted reconstruction with w(Q) = %s.\n',SNRs(ind),func2str(Qweightfun)); drawnow;
end

%% cleanup
clear n1s n2s dcf n12s ks Smat data_noisy noise reconfun optfun DSFTspec TTVspec ans;
