% Test content aware ESPIRiT reconstruction on Cartesian undersampled
% cardiac cine data obtained retrospectively from a clinical scan from
% anonymized DICOMs (now stored in .mat format).
% 
% This test automatically optimizes the regularization parameters for
% wavelet and total variation regularization for a single slice. The tuned
% parameters are used to reconstruct the other slices without further
% tuning. The complete reconstruction PSERs and MSSIMs are shown in Figure
% 9, and reconstructed images are in Figure 10 (for accel = 3, SNR = 16
% dB). The ground truth is in Figure 2(c).
% 
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

rng('default');

%% simulated acquisition parameters
accel = 3; % reasonable range between 3 and 6
vd_power = 1; % variable density
center_lines = 8; % dense center (note: no calibration here)
SNRs = [20,16,13,10]; % dB (reasonable range between 10 and 20 dB)

%% reconstruction parameters
lambda_range = [1e-4,1e-1];
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
filename = 'data/gj8an_cine_exc.mat'; % contains 64/256 slices

load(filename);

series = cellfun(@double,series,'UniformOutput',false);
series = cellfun(@(s) s./clims(2),series,'UniformOutput',false);
Nx = double(Nx); Ny = double(Ny);
clims = clims./clims(2);

%% variable density subsampling
samp = my_zpad_crop(true(center_lines,Nt),[Ny,Nt]);
samp_outside = false(Ny-center_lines,Nt);
psamp = abs(linspace(-1,1,Ny)).';
psamp = repmat((1-psamp./max(psamp(:))).^vd_power,1,Nt);
psamp = psamp(~samp(:));
samp_outside(randpermw(psamp,round(Ny*Nt/accel-center_lines*Nt))) = true;
samp(~samp) = samp_outside;
samp = repmat(shiftdim(samp,-1),[Nx,1,1]);

%% Voronoi-based density correction factors
try
    pool = gcp();
catch
    pool = [];
end
w = cell(1,Nt);
parfor t = 1:Nt
    samp1 = samp(:,:,t);
    subs = cell(1,2);
    [subs{:}] = ind2sub([Nx,Ny],find(samp1(:)));
    subs = [subs{:}];
    w{t} = make_dcf(bsxfun(@minus,subs,1+floor([Nx,Ny]./2)),'max');
end
delete(pool); clear pool;
w = cat(1,w{:});
fprintf(1,'Done computing density correction factors.\n'); drawnow;

%% set up DSFT (DFT operator for Cartesian data)
DSFTspec = make_DSFT_spec(Nx,Ny,samp,w);

%% set up temporal total variation transform
% TTVspec = make_TTV_spec(epsilon);
TVspec = make_TV3D_spec(epsilon);

%% get subsampled k-space data from time series for parameter selection
iselect = floor(Nslices/2) + 1; % use middle by default
data = fftshift(fftshift(fft2(ifftshift(ifftshift(series{iselect},2),1)),2),1)./sqrt(Nx*Ny);
data = data(samp(:));
M = size(data,1);

%% storage
optval_Qs = NaN(size(SNRs));
lambda_Qs = NaN(size(SNRs));
series_Qs = arrayfun(@(SNR) zeros(Nx,Ny,Nt,Nslices),SNRs,'UniformOutput',false);
info_Qs = cell(size(SNRs));
info_sel_Qs = cell(size(SNRs));

optval_RWs = NaN(size(SNRs));
lambda_RWs = NaN(size(SNRs));
series_RWs = arrayfun(@(SNR) zeros(Nx,Ny,Nt,Nslices),SNRs,'UniformOutput',false);
info_RWs = cell(size(SNRs));
info_sel_RWs = cell(size(SNRs));

optval_noQs = NaN(size(SNRs));
lambda_noQs = NaN(size(SNRs));
series_noQs = arrayfun(@(SNR) zeros(Nx,Ny,Nt,Nslices),SNRs,'UniformOutput',false);
info_noQs = cell(size(SNRs));
info_sel_noQs = cell(size(SNRs));

%% single-slice parameter selection

% generate noise for slice
noisescales = (norm(cellfun(@(img) norm(img(:)),series))/sqrt(Nx*Ny*Nt*Nslices)).*10.^(-SNRs./20);
noise = sqrt(1/2)*complex(randn([M,1],class(data)),randn([M,1],class(data)));

% ground truth
if ~save_complex
    series{iselect} = abs(series{iselect});
end

% error metric to use
optfun = @(x) norm(reshape(dispfun(abs(x)-abs(series{iselect})),[],1))/norm(reshape(dispfun(series{iselect}),[],1)); % NRMSE

% run tuning reconstructions for each SNR level
for iSNR = 1:length(SNRs)
    %% tune reconstruction without Q
    Qs = {0};
    data_noisy = data+noisescales(iSNR).*noise;
    initsimplex = log(lambda_range);

    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,1,Qs,exp(params),{TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
    [optparams,optval_noQs(iSNR),info_sel_noQs{iSNR},x,info_noQs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    if ~save_complex, x = abs(x); end
    series_noQs{iSNR}(:,:,:,iselect) = x;
    lambda_noQs(iSNR) = exp(optparams(1));
    fprintf(1,'[SNR = %d] Done tuning unweighted reconstruction.\n',SNRs(iSNR)); drawnow;

    %% tune reweighted reconstruction
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,1,Qs,exp(params),{TVspec},struct('useReweighted',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateReweighted',updateReweighted,'penalties',penalty,'Reweightfun',Reweightfun,'lastupdateReweighted',lastupdateReweighted));
%     [optparams,optval_RWs(iSNR),info_sel_RWs{iSNR},x,info_RWs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    [x,info_RWs{iSNR}] = reconfun(optparams);
    optval_RWs(iSNR) = optfun(x);
    if ~save_complex, x = abs(x); end
    series_RWs{iSNR}(:,:,:,iselect) = x;
    lambda_RWs(iSNR) = exp(optparams(1));
    fprintf(1,'[SNR = %d] Done tuning reweighted reconstruction.\n',SNRs(iSNR)); drawnow;

    %% tune reconstruction with Q
    reconfun = @(params) ContentAwareRecon([],data_noisy,DSFTspec,1,Qs,exp(params),{TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
%     [optparams,optval_Qs(iSNR),info_sel_Qs{iSNR},x,info_Qs{iSNR}] = OptimizeReconParameters(initsimplex,[],[],reconfun,optfun,struct('niters',niters_params,'verbose',true,'tolfun',tolfun_params,'tolparams',tolparams_params));
    [x,info_Qs{iSNR}] = reconfun(optparams);
    optval_Qs(iSNR) = optfun(x);
    if ~save_complex, x = abs(x); end
    series_Qs{iSNR}(:,:,:,iselect) = x;
    lambda_Qs(iSNR) = exp(optparams(1));
    fprintf(1,'[SNR = %d] Done tuning Q-weighted reconstruction with w(Q) = %s.\n',SNRs(iSNR),func2str(Qweightfun)); drawnow;

end

%% recon other slices with tuned parameters
for islice = [1:iselect-1,iselect+1:Nslices]
    %% get subsampled k-space data for this slice
    data = fftshift(fftshift(fft2(ifftshift(ifftshift(series{islice},2),1)),2),1)./sqrt(Nx*Ny);
    data = data(samp(:));
    
    %% generate noise
    noise = sqrt(1/2)*complex(randn([M,1],class(data)),randn([M,1],class(data)));
    
    %% ground truth
    if ~save_complex
        series{islice} = abs(series{islice});
    end
    
    %% error metric to use
    errfun = @(x) norm(reshape(dispfun(abs(x)-abs(series{islice})),[],1))/norm(reshape(dispfun(series{islice}),[],1)); % NRMSE
    
    %% run reconstructions for each SNR level
    for iSNR = 1:length(SNRs)
        %% prep
        Qs = {0};
        data_noisy = data+noisescales(iSNR).*noise;
        
        %% perform reconstruction without Q
        x = ContentAwareRecon([],data_noisy,DSFTspec,1,Qs,lambda_noQs(iSNR),{TVspec},struct('useQ',false,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'penalties',penalty));
        if ~save_complex, x = abs(x); end
        series_noQs{iSNR}(:,:,:,islice) = x;
        fprintf(1,'[Slice %d, SNR = %d] Done with unweighted reconstruction (err = %g).\n',islice,SNRs(iSNR),errfun(x)); drawnow;
        
        %% perform reweighted reconstruction
        x = ContentAwareRecon([],data_noisy,DSFTspec,1,Qs,lambda_RWs(iSNR),{TVspec},struct('useReweighted',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateReweighted',updateReweighted,'penalties',penalty,'Reweightfun',Reweightfun,'lastupdateReweighted',lastupdateReweighted));
        if ~save_complex, x = abs(x); end
        series_RWs{iSNR}(:,:,:,islice) = x;
        fprintf(1,'[Slice %d, SNR = %d] Done tuning reweighted reconstruction (err = %g).\n',islice,SNRs(iSNR),errfun(x)); drawnow;
        
        %% perform reconstruction with Q
        x = ContentAwareRecon([],data_noisy,DSFTspec,1,Qs,lambda_Qs(iSNR),{TVspec},struct('useQ',true,'niters',niters,'alpha0',alpha0,'c0',c0,'alpha_tau',alpha_tau,'tolx',tolx,'tolfun',tolfun,'updateQ',updateQ,'penalties',penalty,'Qweightfun',Qweightfun,'lastupdateQ',lastupdateQ));
        if ~save_complex, x = abs(x); end
        series_Qs{iSNR}(:,:,:,islice) = x;
        fprintf(1,'[Slice %d, SNR = %d] Done with Q-weighted reconstruction with w(Q) = %s (err = %g).\n',islice,SNRs(iSNR),func2str(Qweightfun),errfun(x)); drawnow;
    end
end

%% cleanup
series = single(cat(4,series{:}));
series_noQs = cellfun(@single,series_noQs,'UniformOutput',false);
series_RWs = cellfun(@single,series_RWs,'UniformOutput',false);
series_Qs = cellfun(@single,series_Qs,'UniformOutput',false);
clear samp_outside subs psamp x data_noisy noise reconfun optfun DSFTspec TVspec ans data calibdata;
