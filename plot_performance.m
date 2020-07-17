% Read in recons at different accelerations, SNR levels and plot
% info including PSERs and MSSIMs for results with and without content
% aware reconstruction. The trends versus acceleration are shown in Figures
% 5 and 7 (for Brainweb and MRXCAT), while the trends versus SNR are shown
% in Figures 8 and 10 (for real brain and cardiac recons).
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

[filenames,pathname] = uigetfile('*.mat','Select results files...','.','MultiSelect','on');
if isequal(filenames,0) || isequal(pathname,0), return; end
if ~iscell(filenames), filenames = {filenames}; end
Sread = cellfun(@(filename) load(fullfile(pathname,filename)),filenames,'UniformOutput',false);

%% get acceleration, SNR values
multiSNRs = cellfun(@(S) isfield(S,'SNRs'),Sread); % contains multiple SNRs
accels = cellfun(@(S) S.accel,Sread,'UniformOutput',false);
accels(multiSNRs) = cellfun(@(accel,S) accel.*ones(length(S.SNRs),1),accels(multiSNRs),Sread(multiSNRs),'UniformOutput',false);
SNRs = cell(size(Sread));
SNRs(~multiSNRs) = cellfun(@(S) S.SNR,Sread(~multiSNRs),'UniformOutput',false);
SNRs(multiSNRs) = cellfun(@(S) S.SNRs(:),Sread(multiSNRs),'UniformOutput',false);
[uniquevals,inds] = unique([cat(1,accels{:}),cat(1,SNRs{:})],'rows');
accels = uniquevals(:,1);
SNRs = uniquevals(:,2);

%% get unique accels and SNRs and arrange values into a grid
[accels,~,inds_accel] = unique(accels);
[SNRs,~,inds_SNR] = unique(SNRs);
accelstrs = arrayfun(@(accel) sprintf('%d',accel),accels,'UniformOutput',false);
SNRstrs = arrayfun(@(SNR) sprintf('%d',SNR),SNRs,'UniformOutput',false);
sz_grid = [length(accels),length(SNRs)];
inds_grid = sub2ind(sz_grid,inds_accel,inds_SNR);

%% recover images for unique tests
if isfield(Sread{1},'img')
    img = Sread{1}.img;
elseif isfield(Sread{1},'series')
    img = Sread{1}.series;
else
    error('Ground truth not found.');
end
dispfun = Sread{1}.dispfun;
clims = Sread{1}.clims;

if all(cellfun(@(S) any(isfield(S,{'x_Q','x_Qs'})),Sread))
    x_Qs = cell(size(Sread));
    x_Qs(~multiSNRs) = cellfun(@(S) {S.x_Q},Sread(~multiSNRs),'UniformOutput',false);
    x_Qs(multiSNRs) = cellfun(@(S) S.x_Qs(:),Sread(multiSNRs),'UniformOutput',false);
    x_Qs = cat(1,x_Qs{:});
    x_Qs = x_Qs(inds);
    x_Qs = subsasgn(cell(sz_grid),struct('type','()','subs',{{inds_grid}}),x_Qs);
elseif all(cellfun(@(S) any(isfield(S,{'series_Q','series_Qs'})),Sread))
    x_Qs = cell(size(Sread));
    x_Qs(~multiSNRs) = cellfun(@(S) {S.series_Q},Sread(~multiSNRs),'UniformOutput',false);
    x_Qs(multiSNRs) = cellfun(@(S) S.series_Qs(:),Sread(multiSNRs),'UniformOutput',false);
    x_Qs = cat(1,x_Qs{:});
    x_Qs = x_Qs(inds);
    x_Qs = subsasgn(cell(sz_grid),struct('type','()','subs',{{inds_grid}}),x_Qs);
else
    x_Qs = cell(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'x_RW','x_RWs'})),Sread))
    x_RWs = cell(size(Sread));
    x_RWs(~multiSNRs) = cellfun(@(S) {S.x_RW},Sread(~multiSNRs),'UniformOutput',false);
    x_RWs(multiSNRs) = cellfun(@(S) S.x_RWs(:),Sread(multiSNRs),'UniformOutput',false);
    x_RWs = cat(1,x_RWs{:});
    x_RWs = x_RWs(inds);
    x_RWs = subsasgn(cell(sz_grid),struct('type','()','subs',{{inds_grid}}),x_RWs);
elseif all(cellfun(@(S) any(isfield(S,{'series_RW','series_RWs'})),Sread))
    x_RWs = cell(size(Sread));
    x_RWs(~multiSNRs) = cellfun(@(S) {S.series_RW},Sread(~multiSNRs),'UniformOutput',false);
    x_RWs(multiSNRs) = cellfun(@(S) S.series_RWs(:),Sread(multiSNRs),'UniformOutput',false);
    x_RWs = cat(1,x_RWs{:});
    x_RWs = x_RWs(inds);
    x_RWs = subsasgn(cell(sz_grid),struct('type','()','subs',{{inds_grid}}),x_RWs);
else
    x_RWs = cell(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'x_noQ','x_noQs'})),Sread))
    x_noQs = cell(size(Sread));
    x_noQs(~multiSNRs) = cellfun(@(S) {S.x_noQ},Sread(~multiSNRs),'UniformOutput',false);
    x_noQs(multiSNRs) = cellfun(@(S) S.x_noQs(:),Sread(multiSNRs),'UniformOutput',false);
    x_noQs = cat(1,x_noQs{:});
    x_noQs = x_noQs(inds);
    x_noQs = subsasgn(cell(sz_grid),struct('type','()','subs',{{inds_grid}}),x_noQs);
elseif all(cellfun(@(S) any(isfield(S,{'series_noQ','series_noQs'})),Sread))
    x_noQs = cell(size(Sread));
    x_noQs(~multiSNRs) = cellfun(@(S) {S.series_noQ},Sread(~multiSNRs),'UniformOutput',false);
    x_noQs(multiSNRs) = cellfun(@(S) S.series_noQs(:),Sread(multiSNRs),'UniformOutput',false);
    x_noQs = cat(1,x_noQs{:});
    x_noQs = x_noQs(inds);
    x_noQs = subsasgn(cell(sz_grid),struct('type','()','subs',{{inds_grid}}),x_noQs);
else
    x_noQs = cell(sz_grid);
end

%% get parameters
if all(cellfun(@(S) any(isfield(S,{'beta_noQ','beta_noQs'})),Sread))
    beta_noQs = cell(size(Sread));
    beta_noQs(~multiSNRs) = cellfun(@(S) {S.beta_noQ},Sread(~multiSNRs),'UniformOutput',false);
    beta_noQs(multiSNRs) = cellfun(@(S) S.beta_noQs(:),Sread(multiSNRs),'UniformOutput',false);
    beta_noQs = cat(1,beta_noQs{:});
    beta_noQs = beta_noQs(inds);
    beta_noQs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),beta_noQs);
else
    beta_noQs = NaN(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'beta_RW','beta_RWs'})),Sread))
    beta_RWs = cell(size(Sread));
    beta_RWs(~multiSNRs) = cellfun(@(S) {S.beta_RW},Sread(~multiSNRs),'UniformOutput',false);
    beta_RWs(multiSNRs) = cellfun(@(S) S.beta_RWs(:),Sread(multiSNRs),'UniformOutput',false);
    beta_RWs = cat(1,beta_RWs{:});
    beta_RWs = beta_RWs(inds);
    beta_RWs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),beta_RWs);
else
    beta_RWs = NaN(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'beta_Q','beta_Qs'})),Sread))
    beta_Qs = cell(size(Sread));
    beta_Qs(~multiSNRs) = cellfun(@(S) {S.beta_Q},Sread(~multiSNRs),'UniformOutput',false);
    beta_Qs(multiSNRs) = cellfun(@(S) S.beta_Qs(:),Sread(multiSNRs),'UniformOutput',false);
    beta_Qs = cat(1,beta_Qs{:});
    beta_Qs = beta_Qs(inds);
    beta_Qs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),beta_Qs);
else
    beta_Qs = NaN(sz_grid);
end

%%
if all(cellfun(@(S) any(isfield(S,{'lambda_noQ','lambda_noQs'})),Sread))
    lambda_noQs = cell(size(Sread));
    lambda_noQs(~multiSNRs) = cellfun(@(S) {S.lambda_noQ},Sread(~multiSNRs),'UniformOutput',false);
    lambda_noQs(multiSNRs) = cellfun(@(S) S.lambda_noQs(:),Sread(multiSNRs),'UniformOutput',false);
    lambda_noQs = cat(1,lambda_noQs{:});
    lambda_noQs = lambda_noQs(inds);
    lambda_noQs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),lambda_noQs);
else
    lambda_noQs = NaN(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'lambda_RW','lambda_RWs'})),Sread))
    lambda_RWs = cell(size(Sread));
    lambda_RWs(~multiSNRs) = cellfun(@(S) {S.lambda_RW},Sread(~multiSNRs),'UniformOutput',false);
    lambda_RWs(multiSNRs) = cellfun(@(S) S.lambda_RWs(:),Sread(multiSNRs),'UniformOutput',false);
    lambda_RWs = cat(1,lambda_RWs{:});
    lambda_RWs = lambda_RWs(inds);
    lambda_RWs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),lambda_RWs);
else
    lambda_RWs = NaN(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'lambda_Q','lambda_Qs'})),Sread))
    lambda_Qs = cell(size(Sread));
    lambda_Qs(~multiSNRs) = cellfun(@(S) {S.lambda_Q},Sread(~multiSNRs),'UniformOutput',false);
    lambda_Qs(multiSNRs) = cellfun(@(S) S.lambda_Qs(:),Sread(multiSNRs),'UniformOutput',false);
    lambda_Qs = cat(1,lambda_Qs{:});
    lambda_Qs = lambda_Qs(inds);
    lambda_Qs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),lambda_Qs);
else
    lambda_Qs = NaN(sz_grid);
end

%%
if all(cellfun(@(S) any(isfield(S,{'lambdaT_noQ','lambdaT_noQs'})),Sread))
    lambdaT_noQs = cell(size(Sread));
    lambdaT_noQs(~multiSNRs) = cellfun(@(S) {S.lambdaT_noQ},Sread(~multiSNRs),'UniformOutput',false);
    lambdaT_noQs(multiSNRs) = cellfun(@(S) S.lambdaT_noQs(:),Sread(multiSNRs),'UniformOutput',false);
    lambdaT_noQs = cat(1,lambdaT_noQs{:});
    lambdaT_noQs = lambdaT_noQs(inds);
    lambdaT_noQs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),lambdaT_noQs);
else
    lambdaT_noQs = NaN(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'lambdaT_RW','lambdaT_RWs'})),Sread))
    lambdaT_RWs = cell(size(Sread));
    lambdaT_RWs(~multiSNRs) = cellfun(@(S) {S.lambdaT_RW},Sread(~multiSNRs),'UniformOutput',false);
    lambdaT_RWs(multiSNRs) = cellfun(@(S) S.lambdaT_RWs(:),Sread(multiSNRs),'UniformOutput',false);
    lambdaT_RWs = cat(1,lambdaT_RWs{:});
    lambdaT_RWs = lambdaT_RWs(inds);
    lambdaT_RWs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),lambdaT_RWs);
else
    lambdaT_RWs = NaN(sz_grid);
end

if all(cellfun(@(S) any(isfield(S,{'lambdaT_Q','lambdaT_Qs'})),Sread))
    lambdaT_Qs = cell(size(Sread));
    lambdaT_Qs(~multiSNRs) = cellfun(@(S) {S.lambdaT_Q},Sread(~multiSNRs),'UniformOutput',false);
    lambdaT_Qs(multiSNRs) = cellfun(@(S) S.lambdaT_Qs(:),Sread(multiSNRs),'UniformOutput',false);
    lambdaT_Qs = cat(1,lambdaT_Qs{:});
    lambdaT_Qs = lambdaT_Qs(inds);
    lambdaT_Qs = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),lambdaT_Qs);
else
    lambdaT_Qs = NaN(sz_grid);
end

%% other stuff
Ms = cell(size(Sread));
Ms(~multiSNRs) = cellfun(@(S) {S.M},Sread(~multiSNRs),'UniformOutput',false);
Ms(multiSNRs) = cellfun(@(S) S.M.*ones(length(S.SNRs),1),Sread(multiSNRs),'UniformOutput',false);
Ms = cat(1,Ms{:});
Ms = Ms(inds);
Ms = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),Ms);

sigma2s = cell(size(Sread));
sigma2s(~multiSNRs) = cellfun(@(S) {S.SNR},Sread(~multiSNRs),'UniformOutput',false);
sigma2s(multiSNRs) = cellfun(@(S) S.SNRs(:),Sread(multiSNRs),'UniformOutput',false);
sigma2s = cat(1,sigma2s{:});
sigma2s = sigma2s(inds);
sigma2s = subsasgn(NaN(sz_grid),struct('type','()','subs',{{inds_grid}}),sigma2s);
if isfield(Sread{1},'noisecc')
    sigma2s = trace(real(Sread{1}.noisecc)).*(norm(img(:))^2/numel(img)).*10.^(-sigma2s./10);
else
    sigma2s = (norm(img(:))^2/numel(img)).*10.^(-sigma2s./10);
end

%% compute PSNR, MSSIM values for magnitude images
PSER_Qs = NaN(sz_grid);
MSSIM_Qs = NaN(sz_grid);
if all(cellfun(@(S) any(isfield(S,{'x_Q','x_Qs','series_Q','series_Qs'})),Sread))
    PSER_Qs(inds_grid) = cellfun(@(x) 20*log10(clims(2))-20*log10(norm(reshape(dispfun(abs(abs(x)-abs(img))),[],1))/sqrt(numel(dispfun(img)))),x_Qs(inds_grid));
    MSSIM_Qs(inds_grid) = cellfun(@(x) mssim(dispfun(abs(x)).*255,dispfun(abs(img)).*255),x_Qs(inds_grid));
end
PSER_RWs = NaN(sz_grid);
MSSIM_RWs = NaN(sz_grid);
if all(cellfun(@(S) any(isfield(S,{'x_RW','x_RWs','series_RW','series_RWs'})),Sread))
    PSER_RWs(inds_grid) = cellfun(@(x) 20*log10(clims(2))-20*log10(norm(reshape(dispfun(abs(abs(x)-abs(img))),[],1))/sqrt(numel(dispfun(img)))),x_RWs(inds_grid));
    MSSIM_RWs(inds_grid) = cellfun(@(x) mssim(dispfun(abs(x)).*255,dispfun(abs(img)).*255),x_RWs(inds_grid));
end
PSER_noQs = NaN(sz_grid);
MSSIM_noQs = NaN(sz_grid);
if all(cellfun(@(S) any(isfield(S,{'x_noQ','x_noQs','series_noQ','series_noQs'})),Sread))
    PSER_noQs(inds_grid) = cellfun(@(x) 20*log10(clims(2))-20*log10(norm(reshape(dispfun(abs(abs(x)-abs(img))),[],1))/sqrt(numel(dispfun(img)))),x_noQs(inds_grid));
    MSSIM_noQs(inds_grid) = cellfun(@(x) mssim(dispfun(abs(x)).*255,dispfun(abs(img)).*255),x_noQs(inds_grid));
end

min_PSERs = min([PSER_Qs(:);PSER_RWs(:);PSER_noQs(:)]);
max_PSERs = max([PSER_Qs(:);PSER_RWs(:);PSER_noQs(:)]);
min_MSSIMs = min([MSSIM_Qs(:);MSSIM_RWs(:);MSSIM_noQs(:)]);
max_MSSIMs = max([MSSIM_Qs(:);MSSIM_RWs(:);MSSIM_noQs(:)]);

%% plot SNR trends
figure; subplot(1,2,1); 
hQs = plot(SNRs,PSER_Qs.','--o'); set(hQs,'LineWidth',2); hold on;
hRWs = plot(SNRs,PSER_RWs.','--o'); set(hRWs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
hnoQs = plot(SNRs,PSER_noQs.','-o'); set(hnoQs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
xlabel('SNR (dB)'); ylabel('PSER (dB)'); title('PSERs (dB) vs. SNR (dB)'); xlim([min(SNRs),max(SNRs)]);
set(gca,'XTick',SNRs);
[~,iorder] = sort(max(PSER_Qs,[],2),'descend');
legend(hQs(iorder),cellfun(@(str) ['R = ',str],accelstrs(iorder),'UniformOutput',false),'Location','best');
subplot(1,2,2); 
hQs = plot(SNRs,MSSIM_Qs.','--o'); set(hQs,'LineWidth',2); hold on;
hRWs = plot(SNRs,MSSIM_RWs.','--o'); set(hRWs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
hnoQs = plot(SNRs,MSSIM_noQs.','-o'); set(hnoQs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
xlabel('SNR (dB)'); ylabel('MSSIM'); title('MSSIMs vs. SNR (dB)'); xlim([min(SNRs),max(SNRs)]);
set(gca,'XTick',SNRs);
[~,iorder] = sort(max(MSSIM_Qs,[],2),'descend');
legend(hQs(iorder),cellfun(@(str) ['R = ',str],accelstrs(iorder),'UniformOutput',false),'Location','best');

%% plot acceleration trends
figure; subplot(1,2,1);
hQs = plot(accels,PSER_Qs,'--o'); set(hQs,'LineWidth',2); hold on;
hRWs = plot(accels,PSER_RWs,'--o'); set(hRWs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
hnoQs = plot(accels,PSER_noQs,'-o'); set(hnoQs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
xlabel('Acceleration'); ylabel('PSER (dB)'); title('PSERs (dB) vs. Acceleration'); xlim([min(accels),max(accels)]);
set(gca,'XTick',accels);
[~,iorder] = sort(max(PSER_Qs,[],1),'descend');
legend(hQs(iorder),cellfun(@(str) ['SNR = ',str,' dB'],SNRstrs(iorder),'UniformOutput',false),'Location','best');
subplot(1,2,2); 
hQs = plot(accels,MSSIM_Qs,'--o'); set(hQs,'LineWidth',2); hold on;
hRWs = plot(accels,MSSIM_RWs,'--o'); set(hRWs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
hnoQs = plot(accels,MSSIM_noQs,'-o'); set(hnoQs,'LineWidth',1,{'Color'},get(hQs,{'Color'}));
xlabel('Acceleration'); ylabel('MSSIM'); title('MSSIMs vs. Acceleration'); xlim([min(accels),max(accels)]);
set(gca,'XTick',accels);
[~,iorder] = sort(max(MSSIM_Qs,[],1),'descend');
legend(hQs(iorder),cellfun(@(str) ['SNR = ',str,' dB'],SNRstrs(iorder),'UniformOutput',false),'Location','best');

%% make strings for parameters
if all(cellfun(@(S) any(isfield(S,{'beta_noQ','beta_noQs'})),Sread))
    beta_noQs_str = sprintf('%0.3g\n',beta_noQs(:));
end
if all(cellfun(@(S) any(isfield(S,{'beta_RW','beta_RWs'})),Sread))
    beta_RWs_str = sprintf('%0.3g\n',beta_RWs(:));
end
if all(cellfun(@(S) any(isfield(S,{'beta_Q','beta_Qs'})),Sread))
    beta_Qs_str = sprintf('%0.3g\n',beta_Qs(:));
end
if all(cellfun(@(S) any(isfield(S,{'lambda_noQ','lambda_noQs'})),Sread))
    lambda_noQs_str = sprintf('%0.3g\n',lambda_noQs(:));
end
if all(cellfun(@(S) any(isfield(S,{'lambda_RW','lambda_RWs'})),Sread))
    lambda_RWs_str = sprintf('%0.3g\n',lambda_RWs(:));
end
if all(cellfun(@(S) any(isfield(S,{'lambda_Q','lambda_Qs'})),Sread))
    lambda_Qs_str = sprintf('%0.3g\n',lambda_Qs(:));
end
if all(cellfun(@(S) any(isfield(S,{'lambdaT_noQ','lambdaT_noQs'})),Sread))
    lambdaT_noQs_str = sprintf('%0.3g\n',lambdaT_noQs(:));
end
if all(cellfun(@(S) any(isfield(S,{'lambdaT_RW','lambdaT_RWs'})),Sread))
    lambdaT_RWs_str = sprintf('%0.3g\n',lambdaT_RWs(:));
end
if all(cellfun(@(S) any(isfield(S,{'lambdaT_Q','lambdaT_Qs'})),Sread))
    lambdaT_Qs_str = sprintf('%0.3g\n',lambdaT_Qs(:));
end