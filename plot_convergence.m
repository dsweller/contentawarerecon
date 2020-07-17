% Read in reconstructions and plot the convergence of NLCG algorithm. With
% the pretuned Brainweb reconstruction, this script will generate plots
% used in Figures 3 and 4 of the manuscript.
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

[filenames,pathname] = uigetfile('*.mat','Select results files...','.','MultiSelect','on');
if isequal(filenames,0) || isequal(pathname,0), return; end
if ~iscell(filenames), filenames = {filenames}; end
Sread = cellfun(@(filename) load(fullfile(pathname,filename)),filenames,'UniformOutput',false);

%% get info structures, other parameters
multiSNRs = cellfun(@(S) isfield(S,'SNRs'),Sread); % contains multiple SNRs
accels = cellfun(@(S) S.accel,Sread,'UniformOutput',false);
accels(multiSNRs) = cellfun(@(accel,S) accel.*ones(1,length(S.SNRs)),accels(multiSNRs),Sread(multiSNRs),'UniformOutput',false);
accels = [accels{:}];
SNRs = cell(size(Sread));
SNRs(~multiSNRs) = cellfun(@(S) S.SNR,Sread(~multiSNRs),'UniformOutput',false);
SNRs(multiSNRs) = cellfun(@(S) S.SNRs(:).',Sread(multiSNRs),'UniformOutput',false);
SNRs = [SNRs{:}];

info_Qs = cell(size(Sread));
info_Qs(~multiSNRs) = cellfun(@(S) S.info_Q,Sread(~multiSNRs),'UniformOutput',false);
info_Qs(multiSNRs) = cellfun(@(S) [S.info_Qs{:}],Sread(multiSNRs),'UniformOutput',false);
info_Qs = [info_Qs{:}];
info_noQs = cell(size(Sread));
info_noQs(~multiSNRs) = cellfun(@(S) S.info_noQ,Sread(~multiSNRs),'UniformOutput',false);
info_noQs(multiSNRs) = cellfun(@(S) [S.info_noQs{:}],Sread(multiSNRs),'UniformOutput',false);
info_noQs = [info_noQs{:}];

%% select subset
[sel,ok] = listdlg('Name','Test Selection Dialog...','SelectionMode','multiple','PromptString','Select tests to include:','ListString',arrayfun(@(accel,SNR) sprintf('R = %d, SNR = %g dB',accel,SNR),accels,SNRs,'UniformOutput',false),'InitialValue',1:length(SNRs));
if isempty(sel) || ~ok, return; end
accels = accels(sel);
SNRs = SNRs(sel);
info_Qs = info_Qs(sel);
info_noQs = info_noQs(sel);

%% plot convergence versus iterations, time
cols = sqrt(length(info_Qs));
rows = length(info_Qs)/cols;
if ceil(cols)*ceil(length(info_Qs)/ceil(cols)) > ceil(rows)*ceil(length(info_Qs)/ceil(rows))
    rows = ceil(rows);
    cols = ceil(length(info_Qs)/rows);
else
    cols = ceil(cols);
    rows = ceil(length(info_Qs)/cols);
end

figure;
for ii = 1:length(info_Qs)
    ax = subplot(rows,cols,ii);
    h1 = semilogy(1:info_Qs(ii).niters,info_Qs(ii).xdiffs(1:info_Qs(ii).niters),'r-','LineWidth',2);
    xlim([0,info_Qs(ii).niters]);
    hold on;
    h2 = semilogy(1:info_noQs(ii).niters,info_noQs(ii).xdiffs(1:info_noQs(ii).niters),'b-','LineWidth',1);
    if ii + cols > length(info_Qs), xlabel('Iterations'); end
    if mod(ii-1,cols) == 0, ylabel('|x^{(i)}-x^{(i-1)}|'); end
    title(sprintf('R = %d, SNR = %g dB',accels(ii),SNRs(ii)));
    if ii == cols, legend([h1,h2],{'Content-Aware','Conventional'},'Location','best'); end
    set(ax,'OuterPosition',[mod(ii-1,cols)/cols,1-floor((ii-1)/cols)/rows-1/rows,1/cols,1/rows]);
end

%%
figure;
for ii = 1:length(info_Qs)
    ax = subplot(rows,cols,ii);
    h1 = semilogy(info_Qs(ii).tes(1:info_Qs(ii).niters).',info_Qs(ii).xdiffs(1:info_Qs(ii).niters),'r-','LineWidth',2);
    xlim([0,info_Qs(ii).tes(info_Qs(ii).niters)]);
    hold on;
    h2 = semilogy(info_noQs(ii).tes(1:info_noQs(ii).niters).',info_noQs(ii).xdiffs(1:info_noQs(ii).niters),'b-','LineWidth',1);
    if ii + cols > length(info_Qs), xlabel('Time (s)'); end
    if mod(ii-1,cols) == 0, ylabel('|x^{(i)}-x^{(i-1)}|'); end
    title(sprintf('R = %d, SNR = %g dB',accels(ii),SNRs(ii)));
    if ii == cols, legend([h1,h2],{'Content-Aware','Conventional'},'Location','best'); end
end

%% 
tQs = arrayfun(@(info) diff([0;info.tes(1:info.niters)]),info_Qs,'UniformOutput',false);
tnoQs = arrayfun(@(info) diff([0;info.tes(1:info.niters)]),info_noQs,'UniformOutput',false);
tQs = cat(1,tQs{:});
tnoQs = cat(1,tnoQs{:});
fprintf('Computing time (s/iteration) for Content-Aware: %g +/- %g.\n',mean(tQs),std(tQs));
fprintf('Computing time (s/iteration) for Conventional: %g +/- %g.\n',mean(tnoQs),std(tnoQs));

%% plot function changes due to Q update
figure;
markers = {'o','*','x','s','d','^','v','>','<','p','h'};
markers = repmat(markers,1,ceil(length(sel)/length(markers)));
markers = markers(1:length(sel));
iters = find(isfinite(info_Qs(1).fdiffs));
plot(iters,abs(info_Qs(1).fdiffs(iters)),markers{1},'MarkerSize',8); hold on;
xlim([0,Sread{1}.lastupdateQ]);
xlabel('Iterations');
ylabel('|\Delta f(x^{(i)})|');
title('Change in objective function from updating \kappa''s');

for ii = 2:length(info_Qs)
    iters = find(isfinite(info_Qs(ii).fdiffs));
    plot(iters,abs(info_Qs(ii).fdiffs(iters)),markers{ii},'MarkerSize',8);
end
legend(arrayfun(@(accel,SNR) sprintf('R = %d, SNR = %g dB',accel,SNR),accels,SNRs,'UniformOutput',false),'Location','Best');
