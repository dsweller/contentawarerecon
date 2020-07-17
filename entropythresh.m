function t = entropythresh(x,nbins)
%
% "one-dimensional" entropy-based thresholding
% default nbins = 256
% See Abutaleb, Comp. Vis. Graph. Imag. Proc. 47, p. 22-32, 1989

if ~exist('nbins','var') || isempty(nbins), nbins = min(numel(x),256); end

x = x(:);
minval = min(x);
maxval = max(x);
delta = (maxval-minval)/nbins;

x = max(0,min(nbins-1,floor((x-minval)./delta)));
[x,~,inds] = unique(x,'sorted');
fs = accumarray(inds,1,size(x));
fs = fs./sum(fs);
fs = fs(1:end-1);
Fs = cumsum(fs);

hs = -fs.*log(fs);
Hs = cumsum(hs);

psis = log(Fs.*(1-Fs)) + Hs./Fs + (Hs(end)-Hs)./(1-Fs);

[~,indbest] = max(psis);
t = (x(indbest)+0.5).*delta+minval;

end
