function [mssimval,ssimvals] = mssim(ref,tgt)
% be sure to check that ref, tgt are between 0..255!
% this code is based on Zhou Wang's ssim_index.m code
% (https://ece.uwaterloo.ca/~z70wang/research/ssim/ssim_index.m) 
% and is provided for convenience only. Can use MATLAB's ssim() function
% instead (for 2D) if using a new enough version of MATLAB.
%
% now handles 3D images (as sequence of 2D images)
%
% second output is SSIM map (like Figure 9(c,d))

K1 = 0.01; K2 = 0.03;

C1 = (K1*255).^2;
C2 = (K2*255).^2;

% construct Gaussian window (2D window is separable)
win1 = exp(-(-5:5).^2./(2*1.5^2));
win1 = win1./sum(win1);

% ensure images are floating point, real
if ~isfloat(ref), ref = single(ref); end
if ~isfloat(tgt), tgt = single(tgt); end
if ~isreal(ref), ref = abs(ref); end
if ~isreal(tgt), tgt = abs(tgt); end

% do convolution with window
mu1 = convn(ref,win1(:),'valid');
mu1 = convn(mu1,win1(:).','valid');
mu2 = convn(tgt,win1(:),'valid');
mu2 = convn(mu2,win1(:).','valid');

mu1_sq = mu1.^2;
mu2_sq = mu2.^2;
mu1_mu2 = mu1.*mu2;

sigma1_sq = convn(convn(ref.^2,win1(:),'valid'),win1(:).','valid') - mu1_sq;
sigma2_sq = convn(convn(tgt.^2,win1(:),'valid'),win1(:).','valid') - mu2_sq;
sigma12 = convn(convn(ref.*tgt,win1(:),'valid'),win1(:).','valid') - mu1_mu2;

ssimvals = ((2.*mu1_mu2 + C1).*(2.*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
mssimval = mean(ssimvals(:));

end
