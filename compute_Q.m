function Q = compute_Q(x,cwindow)
% From input image series x, return the content measure corresponding to
% gradient structure Q for window around image voxels.
%
% Copyright (c) 2018, Daniel Weller, University of Virginia. All rights reserved.

if ~exist('cwindow','var') || isempty(cwindow), cwindow = 1; end

% central difference for interior voxels; forward/backward difference for
% boundary voxels
if size(x,3) > 1
    [gy,gx,gt] = gradient(x);
else
    [gy,gx] = gradient(x);
end

% compute elements to form Gram matrix of [gx(:),gy(:),gt(:)]
gxx = abs(gx).^2;
gyy = abs(gy).^2;
gxy = conj(gy).*gx;
if size(x,3) > 1
    gtt = abs(gt).^2;
    gxt = conj(gt).*gx;
    gyt = conj(gt).*gy;
end

h = ones(2*cwindow+1,1);
Gxx = gxx;
Gyy = gyy;
Gxy = gxy;
if size(x,3) > 1
    Gtt = gtt;
    Gxt = gxt;
    Gyt = gyt;
end
for idims = 1:ndims(x)
    Gxx = convn(Gxx,shiftdim(h,1-idims),'same');
    Gyy = convn(Gyy,shiftdim(h,1-idims),'same');
    Gxy = convn(Gxy,shiftdim(h,1-idims),'same');
    if size(x,3) > 1
        Gtt = convn(Gtt,shiftdim(h,1-idims),'same');
        Gxt = convn(Gxt,shiftdim(h,1-idims),'same');
        Gyt = convn(Gyt,shiftdim(h,1-idims),'same');
    end
end

% calculate eigenvalues
if size(x,3) == 1
    % closed form for eigenvalues for 2D matrix
    % l1 = [-(Gxx+Gyy)+sqrt((Gxx+Gyy)^2-4*(Gxx*Gyy-|Gxy|^2))]/2
    % l2 = [-(Gxx+Gyy)-sqrt((Gxx+Gyy)^2-4*(Gxx*Gyy-|Gxy|^2))]/2
    % l1 >= l2
    l2 = Gxx+Gyy;
    l1 = real(sqrt(l2.^2-4.*(Gxx.*Gyy-abs(Gxy).^2))); % should be nonnegative
    l2 = 0.5.*(l2-l1);
    l1 = l2 + l1;
else
    % closed form for eigenvalues for 3D matrix -- via trigonometric
    % identity (should already be sorted l1 >= l2 >= l3)
    [l1,l2,l3] = eig3x3symm(Gxx,Gxy,Gxt,Gyy,Gyt,Gtt);
end

% map eigenvalues to singular values and compute coherence
if size(x,3) == 1
    s1 = real(sqrt(l1)); s2 = real(sqrt(l2));
    sdiff = s1 - s2;
    ssum = s1 + s2;
else
    s1 = real(sqrt(l1)); s2 = real(sqrt(l2)); s3 = real(sqrt(l3));
    sdiff = s1+s2 - 2*s3;
    ssum = s1 + s2 + s3;
end

% threshold singular value difference to differentiate strong from weak
% structure - use entropy-based approach 
sdiff_threshold = entropythresh(sdiff);
Q = min(1,sdiff./sdiff_threshold).*(sdiff./ssum);

Q(~isfinite(Q)) = 0;

end
