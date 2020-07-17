function maps = ESPIRiT_kernels_calibrate(calibData,Nx,Ny,calib_specs)
% Learn coil sensitivity profiles from k-space calibration data via ESPIRiT
% method. This approach uses code from Michael Lustig's ESPIRiT toolbox
% (call setup_ESPIRiT() first).
%
% Inputs: calibData [Ncalx,Ncaly,Ncha]
%         Nx, Ny - image size
%         calib_specs structure
%
% Outputs: set of coil sensitivity map(s) for each coil channel

if ~isfield(calib_specs,'kernelSize') || isempty(calib_specs.kernelSize), calib_specs.kernelSize = [7,7]; end
if ~isfield(calib_specs,'eigThresh_k') || isempty(calib_specs.eigThresh_k), calib_specs.eigThresh_k = 0.02; end
if ~isfield(calib_specs,'Nmaps') || isempty(calib_specs.Nmaps), calib_specs.Nmaps = 1; end

% check for enough calibration data
Nfits = (size(calibData,1)-calib_specs.kernelSize(1)+1)*(size(calibData,2)-calib_specs.kernelSize(2)+1);
Ncoeffs = calib_specs.kernelSize(1)*calib_specs.kernelSize(2)*size(calibData,3)-1;
if Nfits < Ncoeffs
    warning('Not enough calibration data for ESPIRiT kernel calibration (%d fits < %d coefficients)',Nfits,Ncoeffs);
else
    fprintf(1,'[%s] Calibrating ESPIRiT kernels (%d coeffs) from %d fits.\n',mfilename,Ncoeffs,Nfits); drawnow;
end

[k,S] = dat2Kernel(calibData,calib_specs.kernelSize);
Nks = find(S >= S(1)*calib_specs.eigThresh_k,1,'last');
k = k(:,:,:,1:Nks);

[M,W] = kernelEig(k,[Nx,Ny]);

M = M(:,:,:,end-calib_specs.Nmaps+1:end); % maps arranged from smallest to largest eigenvalues
W = reshape(W(:,:,end-calib_specs.Nmaps+1:end),[Nx,Ny,1,calib_specs.Nmaps]);

maps = bsxfun(@times, M, sqrt(W));

end
