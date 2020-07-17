# contentawarerecon
 Code and data for the paper "Content-aware compressive magnetic resonance image reconstruction"

Created by Daniel Weller, University of Virginia.

Please direct questions and comments to Daniel Weller at d.s.weller@ieee.org.

The contents of this program are subject to copyright as set forth in the file [LICENSE](LICENSE) accompanying this software. If you did not receive a copy of LICENSE, please email d.s.weller@ieee.org.

## I. Introduction

The enclosed program implements the content-aware reconstruction algorithm described in the manuscript 

D. S. Weller, M. Salerno, and C. H. Meyer, "Content-Aware Compressive Magnetic Resonance Image Reconstruction," *Magnetic Resonance Imaging*, vol. 52, pp. 118-130, October 2018. DOI: [10.1016/j.mri.2018.06.008](https://doi.org/10.1016/j.mri.2018.06.008)

We ask that any user please cite this work in any publications that use or modify this software in any way.

## II. Installation Instructions

This code has been tested with MATLAB (The Mathworks, Inc.) version R2015b. Please note that changes in more recent releases may affect the function of this software. This software is provided "AS IS", and the authors are not responsible for maintaining compatibility with other releases of MATLAB or other related software.

This code uses elements of version 0.3 of the ESPIRiT toolbox (see M. Uecker et al., MRM, 2014) provided from 

https://people.eecs.berkeley.edu/~mlustig/Software.html

as well as version 850 of Wavelab available from 

http://www-stat.stanford.edu/~wavelab

Installation instructions for these toolboxes can be found on their respective websites or in the readme files included with that software. Once these toolboxes are installed, either include the toolbox directories on your MATLAB path at startup, or be sure to call the setup_ESPIRiT.m and setup_Wavelab.m functions in this software before reconstructing. To make the best use of this software, some MEX files may need to be compiled for your software platform - see the toolbox documentation for more info.

Unzip the contents of this directory to your own computer and redirect the MATLAB working directory to this folder. You should now be able to run the reconstruction and data plotting scripts. 

## III. Running the Code

There are seven different reconstruction scripts, two for the Brainweb data (used for simulating convergence plots), two for the spiral brain data (with different parameter tuning approaches), and one for each of the other three datasets. The simulated acquisition parameters near the top of each script allow changing items like acceleration and SNR levels. These files are:

recon_brain_cartesian.m - reconstruct 13-channel 3D Cartesian brain volume  
recon_brain_spiral.m - reconstruct 8-channel T1-weighted brain data from retrospective undersampled spiral data with true error-based tuning  
recon_brain_spiral_WSURE.m - reconstruct 8-channel T1-weighted brain data from retrospective undersampled spiral data with WSURE-based tuning  
recon_brainweb_cartesian.m - reconstruct 8-channel Brainweb simulation  
recon_brainweb_cartesian_notune.m - reconstruct 8-channel Brainweb simulation using pretuned beta, lambda parameter values (such as those previously optimized by previously running recon_brainweb_cartesian.m)  
recon_cardiac_spiral.m - reconstruct seven frames of 3-channel spiral cardiac data from Xue Feng's 2016 MRM paper  
recon_cine_cartesian.m - reconstruct three-slice clinical cine video dataset

There are also two scripts for plotting results saved from running the above reconstruction scripts: plot_convergence.m and plot_performance.m. In the paper, Figs. 3 and 4 are produced from plot_convergence.m using the outputs of recon_brainweb_cartesian_notune.m with accel = 6 and 8, and selecting SNR = 20, 13 dB data. Figs. 5, 7, 9, and 11, and supplementary figure S5 are generated via plot_performance.m with the Cartesian brain, spiral brain, cardiac cine, and cardiac spiral reconstructions. The images in Figs. 2, 6, 8, 10, and 12, and supplementary figures S1-S4 are loaded from the saved reconstruction files. The Q-map in Fig. 1 is generated using the function compute_Q.m (with the ground truth image as input) and the default "window". The sparsity reweighting and content-aware reweighting use the appropriate functions for each transform. 

## IV. Adapting this Code

The ContentAwareRecon.m function contains the core algorithm of the reconstruction, and relies on a set of helper functions for setup and other purposes. Here is a list of helper functions from the provided code:

compute_Q.m - computes Q-map for an image (or image estimate)  
ContentAwareRecon.m - the core algorithm, called directly or via another helper function  
eig3x3symm.m - computes eigenvalues of 3x3 Hermitian matrix for fast Q computation  
entropythresh.m - one-dimensional entropy-based two-class thresholding used to compute Q-map  
ESPIRiT_kernels_calibrate.m - call appropriate ESPIRiT toolbox functions to calibrate SENSE maps  
estimate_WSURE_MC.m - compute WSURE estimate for a given reconstruction  
generate_WSURE_MC_Lambda.m - generate data-weighted optional Lambda matrix for WSURE  
GridSearchReconParameters.m - coarse-to-fine iterative grid-based optimization of regularization tuning parameters  
make_dcf.m - Voronoi-based density compensation factor computation (based on Fessler's Image Reconstruction Toolbox)  
make_DSFT_spec.m - construct operators for Cartesian and non-Cartesian Fourier transform, and circulant extension (if applicable)  
make_TTV_spec.m - construct temporal total variation operator for time series recons  
make_TV_spec.m - construct 2D spatial total variation operator  
make_TV3D_spec.m - construct 3D spatial or spatiotemporal total variation operator  
make_WAV_spec.m - construct orthogonal wavelet transform operator (requires Wavelab toolbox)  
mssim.m - compute mean structural similarity image index  
my_zpad_crop.m - crops or zero-pads an image, useful for dealing with Wavelab or selecting ESPIRiT calibration data  
OptimizeReconParameters.m - Nelder-Mead simplex-based optimization of regularization tuning parameters  
randpermw.m - weighted random permutation extension of MATLAB's randperm()  
setup_ESPIRiT.m - locate ESPIRiT toolbox if not on path  
setup_Wavelab.m - locate Wavelab toolbox if not on path  
spblkdiag.m - construct sparse block diagonal matrix from sparse submatrices  

Each file contains a brief description of its functionality.

Enjoy!

Daniel Weller  
Assistant Professor  
University of Virginia  
d.s.weller@ieee.org

June 7, 2018