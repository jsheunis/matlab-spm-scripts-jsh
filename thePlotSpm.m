% thePlotSpm
% 
% This is a Matlab script that uses custom code and spm12 routines to plot
% a version of THE PLOT from Jonathan Power. See:
% - https://www.jonathanpower.net/2017-ni-the-plot.html
% - https://www.ncbi.nlm.nih.gov/pubmed/27510328
% (code has been partly adopted from Power's script)
% 
% THE PLOT is a useful visual quality checking tool for fMRI timeseries
% data. It is a 2D plot of voxel intensities over time, with voxels ordered
% in bins (GM, WM, CSF) to indicate some directionality in the data. Motion
% and physiological noise are more easily visually inspected with this
% tool, especially if viewed next to other quality traces, like FD, DVARS
% or physiological recordings.

% My script does the following: 
%   - Run preprocessing on input data:
%       - Checks if T1w segments already exist.
%       - If not, realigns all images in time series to first functional
%       image, coregisters (estimate) T1w image to first functional image,
%       segments the coregistered T1w image into GM, WM and CSF
%       compartments, reslices all segments into functional image grid,
%       smooths realigned functional images, smooths unprocessed functional
%       images.
%   - Calculates framewise displacement
%   - Creates brain mask (GM+WM+CSF)
%   - Calculates BOLD percentage signal changes for masked brain voxels
%   - Creates figure with FD and THEPLOT (with blue and red lines
%     separating GM (top), WM (middle), and CSF(bottom)) 
% 
% INPUTS:
% structural_fn     - filename of T1-weighted structural scan, a 4D nifti
%                     file is expected (e.g. 'mprage.nii')
% funcional4D_fn    - filename of main functional scan, a 4D nifti file is
%                     expected (e.g. 'rest.nii')
% 
% DEPENDENCIES:
% thePlotSpmPreproc.m     - preprocesses data (functional realignment, coregistration of
%                         T1w to first functional image, segmentation into
%                         tissue types, reslicing to functional grid) 
% createBinarySegments.m  - constructs 3D binary images for GM, WM and CSF based
%                         on the relative value of GM/WM/CSF tissue probability
%                         maps per voxel.
% 
% NOTES/ASSUMPTIONS:
% - This script makes use of spm12 batch scripting
% - This script assumes that image files are in the Nifti file format
% - This script was developed in Matlab 2016b on a Mac OS Sierra, and not
% tested on other platforms (it wasn't even tested much on my Mac)
% - This is a simplified, Matlab-SPM-only version of Power's THEPLOT and
% excludes things like erosion of masks and segmentation into more bins
% (done in freesurfer). 
% - I might have made some mistakes in the code, or it could be improved in
% many ways. Please contribute!
%__________________________________________________________________________
% Copyright (C) 2018 - Stephan Heunis


% thePlotSpm

% User defined variables:
% -------------------------------------------------------------------------
data_dir = ''; % e.g. '/users/me/matlab/data/subj1'
functional4D_fn = [data_dir filesep '']; % e.g. [data_dir filesep 'rest.nii']
structural_fn = [data_dir filesep '']; % e.g. [data_dir filesep 'mprage.nii']
spm_dir = ''; % e.g. '/users/me/matlab/spm12'
fwhm = 6; % for preprocessing smoothing steps
use_processed = 0;  % use realigned data for the plot? yes = 1
intensity_scale = [-6 6]; % scaling for plot image intensity, see what works
% -------------------------------------------------------------------------

% Get image information
func_spm = spm_vol(functional4D_fn);
tsize = size(func_spm);
Nt = tsize(1);
Ni= func_spm(1).dim(1);
Nj= func_spm(1).dim(2);
Nk= func_spm(1).dim(3);

% Preprocess structural and functional data. First check if a resliced grey
% matter image is present in the specified data directory. If true, assume
% that all preprocessing has been completed by thePlotSpmPreproc, declare
% appropriate variables and 
[d, f, e] = fileparts(structural_fn);
if exist([d filesep 'rc1' f e], 'file')
    disp('Preproc done!')
    preproc_data = struct;
    preproc_data.forward_transformation = [d filesep 'y_' f e];
    preproc_data.inverse_transformation = [d filesep 'iy_' f e];
    preproc_data.gm_fn = [d filesep 'c1' f e];
    preproc_data.wm_fn = [d filesep 'c2' f e];
    preproc_data.csf_fn = [d filesep 'c3' f e];
    preproc_data.bone_fn = [d filesep 'c4' f e];
    preproc_data.soft_fn = [d filesep 'c5' f e];
    preproc_data.air_fn = [d filesep 'c6' f e];
    preproc_data.rstructural_fn = [d filesep 'r' f e];
    preproc_data.rgm_fn = [d filesep 'rc1' f e];
    preproc_data.rwm_fn = [d filesep 'rc2' f e];
    preproc_data.rcsf_fn = [d filesep 'rc3' f e];
    preproc_data.rbone_fn = [d filesep 'rc4' f e];
    preproc_data.rsoft_fn = [d filesep 'rc5' f e];
    preproc_data.rair_fn = [d filesep 'rc6' f e];
    
    [d, f, e] = fileparts(functional4D_fn);
    preproc_data.rfunctional_fn = [d filesep 'r' f e];
    preproc_data.srfunctional_fn = [d filesep 'sr' f e];
    preproc_data.sfunctional_fn = [d filesep 's' f e];
    preproc_data.mp_fn = [d filesep 'rp_' f '.txt'];
    preproc_data.MP = load(preproc_data.mp_fn);
else
    preproc_data = thePlotSpmPreproc(functional4D_fn, structural_fn, fwhm, spm_dir);
end


% Calculate Framewise displacement:
% - First demean and detrend the movement parameters and 
% - then calculate FD.
MP2 = preproc_data.MP - repmat(mean(preproc_data.MP, 1), Nt,1);
X = (1:Nt)';
X2 = X - mean(X);
X3 = [ones(Nt,1) X2];
b = X3\MP2;
MP_corrected = MP2 - X3(:, 2)*b(2,:);
MP_mm = MP_corrected; 
MP_mm(:,4:6) = MP_mm(:,4:6)*50; % 50mm = assumed brain radius; from article
MP_diff = [zeros(1, 6); diff(MP_mm)];
FD = sum(abs(MP_diff),2);

% Calculate BOLD percentage signal change:
% - First create binary masks for GM, WM and CSF,
% - Then combine them to get a brain mask for which calculations are run
% - Then compute PSC for unprocessed or realigned data (user selection)
% - Then plot everything
[GM_img_bin, WM_img_bin, CSF_img_bin] = createBinarySegments(preproc_data.rgm_fn, preproc_data.rwm_fn, preproc_data.rcsf_fn, 0.1);
I_GM = find(GM_img_bin);
I_WM = find(WM_img_bin);
I_CSF = find(CSF_img_bin);
mask_reshaped = GM_img_bin | WM_img_bin | CSF_img_bin;
I_mask = find(mask_reshaped);
Nmaskvox = numel(I_mask);

% Use realigned or unprocessed data for plot, user-selected
if use_processed == 1
    F_img = spm_read_vols(spm_vol(preproc_data.srfunctional_fn));
else
    F_img = spm_read_vols(spm_vol(preproc_data.sfunctional_fn));
end

% Remove linear and polynomial trends from data
F_2D = reshape(F_img, Ni*Nj*Nk, Nt);
X_design = [ (1:Nt)' ((1:Nt).^2/(Nt^2))' ((1:Nt).^3/(Nt^3))'];
X_design = X_design - mean(X_design);
X_design = [X_design ones(Nt,1)];
betas = X_design\F_2D';
F_detrended = F_2D' - X_design(:, 1:(end-1))*betas(1:(end-1), :);
F_detrended = F_detrended';

% Mask, mean and PSC
F_masked = F_detrended(I_mask, :);
F_mean = mean(F_masked, 2);
F_masked_psc = 100*(F_masked./repmat(F_mean, 1, Nt)) - 100;
F_masked_psc(isnan(F_masked_psc))=0;
F_psc_img = zeros(Ni, Nj, Nk, Nt);
F_2D_psc = reshape(F_psc_img, Ni*Nj*Nk, Nt);
F_2D_psc(I_mask, :) = F_masked_psc;
F_psc_img = reshape(F_2D_psc, Ni, Nj, Nk, Nt);

% Figure
GM_img = F_2D_psc(I_GM, :);
WM_img = F_2D_psc(I_WM, :);
CSF_img = F_2D_psc(I_CSF, :);
all_img = [GM_img; WM_img; CSF_img];
line1_pos = numel(I_GM);
line2_pos = numel(I_GM) + numel(I_WM);
figure; 
fontsizeL = 17;
fontsizeM = 15;
ax1 = subplot(5,1,[2:5]);
imagesc(ax1, all_img); colormap(gray); caxis(intensity_scale);
title(ax1, 'thePlotSpm','fontsize',fontsizeL)
ylabel(ax1, 'Voxels','fontsize',fontsizeM)
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
hold on; line([1 Nt],[line1_pos line1_pos],  'Color', 'b', 'LineWidth', 2 )
line([1 Nt],[line2_pos line2_pos],  'Color', 'r', 'LineWidth', 2 )
hold off;
ax2 = subplot(5,1,1);
plot(ax2, FD, 'LineWidth', 2); grid;
set(ax2,'Xticklabel',[]);
title(ax2, 'FD','fontsize',fontsizeL)
ylabel(ax2, 'mm','fontsize',fontsizeM)




