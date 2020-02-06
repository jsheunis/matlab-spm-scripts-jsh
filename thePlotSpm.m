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
%
% 2018-08-17 Remi Gau 
% This seems to work on windows 10 with matlab 2017a.
% While looking around I realized that Cyril Pernet has a function 
% (spmup_timeseriesplot) to plot this (and do quite a few other things as 
% well) in his spmup toolbox: % https://github.com/CPernet/spmup
% But I like your script a lot as well as it is more stand-alone and
% therefore more pedagogical.
% 
%__________________________________________________________________________
% Copyright (C) 2018 - Stephan Heunis

clear
clc

%% User defined variables:
% -------------------------------------------------------------------------
data_dir = ''; % e.g. '/users/me/matlab/data/subj1'
functional4D_fn = [data_dir filesep '']; % e.g. [data_dir filesep 'rest.nii']
structural_fn = [data_dir filesep '']; % e.g. [data_dir filesep 'mprage.nii']
spm_dir = ''; % e.g. '/users/me/matlab/spm12'
fwhm = 6; % for preprocessing smoothing steps
use_processed = 0;  % use realigned data for the plot? yes = 1
intensity_scale = [-6 6]; % scaling for plot image intensity, see what works
davg_def = 50; %Average distance to the surface of the brain ; 50mm = assumed brain radius; from article
% -------------------------------------------------------------------------


% RG
% -------------------------------------------------------------------------
fwhm = 6;
[spm_dir, ~, ~] = fileparts(which('spm'));
data_dir = 'D:';
subj = 'sub-01';
% structural image
structural_fn = spm_select('FPList', fullfile(data_dir , subj , 'anat'), ['^' subj '_T1w.*nii$']);
% 4D BOLD time series
functional4D_fn = spm_select('FPList', fullfile(data_dir , subj , 'func'), ['^' subj '_task-.*_bold*.nii$']);
use_processed = 0;  % use realigned data for the plot? yes = 1
intensity_scale = [-6 6]; % scaling for plot image intensity, see what works
davg_def = 65;
% -------------------------------------------------------------------------


%% Get image information
func_spm = spm_vol(functional4D_fn);
tsize = size(func_spm);
Nt = tsize(1);
Ni= func_spm(1).dim(1);
Nj= func_spm(1).dim(2);
Nk= func_spm(1).dim(3);


%% Data check
% Preprocess structural and functional data. First check if a resliced grey
% matter image is present in the specified data directory. If true, assume
% that all preprocessing has been completed by thePlotSpmPreproc, declare
% appropriate variables and 
[dir, file, ext] = fileparts(structural_fn);
if exist([dir filesep 'rc1' file ext], 'file')
    disp('Preproc done!')
    preproc_data = struct;
    preproc_data.rstructural_fn = fullfile(dir , ['r' file ext]);
    preproc_data.rgm_fn = fullfile(dir , ['rc1' file ext]);
    preproc_data.rwm_fn = fullfile(dir , ['rc2' file ext]);
    preproc_data.rcsf_fn = fullfile(dir , ['rc3' file ext]);
    
    [dir, file, ext] = fileparts(functional4D_fn);
    preproc_data.rfunctional_fn = fullfile(dir , ['r' file ext]);
    preproc_data.srfunctional_fn = fullfile(dir , ['sr' file ext]);
    preproc_data.sfunctional_fn = fullfile(dir , ['s' file ext]);
    preproc_data.mp_fn = fullfile(dir , ['rp_' file '.txt']);
    preproc_data.MP = load(preproc_data.mp_fn);
else
    preproc_data = thePlotSpmPreproc(functional4D_fn, structural_fn, fwhm, spm_dir);
end

% Additional check to make sure that the functional images and the resliced
% TPMs have same dimensions
files = {};
files{1,1} = preproc_data.rgm_fn;
files{2,1} = [preproc_data.rfunctional_fn ',1'];   
files = char(files);
spm_check_orientations(spm_vol(files));


%% Try to compute the average surface to the brain
% from motion finger print functions and script (mw_mfp.m) from Marko Wilke
% https://www.medizin.uni-tuebingen.de/kinder/en/research/neuroimaging/software/
% see http://www.dx.doi.org/10.1016/j.neuroimage.2011.10.043
% and http://www.dx.doi.org/10.1371/journal.pone.0106498
[dir, file, ext] = fileparts(structural_fn);
try
if isempty(spm_select('FPList',dir,'^rc1.*\.surf\.gii$'))
    % create a surface of the brain using the TPM  if we don't have one
    files = spm_select('FPList',dir,'^rc[12].*\.nii$');
    spm_surf(files,2);
    clear files
end
% compute the average surface to the brain
FV = gifti(spm_select('FPList',dir,'^rc1.*\.surf\.gii$'));
center = FV.vertices(FV.faces(:,:),:);
center = reshape(center, [size(FV.faces,1) 3 3]);
center = squeeze(mean(center,2));
ori_dist = sqrt(sum((center.*-1).^2,2))';
davg = mean(ori_dist);
clear ori_dist center FV

catch
    davg = davg_def;
end
fprintf(' Average distance to the cortex surface: %0.2f ', davg)


%% Calculate Framewise displacement:
% - First demean and detrend the movement parameters and 
% - then calculate FD.
% See also spmup_FD from https://github.com/CPernet/spmup
MP2 = preproc_data.MP - repmat(mean(preproc_data.MP, 1), Nt,1);
X = (1:Nt)';
X2 = X - mean(X);
X3 = [ones(Nt,1) X2];
b = X3\MP2;
MP_corrected = MP2 - X3(:, 2)*b(2,:);
MP_mm = MP_corrected; 
MP_mm(:,4:6) = MP_mm(:,4:6)*davg; 
MP_diff = [zeros(1, 6); diff(MP_mm)];
FD = sum(abs(MP_diff),2);


%% Calculate BOLD percentage signal change (PSC)
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


%% Load and detrend data
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


%% Figure
% Create image to plot by concatenation
GM_img = F_2D_psc(I_GM, :);
WM_img = F_2D_psc(I_WM, :);
CSF_img = F_2D_psc(I_CSF, :);
all_img = [GM_img; WM_img; CSF_img]; 
% Identif limits between the different tissue compartments
line1_pos = numel(I_GM);
line2_pos = numel(I_GM) + numel(I_WM);

figure
fontsizeL = 17;
fontsizeM = 15;

% plot functional data
ax1 = subplot(5,1,[2:5]);
imagesc(ax1, all_img); 
colormap(gray); 
caxis(intensity_scale);
title(ax1, 'thePlotSpm','fontsize',fontsizeL)
ylabel(ax1, 'Voxels','fontsize',fontsizeM)
xlabel(ax1, 'fMRI volumes','fontsize',fontsizeM)
hold on; 
line([1 Nt],[line1_pos line1_pos],  'Color', 'b', 'LineWidth', 2 )
line([1 Nt],[line2_pos line2_pos],  'Color', 'r', 'LineWidth', 2 )
hold off;

% plot Framewise displacement
ax2 = subplot(5,1,1);
plot(ax2, FD, 'LineWidth', 2); grid;
axis tight
set(ax2,'Xticklabel',[]);
title(ax2, 'FD','fontsize',fontsizeL)
ylabel(ax2, 'mm','fontsize',fontsizeM)




