function output = spm_standardPreproc_jsh(functional4D_fn, structural_fn, fwhm, spm_dir)
% Function to complete preprocessing of structural and functional data from
% a single subject for use in any other Matlab/SPM12 script.

% Steps include coregistering structural image to first functional image,
% segmenting the coregistered structural image into tissue types, and
% reslicing the segments to the functional resolution image grid. 
% Makes use of spm12 batch routines. If spm12 batch parameters are not
% explicitly set, defaults are assumed. 
%
% INPUT:
% funcional4D_fn     - filename of pre-real-time functional scan
% structural_fn      - filename of T1-weighted structural scan
% fwhm               - kernel size for smoothing operations
% spm_dir            - SPM12 directory location
% 
% OUTPUT: 
% output            - structure with filenames and data
%__________________________________________________________________________
% Copyright (C) Stephan Heunis 2018

% Load data
f4D_spm = spm_vol(functional4D_fn);
spm_size = size(f4D_spm);
Nt = spm_size(1);
% Declare output structure
output = struct;

% STEP 1 -- Realign (estimate and reslice) all functionals to first functional
disp('Step 1 -- Realign all volumes to first functional volume');
spm('defaults','fmri');
spm_jobman('initcfg');
realign_estimate_reslice = struct;
% Data
fnms={};
for i = 1:Nt
    fnms{i} = [functional4D_fn ',' num2str(i) ];
end
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.data={fnms'};
% Eoptions
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
% Roptions
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
% Run
spm_jobman('run',realign_estimate_reslice.matlabbatch);
[d, f, e] = fileparts(functional4D_fn);
output.rfunctional_fn = [d filesep 'r' f e];
output.mp_fn = [d filesep 'rp_' f '.txt'];
output.MP = load(output.mp_fn);
disp('Step 1 - Done!');

% STEP 2 -- Coregister structural image to first dynamic image (estimate only)
disp('Step 2 -- Coregister structural image to first dynamic image');
spm('defaults','fmri');
spm_jobman('initcfg');
coreg_estimate = struct;
% Ref
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[functional4D_fn ',1']};
% Source
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.source = {structural_fn};
% Eoptions
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
% Run
spm_jobman('run',coreg_estimate.matlabbatch);
disp('Step 2 - Done!');

% STEP 3 -- Segmentation of coregistered structural image into GM, WM, CSF, etc
% (with implicit warping to MNI space, saving forward and inverse transformations)
disp('Step 3 -- Segmentation');
spm('defaults','fmri');
spm_jobman('initcfg');
segmentation = struct;
% Channel
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
segmentation.matlabbatch{1}.spm.spatial.preproc.channel.vols = {structural_fn};
% Tissue
for t = 1:6
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = {[spm_dir filesep 'tpm' filesep 'TPM.nii,' num2str(t)]};
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = t-1;
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [1 0];
    segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
end
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
% Warp
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
segmentation.matlabbatch{1}.spm.spatial.preproc.warp.write=[1 1];
% Run
spm_jobman('run',segmentation.matlabbatch);
% Save filenames
[d, f, e] = fileparts(structural_fn);
output.forward_transformation = [d filesep 'y_' f e];
output.inverse_transformation = [d filesep 'iy_' f e];
output.gm_fn = [d filesep 'c1' f e];
output.wm_fn = [d filesep 'c2' f e];
output.csf_fn = [d filesep 'c3' f e];
output.bone_fn = [d filesep 'c4' f e];
output.soft_fn = [d filesep 'c5' f e];
output.air_fn = [d filesep 'c6' f e];
disp('Step 3 - done!');

% STEP 4 -- Reslice all to functional-resolution image grid
disp('Step 4 -- Reslice all to functional-resolution image grid');
spm('defaults','fmri');
spm_jobman('initcfg');
reslice = struct;
% Ref
reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {[functional4D_fn ',1']};
% Source
source_fns = {};
for i = 1:6
    source_fns{i} = [d filesep 'c' num2str(i) f e];
end
source_fns{7} = structural_fn;
reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_fns';
% Roptions
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
% Run
spm_jobman('run',reslice.matlabbatch);
% Save filenames
[d, f, e] = fileparts(structural_fn);
output.rstructural_fn = [d filesep 'r' f e];
output.rgm_fn = [d filesep 'rc1' f e];
output.rwm_fn = [d filesep 'rc2' f e];
output.rcsf_fn = [d filesep 'rc3' f e];
output.rbone_fn = [d filesep 'rc4' f e];
output.rsoft_fn = [d filesep 'rc5' f e];
output.rair_fn = [d filesep 'rc6' f e];
disp('Step 4 - done!');

% STEP 5 -- Gaussian kernel smoothing of realigned data
disp('STEP 5 -- Gaussian kernel smoothing of realigned data');
spm('defaults','fmri');
spm_jobman('initcfg');
smooth = struct;
% Data
fns={};
for i = 1:Nt
    fns{i} = [output.rfunctional_fn ',' num2str(i) ];
end
smooth.matlabbatch{1}.spm.spatial.smooth.data = fns';
% Other
smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's';
% Run
spm_jobman('run',smooth.matlabbatch);
[d, f, e] = fileparts(output.rfunctional_fn);
output.srfunctional_fn = [d filesep 's' f e];
disp('Step 5 - done!');

