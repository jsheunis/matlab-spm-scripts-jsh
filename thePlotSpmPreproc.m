function output = thePlotSpmPreproc(functional4D_fn, structural_fn, fwhm, spm_dir)
% Function to complete preprocessing of structural and functional data from
% a single subject for use in thePlotSpm.m
%
% Steps include coregistering structural image to first functional image,
% segmenting the coregistered structural image into tissue types, and
% reslicing the segments to the functional resolution image grid.
%
% Makes use of spm12 batch routines.
% If spm12 batch parameters are not explicitly set, defaults are assumed.
%
% INPUT:
% funcional4D_fn     - filename of pre-real-time functional scan
% structural_fn      - filename of T1-weighted structural scan
% fwhm               - kernel size for smoothing operations
%
% OUTPUT:
% output            - structure with filenames and data
%
% EXAMPLE
% fwhm = 6;
% [spm_dir, ~, ~] = fileparts(which('spm'));
% data_dir = 'D:';
% subj = 'sub-01';
% % structural image
% structural_fn = spm_select('FPList', fullfile(data_dir , subj , 'anat'), ['^' subj '_T1w.*nii$']);
% % 4D BOLD time series
% functional4D_fn = spm_select('FPList', fullfile(data_dir , subj , 'func'), ['^' subj '_task-.*_bold*.nii$']);
% _________________________________________________________________________
% Copyright (C) Stephan Heunis

%% Check that we have everything
if nargin<4
    % trying to locate spm
    [spm_dir, ~, ~] = fileparts(which('spm'));
    if isempty(spm_dir)
        error('Tell me where SPM is, please.')
    end
end

if nargin<3
    fwhm = 6;
end

if nargin<2
    error('I need data to work with.')
end

%  check number of volume in this 4D time series
Nt = numel(spm_vol(functional4D_fn));

% Declare output structure
output = struct;

%% STEP 1 -- Realign (estimate and reslice) all functionals to first functional
% Only if this has not been done before
[dir, file, ext] = fileparts(functional4D_fn);
if ~exist(spm_select('FPList', dir, ['^r' file ext]), 'file')
    
    disp('Step 1...');
    realign_estimate_reslice = struct;
    % Data
    fnms={};
    for i = 1:Nt
        fnms{i} = [functional4D_fn ',' num2str(i) ];
    end
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.data={fnms'};
    % Estimate options
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    % Reslice options
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    realign_estimate_reslice.matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    % Run
    cfg_util('run',realign_estimate_reslice.matlabbatch);
    
end

% identify the outputed realigned file
output.rfunctional_fn = fullfile(dir , ['r' file ext]);
% load the data from the realignement parameters
output.mp_fn = fullfile(dir , ['rp_' file '.txt']);
output.MP = load(output.mp_fn);
disp('Step 1 - Done!');


%% STEP 2 -- Coregister structural image to first dynamic image (estimate only)
disp('Step 2...');
spm('defaults','fmri');
spm_jobman('initcfg');
coreg_estimate = struct;
% Ref
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[functional4D_fn ',1']};
% Source
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.source = {structural_fn};
% Estimate options
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
coreg_estimate.matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
% Run
spm_jobman('run',coreg_estimate.matlabbatch);
disp('Step 2 - Done!');


%% STEP 3 -- Segmentation of coregistered structural image into GM, WM, CSF, etc
% (with implicit warping to MNI space, saving forward and inverse transformations)

% Only if this has not been done before
[dir, file, ext] = fileparts(structural_fn);
if ~exist(spm_select('FPList', dir, ['^c1' file ext]), 'file')
    
    disp('Step 3...');
    segmentation = struct;
    % Channel
    segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
    segmentation.matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    segmentation.matlabbatch{1}.spm.spatial.preproc.channel.write = [0 1];
    segmentation.matlabbatch{1}.spm.spatial.preproc.channel.vols = {structural_fn};
    % Tissue
    for t = 1:6
        segmentation.matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm = {fullfile(spm_dir , 'tpm' , ['TPM.nii,' num2str(t)])};
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
    cfg_util('run',segmentation.matlabbatch);
    
end

% Save filenames
output.forward_transformation = fullfile(dir , ['y_' file ext]);
output.inverse_transformation = fullfile(dir , ['iy_' file ext]);
output.gm_fn = fullfile(dir , ['c1' file ext]);
output.wm_fn = fullfile(dir , ['c2' file ext]);
output.csf_fn = fullfile(dir , ['c3' file ext]);
output.bone_fn = fullfile(dir , ['c4' file ext]);
output.soft_fn = fullfile(dir , ['c5' file ext]);
output.air_fn = fullfile(dir , ['c6' file ext]);
disp('Step 3 - Done!');


%% STEP 4 -- Reslice all to functional-resolution image grid
% Only if this has not been done before
[dir, file, ext] = fileparts(structural_fn);
if ~exist(spm_select('FPList', dir, ['^rc1' file ext]), 'file')
    
    disp('Step 4...');
    reslice = struct;
    % Ref
    reslice.matlabbatch{1}.spm.spatial.coreg.write.ref = {[functional4D_fn ',1']};
    % Source: images to reslice
    source_fns = cellstr(spm_select('FPList', dir, ['^c'  '.*nii$']))'; % all the tissue probability maps of that subject
    source_fns{end+1} = structural_fn;
    reslice.matlabbatch{1}.spm.spatial.coreg.write.source = source_fns';
    % Roptions
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    reslice.matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
    % Run
    cfg_util('run',reslice.matlabbatch);
    
end

% Save filenames
[dir, file, ext] = fileparts(structural_fn);
output.rstructural_fn = fullfile(dir , ['r' file ext]);
output.rgm_fn = fullfile(dir , ['rc1' file ext]);
output.rwm_fn = fullfile(dir , ['rc2' file ext]);
output.rcsf_fn = fullfile(dir , ['rc3' file ext]);
output.rbone_fn = fullfile(dir , ['rc4' file ext]);
output.rsoft_fn = fullfile(dir , ['rc5' file ext]);
output.rair_fn = fullfile(dir , ['rc6' file ext]);
disp('Step 4 - Done!');


%% STEP 5 -- Gaussian kernel smoothing of realigned data
% Only if this has not been done before
[dir, file, ext] = fileparts(functional4D_fn);
if ~exist(spm_select('FPList', dir, ['^sr' file ext]), 'file')
    
    disp('Step 5...');
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
    cfg_util('run',smooth.matlabbatch);
    
end

% Save filenames
[dir, file, ext] = fileparts(output.rfunctional_fn);
output.srfunctional_fn = fullfile(dir , ['s' file ext]);
disp('Step 5 - Done!');


%% STEP 6 -- Gaussian kernel smoothing of unprocessed data
% Only if this has not been done before
[dir, file, ext] = fileparts(functional4D_fn);
if ~exist(spm_select('FPList', dir, ['^s' file ext]), 'file')
    
    disp('Step 6...');
    smooth = struct;
    % Data
    fns={};
    for i = 1:Nt
        fns{i} = [functional4D_fn ',' num2str(i) ];
    end
    smooth.matlabbatch{1}.spm.spatial.smooth.data = fns';
    % Other
    smooth.matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm fwhm fwhm];
    smooth.matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    smooth.matlabbatch{1}.spm.spatial.smooth.im = 0;
    smooth.matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    % Run
    cfg_util('run',smooth.matlabbatch);
    
end

% Save filenames
output.sfunctional_fn = fullfile(dir , ['s' file ext]);
disp('Step 6 - Done!');