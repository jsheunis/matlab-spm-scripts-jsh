
% Data
data_dir = '/Users/jheunis/Desktop/sample-data/sub-neufepsetest';
test_rest_fn = [data_dir filesep 'test_rest.nii'];
rest_orig_fn = [data_dir filesep 'rest_orig.nii'];
rest_fn = [data_dir filesep 'rest.nii'];
% fn1 = [data_dir filesep 'dcm2niix_e1.nii'];
% fn2 = [data_dir filesep 'dcm2niix_e2.nii'];
% fn3 = [data_dir filesep 'dcm2niix_e3.nii'];

% Load 2nd echo data to get specs
rest_spm = spm_vol(rest_fn);
tmat = rest_spm.mat;
img_size = size(rest_spm);
rest_vols_1 = spm_read_vols(rest_spm(1));
Nt = img_size(1);
[Nx, Ny, Nz] = size(rest_vols_1);
Nvox = Nx*Ny*Nz;
%%
% Realign (estimate and reslice) all volumes of the 2nd echo timeseries to
% the first volume of the 2nd echo timeseries
spm('defaults','fmri');
spm_jobman('initcfg');
realign_estimate_reslice = struct;
% Data
fnms={};
for i = 1:Nt
    fnms{i} = [rest_fn ',' num2str(i) ];
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

%%
[d, f, e] = fileparts(rest_fn);
rrest_fn = [d filesep 'r' f e];
mp_fn = [d filesep 'rp_' f '.txt'];
MP = load(mp_fn);

rest_orig_fn = [data_dir filesep 'rest_orig.nii'];
rest_fn = [data_dir filesep 'rest.nii'];


rrest_spm = spm_vol(rrest_fn);

rest_spm = spm_vol(rest_fn);


rest_orig_spm = spm_vol(rest_orig_fn);


% rest = spm_read_vols(rest_spm);

rest_orig = spm_read_vols(rest_orig_spm);
% rrest = spm_read_vols(rrest_spm);




% test_rest_orig_spm = spm_vol(test_rest_orig_fn);


%%
figure;
plot(1:Nt, MP)
%%
test_rest_fn = [data_dir filesep 'test_rest.nii'];
test_rest_spm = spm_vol(test_rest_fn);
%%
    % FORMAT [A] = spm_matrix(P [,order])
% P(1)  - x translation
% P(2)  - y translation
% P(3)  - z translation
% P(4)  - x rotation about - {pitch} (radians)
% P(5)  - y rotation about - {roll}  (radians)
% P(6)  - z rotation about - {yaw}   (radians)


for i = 2:Nt
    P = zeros(12,1);
    P(1:6) = MP(i, :);
    orig_mat = rest_orig_spm(i).mat;
    rigid_mat = spm_matrix(P, 'T*R');
    trans_mat = rigid_mat * orig_mat;
    spm_get_space([test_rest_fn ',' num2str(i) ], trans_mat);
    spm_reslice
end


%%
new_ts_vol = spm_read_vols(spm_vol(test_rest_fn));
new_vol = new_ts_vol(:,:,:,1);
new_ts_spm = spm_vol(rest_fn);
new_spm = new_ts_spm(1);
new_spm.mat = trans_mat;
spm_write_vol(new_spm, new_vol);


% and change the name and the transformation matrix in the header
HeaderStructural_Translated.fname = 'Structural_Translated.img';
HeaderStructural_Translated.private.dat.fname = HeaderStructural_Translated.fname;
HeaderStructural_Translated.mat = TransformationMatrix_Translation;

spm_write_vol(HeaderStructural_Translated, Volume);


%%
test_rest_spm = spm_vol(test_rest_fn);



%% Guillaume
% You can determine the affine transformation matrix from a set of
% parameters using spm_matrix and you can apply a transformation to an
% image using spm_get_space. This is what spm_realign does:
%   https://github.com/spm/spm12/blob/r7487/spm_realign.m#L161
% There is a dedicated module in the batch interface to reorient images:
%   SPM > Util > Reorient Images
  

%% Remi
% Applying affine transformations to an image

HeaderStructural = spm_vol(FilesList.name);
TransformationMatrix = HeaderStructural.mat;
Volume = spm_read_vols(HeaderStructural);

% TRANSLATION
% Creates the vector needed for the spm_matrix function
P_translation = zeros(12,1);
P_translation(1) = 14;
P_translation(2) = -58;
P_translation(3) = 65;
% Creates the translation matrix
Translation = spm_matrix(P_translation,'T');
% Applies it to the transformation matrix
TransformationMatrix_Translation = Translation * TransformationMatrix;

% getting the header information of the original image
HeaderStructural_Translated = HeaderStructural;

% and change the name and the transformation matrix in the header
HeaderStructural_Translated.fname = 'Structural_Translated.img';
HeaderStructural_Translated.private.dat.fname = HeaderStructural_Translated.fname;
HeaderStructural_Translated.mat = TransformationMatrix_Translation;

spm_write_vol(HeaderStructural_Translated, Volume);


% ROTATION
% Creates the rotation matrix: spm works with radians, so if you want to 
% rotate the image by 90 degrees you have to transform it to radians first.
P_rotation = zeros(12,1);
P_rotation(4) = 90*pi/180; 
Rotation = spm_matrix(P_rotation,'R');

TransformationMatrix_Rotation = Rotation * TransformationMatrix;

% getting the header information of the original image
HeaderStructural_Rotated = HeaderStructural;

% and change the name and the transformation matrix in the header
HeaderStructural_Rotated.fname = 'Structural_Rotated.img';
HeaderStructural_Rotated.private.dat.fname = HeaderStructural_Rotated.fname;
HeaderStructural_Rotated.mat = TransformationMatrix_Rotation;

spm_write_vol(HeaderStructural_Rotated, Volume);


cd(StartFolder)

