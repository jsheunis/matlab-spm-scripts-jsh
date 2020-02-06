% Data access
data_dir = '/Users/jheunis/Desktop/sample-data/sub-neufepmetest';
fn1 = [data_dir filesep 'dcm2niix_e1.nii'];
fn2 = [data_dir filesep 'dcm2niix_e2.nii'];
fn3 = [data_dir filesep 'dcm2niix_e3.nii'];

% Volume dimensions, and reference image
functional0_fn = [fn2 ',1'];
funcref_spm = spm_vol(functional0_fn);
funcref_3D  = spm_read_vols(funcref_spm);
[Nx, Ny, Nz] = size(funcref_3D);
Nt = 210;
Nvox = Nx*Ny*Nz;
N_skip = 0;
Ne = 3;
Ndyn = Nt - N_skip;
N_start = N_skip+1;
e_ref = 2;

% SPM motion correction settings
flagsSpmRealign = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'rtm',0,'PW','','lkp',1:6);
flagsSpmReslice = struct('quality',.9,'fwhm',5,'sep',4,...
    'interp',4,'wrap',[0 0 0],'mask',1,'mean',0,'which', 2);

% Get motion realignment template data and volume
dicomInfoVox   = sqrt(sum((funcref_spm.mat(1:3,1:3)).^2));
% fwhm = smoothing_kernel ./ dicomInfoVox;
A0=[];x1=[];x2=[];x3=[];wt=[];deg=[];b=[];
% R = struct;
R(1,1).mat = funcref_spm.mat;
R(1,1).dim = funcref_spm.dim;
R(1,1).Vol = funcref_3D;

% Initialize some variables
F = cell(Ne,1);
F_dyn_denoised = F;
rF = F;
realign_params = F;
fdyn_fn = F;
F_dyn_img = F;
currentVol = F;
MP = nan(Ndyn, 6);
% for e = 1:Ne
%     F{e} = nan(Nx*Ny*Nz,Ndyn);
%     F_dyn_denoised{e} = nan(Nx*Ny*Nz,Ndyn);
%     rF{e}  = nan(Nx*Ny*Nz,Ndyn);
% end

e = 2;


for i = 1:Nt
    % Display iteration number
    disp(num2str(i))
    % Skip iteration based on specified amount of initial volumes to skip
    if i <= N_skip
        continue;
    end
    j = i - N_skip;
    
    tic;
    % STEP 1: LOAD CURRENT VOLUME for reference echo   
    
    % Build dynamic filename
    fdyn_fn = [data_dir filesep 'dcm2niix_e' num2str(e) '.nii,' num2str(j)];
    currentVol = spm_vol(fdyn_fn);
    F_dyn_img = spm_read_vols(currentVol); % Read image volume with SPM. This is the raw/unprocessed image
    
    % (NOTE: SLICE TIME CORRECTION IS NOT DONE HERE, BUT OFTEN FORMS PART OF A STANDARD FMRI PROCESSING PIPELINE)
    
    % STEP 2 + 3: REALIGN AND RESLICE TO REFERENCE VOLUME
    % Current best practice is to realign one echo timeseries and then use
    % the resulting realignment paramaters to calculate and apply a single
    % affine transformation matrix to all echo volumes. Refer to tedana
    % docs
    R(2,1).mat = currentVol.mat;
    R(2,1).dim = currentVol.dim;
    R(2,1).Vol = F_dyn_img;

    % realign (FROM OPENNFT: preprVol.m)
    [R, A0, x1, x2, x3, wt, deg, b, nrIter] = spm_realign_rt_jsh(R, flagsSpmRealign, i, N_start, A0, x1, x2, x3, wt, deg, b);
    % get realignment parameters from affine matrices (this does the
    % opposite of spm_matrix, apparently) TODO: check if this does what it
    % is supposed to do
    tmpMCParam = spm_imatrix(R(2,1).mat / R(1,1).mat);
    if (i == N_start)
        offsetMCParam = tmpMCParam(1:6);
    end
    MP(j,:) = tmpMCParam(1:6) - offsetMCParam;
    
    % Reslice (FROM OPENNFT: preprVol.m)
    R(2,1).mat = currentVol.mat;
    R(2,1).dim = currentVol.dim;
    R(2,1).Vol = F_dyn_img;
    reslVol = spm_reslice_rt_jsh(R, flagsSpmReslice);
    
%     for e = 1:Ne
%         disp(num2str(e))
%         if e ~= e_ref
%             fdyn_fn{e} = [data_dir filesep 'dcm2niix_e' num2str(e) '.nii,' num2str(j)];
%             currentVol{e} = spm_vol(fdyn_fn{e});
%             Pm = zeros(12,1);
%             Pm(1:6) = MP(j, :);
%             orig_mat = currentVol{e}.mat;
%             rigid_mat = spm_matrix(Pm, 'T*R');
%             trans_mat = rigid_mat * orig_mat;
%             spm_get_space(fdyn_fn{e}, trans_mat);
%             currentVol{e} = spm_vol(fdyn_fn{e});
%             F_dyn_img{e} = spm_read_vols(currentVol{e});
%         end
%         R{e}(2,1).mat = currentVol{e}.mat;
%         R{e}(2,1).dim = currentVol{e}.dim;
%         R{e}(2,1).Vol = F_dyn_img{e};
%         reslVol{e} = spm_reslice_rt_jsh(R{e}, flagsSpmReslice);
% %         rF{e}(:,j) = realign_params{e}.reslVol(:);
%     end
    
    toc;
    
% %     spm_reslice
%     
%     
%     
%     
%     % STEP 4: MULTI-ECHO PARAMETER ESTIMATION AND COMBINATION (not for sample data)
%     if Ne > 1
%         % If option selected to use all echoes (i.e. use_echo = 0) continue
%         % with multi-echo parameter estimation and combination, else use
%         % specified echo (e.g use_echo = 2)
%         if use_echo == 0
%             % First estimate T2* using log-linear regression
%             S_pv = [rF{1}(I_mask,j)'; rF{2}(I_mask,j)'; rF{3}(I_mask,j)'];
%             S_pv = max(S_pv,1e-11);  %negative or zero signal values not allowed
%             b_pv = X\log(S_pv);
%             S0_pv(I_mask,j)=exp(b_pv(1,:));
%             T2star_pv(I_mask,j)=1./b_pv(2,:);
%             if isnan(b_pv)
%                 disp(['isnan: i = ' j ])
%             end
%             % Now threshold the T2star and S0 values based on expected (yet broad) range of T2star values
%             T2star_pv_corrected(I_mask,j) = T2star_pv(I_mask,j);
%             T2star_pv_corrected((T2star_pv_corrected(:,j)<0), j) = 0;
%             T2star_pv_corrected((T2star_pv_corrected(:,j)>T2star_thresh), j) = 0;
%             S0_pv_corrected(I_mask,j) = S0_pv(I_mask,j);
%             S0_pv_corrected((T2star_pv_corrected(:,j)<0), j) = 0;
%             S0_pv_corrected((T2star_pv_corrected(:,j)>T2star_thresh), j) = 0;
%             
%             % Combine the echoes using a specified weighting scheme (real-time
%             % T2star)
%             F_all_echoes = {rF{1}(:,j), rF{2}(:,j), rF{3}(:,j)};
%             weight_img = reshape(T2star_pv_corrected(:,j), Nx, Ny, Nz);
%             combined_t2s_rt(:, j) = combine_echoes_v2(TE, 1, weight_img, F_all_echoes, I_mask, [Nx, Ny, Nz]);
%             % rf = reshape(combined_t2s_rt(:, j), Nx, Ny, Nz);
%             rf = reshape(T2star_pv_corrected(:, j), Nx, Ny, Nz);
%         else
%             rf = rF{use_echo}(:,j);
%         end
%     else
%         % if single-echo, use first volume in rF cell array
%         rf = rF{1}(:,j);
%     end
%     
%     % STEP 5: SMOOTH REALIGNED VOLUME
%     % Using OpenNFT functionality and SPM
%     srf = zeros(Nx, Ny, Nz);
%     gKernel = smoothing_kernel ./ dicomInfoVox;
%     spm_smooth(rf, srf, gKernel);
%     srF(:,j) = srf(:);
    
end