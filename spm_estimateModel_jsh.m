function spm_estimateModel_jsh(stats_dir)

% SETUP BATCH JOB STRUCTURE
spm('defaults','fmri');
spm_jobman('initcfg');
model_estimation = struct;
% spmmat
model_estimation.matlabbatch{1}.spm.stats.fmri_est.spmmat = {[stats_dir filesep 'SPM.mat']};
% write_residuals
model_estimation.matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
% method
model_estimation.matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
% RUN BATCH JOB
spm_jobman('run',model_estimation.matlabbatch);
