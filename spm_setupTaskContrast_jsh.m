function spm_setupTaskContrast_jsh(stats_dir, params)

% SETUP BATCH JOB STRUCTURE
spm('defaults','fmri');
spm_jobman('initcfg');
contrast = struct;
% spmmat
contrast.matlabbatch{1}.spm.stats.con.spmmat = {[stats_dir filesep 'SPM.mat']};
% consess
contrast.matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = params.name;
contrast.matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = params.weights;
contrast.matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
% delete
contrast.matlabbatch{1}.spm.stats.con.delete = 0;
% RUN BATCH JOB
spm_jobman('run',contrast.matlabbatch);