function spm_runResults_jsh(stats_dir)

spm('defaults','fmri');
spm_jobman('initcfg');
% SETUP BATCH JOB STRUCTURE
results = struct;
% spmmat
results.matlabbatch{1}.spm.stats.results.spmmat = {[stats_dir filesep 'SPM.mat']};
% conspec
results.matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
results.matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
results.matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
results.matlabbatch{1}.spm.stats.results.conspec.thresh = 0.0500;
results.matlabbatch{1}.spm.stats.results.conspec.extent = 0;
results.matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
results.matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
% units
results.matlabbatch{1}.spm.stats.results.units = 1;
% export
results.matlabbatch{1}.spm.stats.results.export{1}.ps = 1;
% RUN BATCH JOB
spm_jobman('run',results.matlabbatch);