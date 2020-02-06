function spm_specify1stlevel_jsh(stats_dir, func4D_fn, multi_reg_fn, params)

spm('defaults','fmri');
spm_jobman('initcfg');
design_stats = struct;
func4D_spm = spm_vol(func4D_fn);
func4D_size = size(func4D_spm);
Nt = func4D_size(1);

% SETUP BATCH JOB STRUCTURE
% dir
design_stats.matlabbatch{1}.spm.stats.fmri_spec.dir = {stats_dir}; 
% timing
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.units = params.timing_units;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.RT = params.timing_RT;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
% sess
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.scans = {};
for i = 1:Nt
    design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.scans{i,1} = [func4D_fn ',' num2str(i) ];
end
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.name = params.cond_name;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.onset = params.cond_onset;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.duration = params.cond_duration;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.tmod = 0;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.pmod = {''};
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.cond.orth = 1;
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.regress = {''};
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {multi_reg_fn};
design_stats.matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
% fact
design_stats.matlabbatch{1}.spm.stats.fmri_spec.fact = {''};
% bases
design_stats.matlabbatch{1}.spm.stats.fmri_spec.bases.hrf = struct('derivs', [0 0]);
% volt
design_stats.matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
% global
design_stats.matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
% mthresh
design_stats.matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8000;
% mask
design_stats.matlabbatch{1}.spm.stats.fmri_spec.mask = {''}; 
% cvi
design_stats.matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% RUN BATCH JOB
spm_jobman('run',design_stats.matlabbatch);