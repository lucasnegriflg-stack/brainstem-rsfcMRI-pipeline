function Denoising_LC_paperready_batch_v9_4_rankrobust_PATCHED_STABLE()
%% ==========================================================
%% Denoising_LC (BATCH, Subject-Matching) — v9.4 PAPER-READY (rank-robust, consistent)
%% + STABILITY-HARDENED PATCH SET (R2024b-friendly)
%%
%% v9.4 LC patches (2026-02) preserved:
%%   (1) T1-primary mask: brain_mask = T1mask & valid (no MeanEPI AND)
%%   (2) StrictCSF: no early zscore; normalize only after bandpass with other continuous regressors
%%   (3) GSR: compute global signal consistently in filtered space from Y_bp (brain_idx)
%%   (4) Output NIfTI: delete existing file before writing + set pinfo to avoid scaling ambiguity
%%
%% STABILITY PATCHES (Reviewer-safe engineering, no methodological upheaval):
%%   (A) STREAMING NIfTI writer: no huge 4D RAM allocation
%%   (B) Deterministic aCompCor: PCA sign-consistency (prevents sign flips across runs)
%%   (C) Robust regression solver: QR with column pivoting (rank-aware, stable)
%%   (D) Reslice hygiene: delete old r*.nii before spm_reslice (avoids stale files across reruns)
%%   (E) Filter sanitizing: replace non-finite values before filtfilt (hardens pipeline)
%%   (F) Optional strict reproducibility mode: rng + single-thread
%%   (G) Fail-fast mask sanity threshold (optional)
%% ==========================================================

clear; clc;
spm('defaults','fmri');
spm_jobman('initcfg');

fprintf('\n=== START: LC-RSFC Denoising (v9.4 + STABILITY PATCHES, native space, batch) ===\n\n');

%% -------------------- TOOLBOX SANITY --------------------
if ~(exist('butter','file')==2 && exist('filtfilt','file')==2)
    error('Signal Processing Toolbox required for butter/filtfilt bandpass.');
end

%% -------------------- CONFIG --------------------
cfg.TR            = 2.21;          % FIXED TR (seconds)
cfg.band          = [0.008 0.09];  % rsfMRI band
cfg.filter_order  = 2;             % Butterworth order

cfg.doGSR_sens    = true;          % ON by default for LC paper sensitivity

cfg.WM_thr        = 0.90;          % WM prob threshold
cfg.CSF_thr       = 0.95;          % CSF prob threshold

cfg.FD_thr        = 0.5;          % scrubbing FD threshold (mm)
cfg.DVARSz_thr    = 5.00;          % scrubbing DVARS z threshold (robust)

cfg.scrub_warnPct = 20;            % warn if >20% scrubbed
cfg.meanFD_fail   = 0.50;          % QC_fail if meanFD > 0.5
cfg.scrub_failPct = 20;            % QC_fail if scrubPct > 20

cfg.writeFloat32  = true;          % output dtype single/double
cfg.writePrefix   = 'dw_';         % output prefix no-GSR
cfg.writePrefixG  = 'dwGSR_';      % output prefix with GSR
cfg.writeMaskPrefix = 'ANALYSISMASK_';

% brainmask parameters (MeanEPI mask)
cfg.brain_thr_prct = 20;           % threshold at Xth percentile of meanEPI within valid candidates
cfg.brain_min_vox  = 1000;         % fallback to valid if too few candidates

% optional intersection with T1 brainmask built from (c1+c2+c3)>thr, resliced NN to EPI
cfg.intersect_with_T1brain = true;
cfg.t1_brain_thr           = 0.20;

% Use T1 brainmask as PRIMARY analysis mask if c1 exists (LC-robust)
cfg.use_T1brain_as_primary = true;

% LC-proofing: optional erosion of tissue masks after reslicing to EPI grid
cfg.doErodeTissues = true;
cfg.erodeIterWM    = 1;
cfg.erodeIterCSF   = 0;

% strict CSF mean confound
cfg.addStrictCSFmean = true;
cfg.CSF_strict_thr   = 0.98;

% optional: save bandpassed design matrix actually used in regression
cfg.saveRegressorsBP = true;

% --- PAPER-READY additions ---
cfg.addIntercept       = true;
cfg.saveBadVols        = true;
cfg.markBadVolsInQC    = true;
cfg.saveCfgSnapshot    = true;
cfg.saveDesignMatMAT   = true;

% ===== STABILITY OPTIONS =====
cfg.regressionSolver   = 'qrp';     % 'qrp' (recommended) | 'lsqminnorm' | 'pinv'
cfg.strictRepro        = false;     % if true: rng + single-thread for extra determinism
cfg.failIfMaskTooSmall = true;      % fail fast if analysis mask is implausibly small
cfg.minMaskVox         = 5000;      % lower bound (tune to your EPI; 5k is conservative)
cfg.sanitizeFinite     = true;      % replace NaN/Inf with 0 before filtfilt
cfg.pcaSignFix         = true;      % enforce consistent PCA component signs
cfg.resliceDeleteOldR  = true;      % delete existing r*.nii before reslicing (avoid stale files)

if cfg.strictRepro
    rng(0,'twister');
    try
        maxNumCompThreads(1);
    catch
    end
    fprintf(' Strict reproducibility mode: rng(0) + single-thread (where supported)\n');
end

%% -------------------- INPUTS (MultiSelect) --------------------
fprintf('Select required files (multi-select allowed):\n');

% Optional T1
[fn_T1,fp_T1] = uigetfile({'*.nii','NIfTI (*.nii)'}, ...
    'Select optional T1 (Cancel = skip)', 'MultiSelect','on');
fn_T1 = normalize_multiselect(fn_T1, true);

% Optional GM segmentation (c1)
[fn_C1,fp_C1] = uigetfile({'c1*.nii','Gray matter (c1*)'}, ...
    'Select optional GM masks (c1*, T1 space) (Cancel = skip)', 'MultiSelect','on');
fn_C1 = normalize_multiselect(fn_C1, true);

% WM / CSF segmentations — required
[fn_C2,fp_C2] = uigetfile({'c2*.nii','White matter (c2*)'}, ...
    'Select WM masks (c2*, T1 space)', 'MultiSelect','on');
if isequal(fn_C2,0), error('Aborted: WM (c2*) is required.'); end
fn_C2 = normalize_multiselect(fn_C2);

[fn_C3,fp_C3] = uigetfile({'c3*.nii','CSF (c3*)'}, ...
    'Select CSF masks (c3*, T1 space)', 'MultiSelect','on');
if isequal(fn_C3,0), error('Aborted: CSF (c3*) is required.'); end
fn_C3 = normalize_multiselect(fn_C3);

% fMRI 4D EPI (native)
[fn_EPI,fp_EPI] = uigetfile({'*.nii','4D EPI (*.nii)'}, ...
    'Select native 4D EPI (e.g., cT1_*.nii)', 'MultiSelect','on');
if isequal(fn_EPI,0), error('Aborted: EPI is required.'); end
fn_EPI = normalize_multiselect(fn_EPI);

% motion parameters (rp_*.txt)
[fn_Mot,fp_Mot] = uigetfile({'rp_*.txt','Motion params (rp_*.txt)'}, ...
    'Select motion parameters (rp_*.txt)', 'MultiSelect','on');
if isequal(fn_Mot,0), error('Aborted: motion rp_*.txt is required.'); end
fn_Mot = normalize_multiselect(fn_Mot);

% optional LC masks (ROI tSNR only)
[fn_LCL,fp_LCL] = uigetfile({'*.nii','LC left mask (*.nii)'}, ...
    'Select optional LC_L masks (Cancel = skip)', 'MultiSelect','on');
fn_LCL = normalize_multiselect(fn_LCL, true);

[fn_LCR,fp_LCR] = uigetfile({'*.nii','LC right mask (*.nii)'}, ...
    'Select optional LC_R masks (Cancel = skip)', 'MultiSelect','on');
fn_LCR = normalize_multiselect(fn_LCR, true);

%% -------------------- BUILD LISTS (SID matching) --------------------
fprintf('\nBuilding subject lists from filenames...\n');

T1_list   = make_list_from_files(fn_T1, fp_T1);
C1_list   = make_list_from_files(fn_C1, fp_C1);
C2_list   = make_list_from_files(fn_C2, fp_C2);
C3_list   = make_list_from_files(fn_C3, fp_C3);
EPI_list  = make_list_from_files(fn_EPI, fp_EPI);
Mot_list  = make_list_from_files(fn_Mot, fp_Mot);
LC_L_list = make_list_from_files(fn_LCL, fp_LCL);
LC_R_list = make_list_from_files(fn_LCR, fp_LCR);

if isempty(C2_list) || isempty(C3_list) || isempty(EPI_list) || isempty(Mot_list)
    error('After Subject-ID parsing, a required list is empty (C2/C3/EPI/Motion).');
end

common_sid = intersect(intersect({C2_list.sid},{C3_list.sid}), intersect({EPI_list.sid},{Mot_list.sid}));
if isempty(common_sid)
    error('No common Subject-IDs found across required modalities (C2/C3/EPI/Motion).');
end
fprintf('\n Matched subjects (required):\n');
disp(common_sid(:));

%% -------------------- PROCESS LOOP --------------------
for s = 1:numel(common_sid)
    sid = common_sid{s};
    fprintf('\n=============================================\n');
    fprintf(' Subject: %s\n', sid);
    fprintf('=============================================\n');

    % required
    C2_in      = get_fname_for_sid(C2_list,  sid);
    C3_in      = get_fname_for_sid(C3_list,  sid);
    EPI        = get_fname_for_sid(EPI_list, sid);
    MotionFile = get_fname_for_sid(Mot_list, sid);

    % optional
    T1    = get_fname_for_sid_or_empty(T1_list,    sid);
    C1_in = get_fname_for_sid_or_empty(C1_list,    sid);
    LC_L  = get_fname_for_sid_or_empty(LC_L_list,  sid);
    LC_R  = get_fname_for_sid_or_empty(LC_R_list,  sid);

    outdir = fileparts(EPI);

    fprintf(['\nFiles:\n T1: %s\n c1(GM): %s\n WM(c2): %s\n CSF(c3): %s\n EPI: %s\n Motion: %s\n' ...
             ' LC_L: %s\n LC_R: %s\n\n'], ...
        tern(~isempty(T1), T1, '(not selected)'), ...
        tern(~isempty(C1_in), C1_in, '(not selected)'), ...
        C2_in, C3_in, EPI, MotionFile, ...
        tern(~isempty(LC_L), LC_L, '(not selected)'), ...
        tern(~isempty(LC_R), LC_R, '(not selected)'));

    %% ---------- Versions ----------
    matlabVersion = version;
    spmVersion = 'unknown';
    try, spmVersion = spm('Ver'); catch, end

    %% ---------- Load EPI header ----------
    V_EPI_all = spm_vol(EPI);
    nt = numel(V_EPI_all);
    if nt < 60
        warning('Short timeseries (nt=%d). Bandpass stability may be limited (%s).', nt, sid);
    end
    V_EPI_mid = V_EPI_all(round(nt/2));

    %% ---------- TR fixed ----------
    TR = cfg.TR;
    fs = 1/TR;
    fprintf('TR (%s): %.4f s (fs=%.4f Hz)\n', sid, TR, fs);

    %% ---------- Bandpass sanity check ----------
    Wn = cfg.band/(fs/2);
    if any(Wn <= 0) || any(Wn >= 1) || Wn(1) >= Wn(2)
        error('Invalid bandpass. band=[%.4f %.4f], fs=%.4f -> Wn=[%.4f %.4f] (%s)', ...
            cfg.band(1), cfg.band(2), fs, Wn(1), Wn(2), sid);
    end

    %% ---------- Ensure WM/CSF masks in EPI grid (NN) ----------
    [C2_epi, C3_epi] = reslice_tissue_masks_to_epi(C2_in, C3_in, V_EPI_mid, sid, cfg);

    %% ---------- Load 4D EPI ----------
    Y4 = spm_read_vols(V_EPI_all); % [nx,ny,nz,nt]
    [nx,ny,nz,nt_check] = size(Y4);
    if nt_check ~= nt
        nt = nt_check;
        fprintf(' Updated nt from data: %d\n', nt);
    end

    %% ---------- Pre-denoise tSNR ----------
    mu = mean(Y4,4);
    sd = std(Y4,0,4);
    sd(sd==0) = NaN;
    tSNR_map = mu ./ sd;

    % valid voxels
    Y2D_raw = reshape(Y4, [], nt)';          % [nt x Nvox]
    clear Y4;

    sdY = std(Y2D_raw,0,1);
    nonzero = any(Y2D_raw~=0,1);
    valid = (sdY > 0) & nonzero;             % [1 x Nvox]

    valid_tsnr = isfinite(tSNR_map(:)) & isfinite(sd(:)) & nonzero(:);
    mean_tSNR = mean(tSNR_map(valid_tsnr),'omitnan');
    fprintf('Global tSNR pre-denoise (%s): %.2f\n', sid, mean_tSNR);

    %% ---------- Build ANALYSIS BRAIN MASK ----------
    [~, mean_mask, mean_lcc] = build_brainmask_from_meanEPI(mu, valid, [nx,ny,nz], cfg, sid);

    t1_intersection_done = false;
    t1_primary_used      = false;
    lcc_applied          = mean_lcc;

    brain_mask = mean_mask(:) & valid(:); % fallback default

    if cfg.intersect_with_T1brain && ~isempty(C1_in) && isfile(C1_in)
        t1mask_vec = get_T1_brainmask_in_epi_grid(C1_in, C2_in, C3_in, V_EPI_mid, [nx,ny,nz], cfg, sid, outdir);

        if ~isempty(t1mask_vec)
            t1_intersection_done = true;

            if isfield(cfg,'use_T1brain_as_primary') && cfg.use_T1brain_as_primary
                % v9.4 PATCH (1): truly primary = T1mask & valid (no MeanEPI AND)
                brain_mask = t1mask_vec(:) & valid(:);
                t1_primary_used = true;

                % Optional: LCC on T1-primary mask
                lcc_applied = false;
                try
                    brain3 = reshape(brain_mask, [nx,ny,nz]);
                    CC = bwconncomp(brain3, 6);
                    if CC.NumObjects > 0
                        [~,imax] = max(cellfun(@numel, CC.PixelIdxList));
                        brain3 = false(nx,ny,nz);
                        brain3(CC.PixelIdxList{imax}) = true;
                        brain_mask = brain3(:);
                        lcc_applied = true;
                    end
                catch
                end
            else
                brain_mask = mean_mask(:) & t1mask_vec(:) & valid(:);
                lcc_applied = mean_lcc;
            end
        end
    end

    brain_idx = find(brain_mask);
    if isempty(brain_idx)
        warning('Analysis mask empty; falling back to MeanEPI mask & valid (%s).', sid);
        brain_mask = mean_mask(:) & valid(:);
        brain_idx = find(brain_mask);
        t1_primary_used = false;
    end

    nMaskVox = numel(brain_idx);
    fprintf(' Analysis brain mask voxels (%s): %d | T1_intersect=%s | T1_primary=%s | LCC=%s\n', ...
        sid, nMaskVox, tern(t1_intersection_done,'ON','OFF'), tern(t1_primary_used,'ON','OFF'), tern(lcc_applied,'ON','OFF'));

    if cfg.failIfMaskTooSmall && (nMaskVox < cfg.minMaskVox)
        error('Analysis mask too small (%d vox) for subject %s. Aborting (fail-fast).', nMaskVox, sid);
    end

    %% ---------- Optional LC ROI tSNR ----------
    LC_L_tSNR = NaN;
    if ~isempty(LC_L) && isfile(LC_L), LC_L_tSNR = compute_roi_tsnr(LC_L, V_EPI_mid, tSNR_map, sid, 'LC_L'); end
    LC_R_tSNR = NaN;
    if ~isempty(LC_R) && isfile(LC_R), LC_R_tSNR = compute_roi_tsnr(LC_R, V_EPI_mid, tSNR_map, sid, 'LC_R'); end

    %% ---------- Write explicit analysis mask ----------
    Vmask = V_EPI_all(1);
    Vmask.fname = fullfile(outdir, [cfg.writeMaskPrefix spm_file(EPI,'basename') '.nii']);
    Vmask.dt = [2 0];         % uint8
    Vmask.pinfo = [1;0;0];
    spm_write_vol(Vmask, uint8(reshape(brain_mask, [nx,ny,nz])));
    fprintf(' Analysis mask written (%s): %s\n', sid, Vmask.fname);

    %% ---------- Tissue masks + optional erosion ----------
    WMvol  = spm_read_vols(spm_vol(C2_epi)) > cfg.WM_thr;
    CSFvol = spm_read_vols(spm_vol(C3_epi)) > cfg.CSF_thr;

    if cfg.doErodeTissues
        [WMvol, CSFvol] = erode_tissues_if_possible(WMvol, CSFvol, cfg, sid);
    end

    WM_mask  = WMvol(:)'  & valid;    % [1 x Nvox]
    CSF_mask = CSFvol(:)' & valid;

    %% ---------- Strict CSF mean confound (no early zscore) ----------
    StrictCSF = [];
    nStrictCSFvox = 0;

    if isfield(cfg,'addStrictCSFmean') && cfg.addStrictCSFmean
        try
            CSFprob = spm_read_vols(spm_vol(C3_epi));  % probabilities, EPI grid
            strict_mask = (CSFprob > cfg.CSF_strict_thr);
            strict_vec  = strict_mask(:)' & valid;

            nStrictCSFvox = nnz(strict_vec);
            if nStrictCSFvox > 10
                StrictCSF = mean(Y2D_raw(:, strict_vec), 2, 'omitnan'); % raw scale here
                fprintf(' Strict CSF mean added (%s): thr=%.2f, vox=%d\n', sid, cfg.CSF_strict_thr, nStrictCSFvox);
            else
                fprintf(' Strict CSF mean skipped (%s): too few voxels (%d)\n', sid, nStrictCSFvox);
                StrictCSF = [];
            end
        catch ME
            warning('Strict CSF mean failed (%s): %s', sid, ME.message);
            StrictCSF = [];
            nStrictCSFvox = 0;
        end
    end

    %% ---------- aCompCor (deterministic sign) ----------
    [WM_comp, kW, varWM]   = acompcor_from_mask(Y2D_raw, WM_mask, nt, sid, 'WM', cfg);
    [CSF_comp,kC, varCSF]  = acompcor_from_mask(Y2D_raw, CSF_mask, nt, sid, 'CSF', cfg);
    CompCor = zscore([WM_comp CSF_comp],0,1);

    fprintf(' aCompCor (%s): WM %d PCs (%.1f%%), CSF %d PCs (%.1f%%)\n', ...
        sid, kW, varWM, kC, varCSF);

    %% ---------- Motion + Scrubbing ----------
    Motion = load(MotionFile);
    Motion = Motion(1:min(size(Motion,1),nt),:);
    if size(Motion,1) < nt
        Motion(end+1:nt,:) = repmat(Motion(end,:), nt-size(Motion,1), 1);
        warning('Motion shorter than fMRI; padded last row (%s).', sid);
    end

    MotionReg = friston36(Motion);
    fprintf(' Motion regressors (%s): %d columns (Friston-36)\n', sid, size(MotionReg,2));

    [FD, meanFD, medianFD, maxFD] = compute_FD_power(Motion);

    % DVARS z on ANALYSIS MASK
    DV_z = compute_DVARSz(Y2D_raw(:,brain_idx));
    % --- DVARS summary statistics (robust z) ---
    meanDVARSz   = mean(DV_z, 'omitnan');
    medianDVARSz = median(DV_z, 'omitnan');
    maxDVARSz    = max(DV_z);
    p95DVARSz    = prctile(DV_z, 95);


    badVols = find(FD > cfg.FD_thr | DV_z > cfg.DVARSz_thr);
    nScrub = numel(badVols);
    scrubPerc = 100*nScrub/nt;

    badVolsFile = '';
    if cfg.saveBadVols
        badVolsFile = fullfile(outdir, ['BadVols_' spm_file(EPI,'basename') '.txt']);
        fidbv = fopen(badVolsFile,'w');
        if fidbv ~= -1
            fprintf(fidbv, '%d\n', badVols);
            fclose(fidbv);
        else
            warning('Could not write badVols file (%s).', sid);
        end
    end

    if scrubPerc > cfg.scrub_warnPct
        warning(' >%.0f%% volumes flagged for scrubbing (%.1f%%, %s)', cfg.scrub_warnPct, scrubPerc, sid);
    end
    fprintf('Scrubbing candidates (%s): %d/%d (%.1f%%)\n', sid, nScrub, nt, scrubPerc);

    %% ---------- Design matrix pieces ----------
    X_motion = zscore(MotionReg,0,1);
    X_comp   = zscore(CompCor,0,1);

    if ~isempty(StrictCSF)
        X_cont = [X_motion X_comp StrictCSF]; % StrictCSF unscaled here
    else
        X_cont = [X_motion X_comp];
    end

    scrubMat = [];
    if ~isempty(badVols)
        scrubMat = zeros(nt, nScrub);
        for i=1:nScrub
            scrubMat(badVols(i),i) = 1;
        end
    end

    Intercept = [];
    if cfg.addIntercept
        Intercept = ones(nt,1);
    end

    %% ---------- Save regressors CSV (UNFILTERED base) ----------
    csv_regs = fullfile(outdir, ['Regressors_' spm_file(EPI,'basename') '.csv']);
    X_base_unf = [X_cont scrubMat Intercept];

    if ~isempty(X_base_unf)
        sdx = std(X_base_unf,0,1);
        keep = (sdx > 0);
    else
        keep = [];
    end
    if cfg.addIntercept && ~isempty(keep)
        keep(end) = true;
    end
    X_base_unf = X_base_unf(:, keep);

    writematrix(X_base_unf, csv_regs);
    fprintf(' Regressors CSV (no GSR, unfiltered) saved (%s): %s\n', sid, csv_regs);

    %% ---------- Filtering (data + continuous regressors only) ----------
    [b,a] = butter(cfg.filter_order, Wn);

    % Data: detrend linear then bandpass
    Y_dt = detrend(Y2D_raw, 'linear');
    if cfg.sanitizeFinite
        Y_dt(~isfinite(Y_dt)) = 0;
    end
    Y_bp = filtfilt(b,a,double(Y_dt));

    % Continuous regressors: demean then bandpass
    Xc_dt = detrend(X_cont, 'constant');
    if cfg.sanitizeFinite
        Xc_dt(~isfinite(Xc_dt)) = 0;
    end
    Xc_bp = filtfilt(b,a,double(Xc_dt));
    Xc_bp = zscore(Xc_bp,0,1); % includes StrictCSF normalization here

    % Final design (no GSR): bandpassed continuous + unfiltered spikes + unfiltered intercept
    X1 = [Xc_bp scrubMat Intercept];

    sdx1 = std(X1,0,1);
    keep1 = (sdx1 > 0);
    if cfg.addIntercept && ~isempty(keep1)
        keep1(end) = true;
    end
    X1 = X1(:, keep1);

    rX1 = rank(X1);
    nReg_noGSR = size(X1,2);
    dof_noGSR = nt - rX1;

    if rX1 < size(X1,2)
        fprintf(' Design matrix (no GSR) not full rank (rank=%d/%d, %s)\n', rX1, size(X1,2), sid);
    end

    if cfg.saveRegressorsBP
        csv_regs_bp = fullfile(outdir, ['RegressorsBP_' spm_file(EPI,'basename') '.csv']);
        writematrix(X1, csv_regs_bp);
        fprintf(' RegressorsBP CSV (design used, no GSR) saved (%s): %s\n', sid, csv_regs_bp);
    end

    if cfg.saveDesignMatMAT
        matX1 = fullfile(outdir, ['DesignMatrix_noGSR_' spm_file(EPI,'basename') '.mat']);
        save(matX1, 'X1', 'badVols', 'nt', 'TR', 'sid');
    end

    %% ---------- Regression: No-GSR pipeline ----------
    Yclean1 = regress_out_rankrobust(Y_bp, X1, cfg.writeFloat32, cfg.regressionSolver);

    out_nii1 = fullfile(outdir, [cfg.writePrefix spm_file(EPI,'filename')]);
    if exist(out_nii1,'file'), delete(out_nii1); end
    write_4d_nii_from_2d(Yclean1, V_EPI_all, out_nii1, cfg.writeFloat32);
    fprintf(' Denoised EPI saved (no GSR) (%s): %s\n', sid, out_nii1);

    %% ---------- (Optional) GSR sensitivity pipeline ----------
    rX2 = NaN; nReg_GSR = NaN; dof_GSR = NaN;
    if cfg.doGSR_sens
        % v9.4 PATCH (3): compute global signal in filtered space (analysis mask)
        Gsig_bp = mean(Y_bp(:, brain_idx), 2, 'omitnan');
        Gsig_bp = zscore(Gsig_bp,0,1);

        X2 = [Xc_bp scrubMat Gsig_bp Intercept];

        sdx2 = std(X2,0,1);
        keep2 = (sdx2 > 0);
        if cfg.addIntercept && ~isempty(keep2)
            keep2(end) = true;
        end
        X2 = X2(:, keep2);

        rX2 = rank(X2);
        nReg_GSR = size(X2,2);
        dof_GSR = nt - rX2;

        if rX2 < size(X2,2)
            fprintf(' Design matrix (with GSR) not full rank (rank=%d/%d, %s)\n', rX2, size(X2,2), sid);
        end

        if cfg.saveRegressorsBP
            csv_regs_bp_gsr = fullfile(outdir, ['RegressorsBP_GSR_' spm_file(EPI,'basename') '.csv']);
            writematrix(X2, csv_regs_bp_gsr);
            fprintf(' RegressorsBP CSV (design used, with GSR) saved (%s): %s\n', sid, csv_regs_bp_gsr);
        end

        if cfg.saveDesignMatMAT
            matX2 = fullfile(outdir, ['DesignMatrix_GSR_' spm_file(EPI,'basename') '.mat']);
            save(matX2, 'X2', 'badVols', 'nt', 'TR', 'sid');
        end

        Yclean2 = regress_out_rankrobust(Y_bp, X2, cfg.writeFloat32, cfg.regressionSolver);
        out_nii2 = fullfile(outdir, [cfg.writePrefixG spm_file(EPI,'filename')]);
        if exist(out_nii2,'file'), delete(out_nii2); end
        write_4d_nii_from_2d(Yclean2, V_EPI_all, out_nii2, cfg.writeFloat32);
        fprintf(' Denoised EPI saved (with GSR, sensitivity) (%s): %s\n', sid, out_nii2);
    end

    %% ---------- QC figure ----------
    fig = figure('Color','w');
    tvec = (1:nt)';

    subplot(3,1,1);
    plot(tvec,FD,'k'); hold on; yline(cfg.FD_thr,'r--');
    title('Framewise Displacement (Power)'); grid on;
    xlabel('Volume'); ylabel('FD (mm)');
    
    if cfg.markBadVolsInQC && ~isempty(badVols)
    yl = ylim;
    for bv = badVols(:)'
        line([bv bv], yl, ...
             'Color', [0.7 0.7 0.7], ...
             'LineStyle', '-', ...
             'LineWidth', 0.5, ...
             'HandleVisibility','off');
    end
    ylim(yl);
end


    subplot(3,1,2);
    plot(tvec, DV_z,'k'); hold on; yline(cfg.DVARSz_thr,'r--');
    title('DVARS (robust z)'); grid on;
    xlabel('Volume'); ylabel('z(DVARS)');
   
    if cfg.markBadVolsInQC && ~isempty(badVols)
    yl = ylim;
    for bv = badVols(:)'
        line([bv bv], yl, ...
             'Color', [0.7 0.7 0.7], ...
             'LineStyle', '-', ...
             'LineWidth', 0.5, ...
             'HandleVisibility','off');
    end
    ylim(yl);
end

  

    subplot(3,1,3);
    plot(tvec, mean(Y2D_raw(:,brain_idx),2,'omitnan'),'k'); hold on;
    plot(tvec, mean(Yclean1(:,brain_idx),2,'omitnan'),'r');
    legend({'Raw (analysis mask)','Clean (no GSR)'},'Location','best');
    title('Mean signal within analysis mask'); grid on;
    xlabel('Volume');
    
     if cfg.markBadVolsInQC && ~isempty(badVols)
    yl = ylim;
    for bv = badVols(:)'
        line([bv bv], yl, ...
             'Color', [0.7 0.7 0.7], ...
             'LineStyle', '-', ...
             'LineWidth', 0.5, ...
             'HandleVisibility','off');
    end
    ylim(yl);
end
   

    qc_png = fullfile(outdir, ['QC_' spm_file(EPI,'basename') '.png']);
    saveas(fig, qc_png);
    close(fig);

    %% ---------- QC table ----------
    qc_csv = fullfile(outdir, ['QC_' spm_file(EPI,'basename') '_stats.csv']);
    Subject = string(spm_file(EPI,'basename'));
    qc_fail = (meanFD > cfg.meanFD_fail) || (scrubPerc > cfg.scrub_failPct);

    solverUsed   = string(cfg.regressionSolver);
    didLCC       = logical(lcc_applied);
    didIntercept = logical(cfg.addIntercept);
    badVolsPath  = string(tern(~isempty(badVolsFile), badVolsFile, ''));

    didT1Primary = logical(t1_primary_used);
    didStrictCSF = ~isempty(StrictCSF);
    strictCSFvox = nStrictCSFvox;

    Ttab = table(Subject, string(sid), nt, TR, ...
                 nMaskVox, didLCC, t1_intersection_done, didT1Primary, ...
                 didStrictCSF, strictCSFvox, ...
                 didIntercept, solverUsed, badVolsPath, ...
                 meanFD, medianFD, maxFD, nScrub, scrubPerc, ...
                 meanDVARSz, medianDVARSz, maxDVARSz, p95DVARSz, ...
                 mean_tSNR, LC_L_tSNR, LC_R_tSNR, ...
                 kW, kC, varWM, varCSF, ...
                 nReg_noGSR, rX1, dof_noGSR, ...
                 rX2, nReg_GSR, dof_GSR, ...
                 qc_fail, ...
                 'VariableNames', {'Subject','SubID','nVols','TR', ...
                                   'nMaskVox','LCC_applied','MaskIntersect_T1brain','T1brain_primary', ...
                                   'StrictCSFmean_in_model','StrictCSF_nVox', ...
                                   'Intercept_in_GLM','Solver','BadVolsFile', ...
                                   'meanFD','medianFD','maxFD','nScrub','Scrubbed_percent', ...
                                   'meanDVARSz','medianDVARSz','maxDVARSz','p95DVARSz', ...
                                   'tSNR_global','tSNR_LC_L','tSNR_LC_R', ...
                                   'WM_PCs','CSF_PCs','VarWM_percent','VarCSF_percent', ...
                                   'nReg_noGSR','rankX_noGSR','DOF_noGSR', ...
                                   'rankX_GSR','nReg_GSR','DOF_GSR', ...
                                   'QC_fail'});
    writetable(Ttab, qc_csv);
    fprintf(' QC CSV saved (%s): %s\n', sid, qc_csv);

    %% ---------- Methods log ----------
    methods_txt = fullfile(outdir, ['METHODS_' spm_file(EPI,'basename') '.txt']);
    write_methods_log(methods_txt, cfg, sid, EPI, MotionFile, C1_in, C2_in, C3_in, TR, ...
                      kW, kC, varWM, varCSF, meanFD, medianFD, maxFD, nScrub, scrubPerc, ...
                      nReg_noGSR, rX1, dof_noGSR, t1_intersection_done, didT1Primary, ...
                      nMaskVox, lcc_applied, matlabVersion, spmVersion, didStrictCSF, strictCSFvox, meanDVARSz, medianDVARSz, maxDVARSz, p95DVARSz);

    %% ---------- cfg snapshot ----------
    if cfg.saveCfgSnapshot
        cfgFile = fullfile(outdir, ['CFG_' spm_file(EPI,'basename') '.mat']);
        save(cfgFile, 'cfg', 'sid', 'EPI', 'MotionFile', 'C1_in', 'C2_in', 'C3_in', ...
             'TR', 'nt', 'matlabVersion', 'spmVersion', 'nMaskVox', 'lcc_applied', 't1_intersection_done', ...
             't1_primary_used', 'nStrictCSFvox', ...
             'badVols', 'meanFD', 'scrubPerc');
    end

    fprintf('\nDONE — Subject %s\n', sid);
end

fprintf('\n=== DONE: all matched subjects processed ===\n');
end

%% ============================ HELPERS ============================

function c = normalize_multiselect(x, allowCancel)
if nargin < 2, allowCancel = false; end
if isequal(x,0)
    if allowCancel
        c = {};
        return;
    else
        error('Selection canceled but this input is required.');
    end
end
if ischar(x), c = {x}; else, c = x; end
end

function L = make_list_from_files(fn_cell, fp)
L = struct('sid',{},'fname',{});
if isempty(fn_cell), return; end
for i = 1:numel(fn_cell)
    sid = get_subid_from_name(fn_cell{i});
    if isempty(sid)
        fprintf('    Skipped (no Sub-ID): %s\n', fn_cell{i});
        continue;
    end
    L(end+1).sid   = sid; %#ok<SAGROW>
    L(end).fname   = fullfile(fp, fn_cell{i});
end
end

function sid = get_subid_from_name(fname)
[~,name,~] = fileparts(fname);
tok = regexp(name, '(?i)(?:sub|ausub)[-_]?(\d+)', 'tokens');
if isempty(tok), sid = ''; else, sid = tok{1}{1}; end
end

function fname = get_fname_for_sid(list, sid)
idx = find(strcmp({list.sid}, sid));
if isempty(idx)
    error('No file found for Subject %s', sid);
elseif numel(idx) > 1
    warning('Multiple files for Subject %s found; using the first.', sid);
end
fname = list(idx(1)).fname;
end

function fname = get_fname_for_sid_or_empty(list, sid)
if isempty(list), fname = ''; return; end
idx = find(strcmp({list.sid}, sid), 1, 'first');
if isempty(idx), fname = ''; else, fname = list(idx).fname; end
end

function out = tern(cond,a,b)
if cond, out=a; else, out=b; end
end

function [C2_epi, C3_epi] = reslice_tissue_masks_to_epi(C2_in, C3_in, V_EPI_mid, sid, cfg)
V_WM  = spm_vol(C2_in);
V_CSF = spm_vol(C3_in);

sameWM  = all(V_WM.dim==V_EPI_mid.dim)  && max(abs(V_WM.mat(:)-V_EPI_mid.mat(:)))<1e-6;
sameCSF = all(V_CSF.dim==V_EPI_mid.dim) && max(abs(V_CSF.mat(:)-V_EPI_mid.mat(:)))<1e-6;

flags.interp=0; flags.wrap=[0 0 0]; flags.mask=0; flags.which=1; flags.mean=0;

if ~sameWM
    if isfield(cfg,'resliceDeleteOldR') && cfg.resliceDeleteOldR
        rwm = spm_file(C2_in,'prefix','r');
        if exist(rwm,'file'), try, delete(rwm); catch, end, end
    end
    spm_reslice(char(V_EPI_mid.fname, V_WM.fname), flags);
    C2_epi = spm_file(C2_in,'prefix','r');
else
    C2_epi = C2_in;
end

if ~sameCSF
    if isfield(cfg,'resliceDeleteOldR') && cfg.resliceDeleteOldR
        rcsf = spm_file(C3_in,'prefix','r');
        if exist(rcsf,'file'), try, delete(rcsf); catch, end, end
    end
    spm_reslice(char(V_EPI_mid.fname, V_CSF.fname), flags);
    C3_epi = spm_file(C3_in,'prefix','r');
else
    C3_epi = C3_in;
end

if ~sameWM || ~sameCSF
    fprintf(' WM/CSF masks resliced to EPI grid (NN) (%s).\n', sid);
end
end

function [brain_idx, brain_mask, lcc_applied] = build_brainmask_from_meanEPI(mu3, valid1xN, dim3, cfg, sid)
mu_vec = mu3(:);
valid_vec = valid1xN(:);

cand = isfinite(mu_vec) & (mu_vec > 0) & valid_vec;
if nnz(cand) < cfg.brain_min_vox
    warning('Brainmask fallback: too few candidates (%d). Using valid mask (%s).', nnz(cand), sid);
    brain_mask = valid_vec;
else
    vals = mu_vec(cand);
    thrP = prctile(vals, cfg.brain_thr_prct);
    thrM = 0.10 * median(vals);
    thr  = max(thrP, thrM);
    brain_mask = cand & (mu_vec > thr);
end

lcc_applied = false;
try
    nx = dim3(1); ny = dim3(2); nz = dim3(3);
    brain3 = reshape(brain_mask, [nx,ny,nz]);
    CC = bwconncomp(brain3, 6);
    if CC.NumObjects > 0
        [~,imax] = max(cellfun(@numel, CC.PixelIdxList));
        brain3 = false(nx,ny,nz);
        brain3(CC.PixelIdxList{imax}) = true;
        brain_mask = brain3(:);
        lcc_applied = true;
    end
catch
end

brain_idx = find(brain_mask);
if isempty(brain_idx)
    warning('Brainmask empty after processing; using valid mask (%s).', sid);
    brain_mask = valid_vec;
    brain_idx = find(brain_mask);
end
end

function t1mask_vec = get_T1_brainmask_in_epi_grid(C1_in, C2_in, C3_in, V_EPI_mid, dim3, cfg, sid, outdir)
t1mask_vec = [];
try
    V1 = spm_vol(C1_in); V2 = spm_vol(C2_in); V3 = spm_vol(C3_in);
    GM = spm_read_vols(V1);
    WM = spm_read_vols(V2);
    CSF= spm_read_vols(V3);

    brain_t1 = (GM + WM + CSF) > cfg.t1_brain_thr;

    Vbm = V1;
    bm_t1 = fullfile(outdir, sprintf('T1brainmask_%s_T1space.nii', sid));
    Vbm.fname = bm_t1;
    Vbm.pinfo = [1;0;0];
    spm_write_vol(Vbm, double(brain_t1));

    flags.interp=0; flags.wrap=[0 0 0]; flags.mask=0; flags.which=1; flags.mean=0;

    if isfield(cfg,'resliceDeleteOldR') && cfg.resliceDeleteOldR
        rbm = spm_file(bm_t1,'prefix','r');
        if exist(rbm,'file'), try, delete(rbm); catch, end, end
    end

    spm_reslice(char(V_EPI_mid.fname, bm_t1), flags);
    bm_epi = spm_file(bm_t1,'prefix','r');

    Vm = spm_vol(bm_epi);
    m  = spm_read_vols(Vm) > 0.5;

    nx = dim3(1); ny = dim3(2); nz = dim3(3);
    if any(size(m) ~= [nx,ny,nz])
        warning('T1 brainmask reslice size mismatch (%s). Ignoring.', sid);
        return;
    end
    t1mask_vec = m(:);

    fprintf(' T1 brainmask built+resliced (%s): thr=%.2f\n', sid, cfg.t1_brain_thr);
catch ME
    warning('T1 brainmask failed (%s): %s', sid, ME.message);
    t1mask_vec = [];
end
end

function [WMvol_out, CSFvol_out] = erode_tissues_if_possible(WMvol, CSFvol, cfg, sid)
WMvol_out = WMvol;
CSFvol_out = CSFvol;

try
    se = strel('sphere', 1);

    if cfg.erodeIterWM > 0
        tmp = WMvol_out;
        for i=1:cfg.erodeIterWM
            tmp = imerode(tmp, se);
        end
        WMvol_out = tmp;
    end

    if cfg.erodeIterCSF > 0
        tmp = CSFvol_out;
        for i=1:cfg.erodeIterCSF
            tmp = imerode(tmp, se);
        end
        CSFvol_out = tmp;
    end

    fprintf(' Tissue erosion applied (%s): WM iter=%d, CSF iter=%d\n', sid, cfg.erodeIterWM, cfg.erodeIterCSF);
catch
    warning('Image Processing Toolbox not available -> skipping tissue erosion (%s).', sid);
end
end

function [comp, k, varExpl] = acompcor_from_mask(Y2D_raw, mask1xNvox, nt, sid, label, cfg)
idx = find(mask1xNvox);
if isempty(idx)
    warning('%s mask has 0 valid voxels — using dummy (%s).', label, sid);
    comp = zeros(nt,1); k=1; varExpl = 100;
    return;
end

Y = zscore(double(Y2D_raw(:, idx)),0,1);
if size(Y,2) == 1
    comp = Y; k=1; varExpl = 100;
    return;
end

% Use SVD-based PCA (stable)
[coeff,score,latent] = pca(Y,'Algorithm','svd','Economy',true);

% STABILITY PATCH: fix sign ambiguity deterministically
if isfield(cfg,'pcaSignFix') && cfg.pcaSignFix
    try
        [~,imax] = max(abs(coeff),[],1);
        sgn = sign(coeff(sub2ind(size(coeff), imax, 1:size(coeff,2))));
        sgn(sgn==0) = 1;
        score = score .* sgn;
    catch
        % fallback: mean-based
        for j=1:size(score,2)
            if mean(score(:,j)) < 0
                score(:,j) = -score(:,j);
            end
        end
    end
end

expPct = 100*latent/sum(latent);
cumPct = cumsum(expPct);
k = find(cumPct>=90, 1, 'first');
k = min(max(k,3),5);
k = min(k, nt-1);
comp = score(:,1:k);
varExpl = sum(expPct(1:k));
end

function M36 = friston36(R)
R_deriv      = [zeros(1,6); diff(R)];
R_sq         = R.^2;
R_deriv_sq   = R_deriv.^2;
R_lag1       = [zeros(1,6); R(1:end-1,:)];
R_deriv_lag1 = [zeros(1,6); R_deriv(1:end-1,:)];
M36 = [R R_deriv R_sq R_deriv_sq R_lag1 R_deriv_lag1];
end

function [FD, meanFD, medianFD, maxFD] = compute_FD_power(R)
R_head = 50; % mm radius
dTrans = [zeros(1,3); abs(diff(R(:,1:3)))];
dRot   = [zeros(1,3); abs(diff(R(:,4:6)))];
FD = sum(dTrans,2) + R_head * sum(dRot,2);
meanFD = mean(FD);
medianFD = median(FD);
maxFD = max(FD);
end

function DV_z = compute_DVARSz(Yvalid)
dY = [zeros(1,size(Yvalid,2)); diff(Yvalid)];
DV = sqrt(mean(dY.^2, 2));  % RMS
medDV = median(DV);
madDV = median(abs(DV - medDV));
robust_sigma = 1.4826 * madDV;

if robust_sigma <= 0 || ~isfinite(robust_sigma)
    robust_sigma = std(DV);
    if robust_sigma == 0 || ~isfinite(robust_sigma)
        robust_sigma = 1;
    end
end


DV_z = (DV - medDV) / robust_sigma;
DV_z(1) = 0;
end


function Yclean = regress_out_rankrobust(Y, X, writeFloat32, solver)
if nargin < 4 || isempty(solver), solver = 'qrp'; end

switch lower(solver)
    case 'qrp'
        % STABILITY PATCH: QR with column pivoting (rank-aware, stable)
        [Q,R,E] = qr(X,0);
        tol = max(size(X))*eps(norm(R,inf));
        r = sum(abs(diag(R)) > tol);

        if r == 0
            Yres = Y;
        else
            b = zeros(size(X,2), size(Y,2));
            b(E(1:r),:) = R(1:r,1:r) \ (Q(:,1:r)' * Y);
            Yres = Y - X*b;
        end

    case 'lsqminnorm'
        if exist('lsqminnorm','file')==2
            beta = lsqminnorm(X, Y);
            Yres = Y - X*beta;
        else
            warning('lsqminnorm not available; falling back to pinv.');
            beta = pinv(X) * Y;
            Yres = Y - X*beta;
        end

    case 'pinv'
        beta = pinv(X) * Y;
        Yres = Y - X*beta;

    otherwise
        warning('Unknown solver "%s". Falling back to qrp.', solver);
        [Q,R,E] = qr(X,0);
        tol = max(size(X))*eps(norm(R,inf));
        r = sum(abs(diag(R)) > tol);
        if r == 0
            Yres = Y;
        else
            b = zeros(size(X,2), size(Y,2));
            b(E(1:r),:) = R(1:r,1:r) \ (Q(:,1:r)' * Y);
            Yres = Y - X*b;
        end
end

if writeFloat32
    Yclean = single(Yres);
else
    Yclean = Yres;
end
end

function write_4d_nii_from_2d(Y2D, V_EPI_all, outname, writeFloat32)
% STABILITY PATCH: streaming writer (no full 4D RAM allocation)
if exist(outname,'file')
    try, delete(outname); catch, end
end

nt = size(Y2D,1);
Vt = V_EPI_all(1);

if writeFloat32
    Vt.dt = [16 0]; % float32
else
    Vt.dt = [64 0]; % float64
end
Vt.pinfo = [1;0;0];

nx = V_EPI_all(1).dim(1);
ny = V_EPI_all(1).dim(2);
nz = V_EPI_all(1).dim(3);

for t=1:nt
    Vt.n = [t 1];
    Vt.fname = outname;
    vol = reshape(Y2D(t,:).', [nx,ny,nz]);
    spm_write_vol(Vt, vol);
end
end

function roi_tsnr = compute_roi_tsnr(maskfile, V_ref, tSNR_map, sid, label)
V_roi = spm_vol(maskfile);
if max(abs(V_roi.mat(:)-V_ref.mat(:)))>1e-6 || any(V_roi.dim~=V_ref.dim)
    flags.interp=0; flags.wrap=[0 0 0]; flags.mask=0; flags.which=1; flags.mean=0;
    spm_reslice(char(V_ref.fname, maskfile), flags);
    maskfile = spm_file(maskfile,'prefix','r');
    V_roi  = spm_vol(maskfile);
end
M = spm_read_vols(V_roi) > 0.5;
if ~any(M(:))
    warning('ROI %s has no voxels >0.5 (Sub %s). tSNR=NaN.', label, sid);
    roi_tsnr = NaN;
else
    roi_tsnr = mean(tSNR_map(M), 'omitnan');
    fprintf('%s tSNR (%s): %.2f\n', label, sid, roi_tsnr);
end
end

function write_methods_log(fname, cfg, sid, EPI, MotionFile, C1_in, C2_in, C3_in, TR, ...
                           kW, kC, varWM, varCSF, meanFD, medianFD, maxFD, nScrub, scrubPerc, ...
                           nReg_noGSR, rankX_noGSR, dof_noGSR, didT1Intersect, didT1Primary, ...
                           nMaskVox, didLCC, matlabVersion, spmVersion, didStrictCSF, strictCSFvox, meanDVARSz, medianDVARSz, maxDVARSz, p95DVARSz)
fid = fopen(fname,'w');
if fid==-1
    warning('Could not write methods log: %s', fname);
    return;
end

fprintf(fid,'LC-RSFC DENOISING — METHODS LOG (v9.4 + STABILITY PATCHES)\n');
fprintf(fid,'Subject: %s\n', sid);

fprintf(fid,'\nSoftware:\n');
fprintf(fid,'  MATLAB: %s\n', matlabVersion);
fprintf(fid,'  SPM: %s\n', spmVersion);

fprintf(fid,'\nInputs:\n');
fprintf(fid,'  EPI: %s\n', EPI);
fprintf(fid,'  Motion: %s\n', MotionFile);
fprintf(fid,'  GM(c1): %s\n', tern(~isempty(C1_in), C1_in, '(not provided)'));
fprintf(fid,'  WM(c2): %s\n', C2_in);
fprintf(fid,'  CSF(c3): %s\n', C3_in);

fprintf(fid,'\nSpace:\n');
fprintf(fid,'  Native EPI grid\n');
fprintf(fid,'  TR: %.6f s (fixed)\n', TR);

fprintf(fid,'\nAnalysis mask:\n');
if didT1Primary
    fprintf(fid,'  PRIMARY: T1 brainmask from (c1+c2+c3)>%.2f, resliced NN to EPI, AND valid voxels\n', cfg.t1_brain_thr);
    fprintf(fid,'  Largest connected component applied (if available): %s\n', tern(didLCC,'ON','OFF'));
else
    fprintf(fid,'  MeanEPI-based within valid voxels; threshold = max(P%d, 0.1*median)\n', cfg.brain_thr_prct);
    fprintf(fid,'  Largest connected component applied (if available): %s\n', tern(didLCC,'ON','OFF'));
    fprintf(fid,'  T1 brainmask intersect applied: %s\n', tern(didT1Intersect,'ON','OFF'));
end
fprintf(fid,'  Mask voxel count (final): %d\n', nMaskVox);

fprintf(fid,'\nTemporal filtering:\n');
fprintf(fid,'  Bandpass: [%.4f %.4f] Hz, Butterworth order %d, filtfilt (zero-phase)\n', cfg.band(1), cfg.band(2), cfg.filter_order);
fprintf(fid,'  Data detrending: linear detrend before bandpass\n');
fprintf(fid,'  Continuous regressors: demean before bandpass, then z-score\n');
fprintf(fid,'  Spike regressors (scrubbing): unfiltered one-hot regressors\n');
fprintf(fid,'  Intercept in GLM (unfiltered): %s\n', tern(cfg.addIntercept,'ON','OFF'));

fprintf(fid,'\nConfound regression:\n');
fprintf(fid,'  Model: Motion (Friston-36) + aCompCor (WM/CSF) + strictCSFmean + scrubbing spikes (+ intercept)\n');
fprintf(fid,'  Solver (rank-aware): %s\n', cfg.regressionSolver);
fprintf(fid,'  Strict CSF mean included: %s (nVox=%d, thr=%.2f)\n', ...
    tern(didStrictCSF,'ON','OFF'), strictCSFvox, cfg.CSF_strict_thr);
fprintf(fid,'  GSR sensitivity pipeline computed: %s (global mean in filtered space, analysis mask)\n', tern(cfg.doGSR_sens,'ON','OFF'));

fprintf(fid,'\naCompCor:\n');
fprintf(fid,'  WM thr=%.2f -> %d PCs, variance explained %.1f%% (clamped 3–5)\n', cfg.WM_thr, kW, varWM);
fprintf(fid,'  CSF thr=%.2f -> %d PCs, variance explained %.1f%% (clamped 3–5)\n', cfg.CSF_thr, kC, varCSF);

fprintf(fid,'\nTissue masks (after reslice to EPI grid):\n');
fprintf(fid,'  Threshold: WM>%.2f, CSF>%.2f\n', cfg.WM_thr, cfg.CSF_thr);
fprintf(fid,'  Optional erosion: %s (WM iter=%d, CSF iter=%d)\n', ...
    tern(cfg.doErodeTissues,'ON','OFF'), cfg.erodeIterWM, cfg.erodeIterCSF);

fprintf(fid,'\nMotion model:\n');
fprintf(fid,'  Friston-36: [R, dR, R^2, dR^2, lag(R), lag(dR)]\n');

fprintf(fid,'\nScrubbing:\n');
fprintf(fid,'  Criteria: FD>%.2f OR DVARS(robust z)>%.2f\n', cfg.FD_thr, cfg.DVARSz_thr);
fprintf(fid,'  nScrub=%d (%.1f%%)\n', nScrub, scrubPerc);
fprintf(fid,'\nQC summary:\n');
fprintf(fid,'  meanFD=%.4f, medianFD=%.4f, maxFD=%.4f\n', meanFD, medianFD, maxFD);
fprintf(fid,'  DVARS (robust z): mean=%.2f, median=%.2f, max=%.2f, p95=%.2f\n', ...
        meanDVARSz, medianDVARSz, maxDVARSz, p95DVARSz);

fprintf(fid,'  Design matrix (no GSR): nReg=%d, rank=%d, DOF=%d (nt - rank)\n', nReg_noGSR, rankX_noGSR, dof_noGSR);

fprintf(fid,'\nStability patches enabled:\n');
fprintf(fid,'  Streaming NIfTI writer: ON\n');
fprintf(fid,'  PCA sign fix: %s\n', tern(cfg.pcaSignFix,'ON','OFF'));
fprintf(fid,'  Reslice delete old r*.nii: %s\n', tern(cfg.resliceDeleteOldR,'ON','OFF'));
fprintf(fid,'  Sanitize finite before filtfilt: %s\n', tern(cfg.sanitizeFinite,'ON','OFF'));
fprintf(fid,'  Strict reproducibility (rng + single-thread): %s\n', tern(cfg.strictRepro,'ON','OFF'));

fclose(fid);
end
