function LC_seedLR_to_voxel_native_batch_v2_3_1()
%% ==========================================================
%% LC_seedLR_to_voxel_native_batch_v2_3_1.m
%% PAPER-READY v2.3.1 — RSFC Seed-to-Voxel (native EPI space, L/R, BATCH)
%% ==========================================================
% v2.3.1 patches vs v2.3.0:
%   PATCH 1: valid-mask not overly strict (finite policy handled later)
%   PATCH 2: BandAmp metric no longer removes timepoints inside FFT (no ts(isfinite))
%   PATCH 3: tmp-reslice is now FAIL-FAST if copy to tmp fails (prevents r* collisions on shares)
%   PATCH 4: nValidCorrVox computed from final r_mask (nnz(isfinite(r_mask))) not nMask_good proxy
%
% Core properties:
%   • Native EPI space workflow
%   • Preferred masking via ANALYSISMASK_*.nii from denoising
%   • Seed time series from Seed ∩ AnalysisMask
%   • Correlations computed only within analysis mask
%   • Fisher-z maps for group level
%
% Requirements: SPM12 on path
%% ==========================================================

clear; clc;
spm('defaults','fmri'); spm_jobman('initcfg');
fprintf('\n=== START: LC Seed-to-Voxel (L/R, native space, v2.3.1 PAPER-READY) ===\n\n');

%% ---------- provenance (MATLAB/SPM versions) ----------
matlabVersion = version;
spmVersion = 'unknown';
try, spmVersion = spm('Ver'); catch, end

%% -------------------- CONFIG --------------------
cfg = struct();

cfg.TR = 2.21;                           % fixed TR
cfg.use_analysisMask_input = true;       % recommended: use ANALYSISMASK_*.nii
cfg.min_seed_vox = 1;                    % minimum voxels required in Seed∩Mask
cfg.seed_winsor_prc = [1 99];            % set [] to disable winsorization
cfg.lc_inmask_warnPct = 90;              % warn if <90% seed voxels inside analysis mask (pre-reslice)
cfg.compute_t_p = false;                 % OFF by default (reviewer-safe)
cfg.hist_bins = 80;                      % QC histogram bins

% Band-amplitude settings (seed) — ALFF-like metric, explicitly defined below.
cfg.band_amp_name = 'BandAmp_RS';
cfg.band_amp_band = [0.008 0.09];

% Optional brainmask fallback from T1 segments (if ANALYSISMASK not supplied)
cfg.make_brainmask_from_segments = true;
cfg.seg_thr = 0.2;                       % (c1+c2+c3) > 0.2

% Last-resort fallback mask from mean EPI
cfg.meanEPI_thr_prct = 20;               % threshold at 20th percentile within candidates
cfg.meanEPI_min_vox = 1000;              % if fewer than this -> fall back to valid voxels

% Optional: keep largest connected component (needs Image Processing Toolbox)
cfg.keep_largest_cc = true;

% Optional intersection with T1 brainmask — applied ONLY for fallback masks (not for ANALYSISMASK_input)
cfg.intersect_with_T1_brainmask = true;
cfg.t1_brain_thr = 0.2;                  % (c1+c2+c3) > 0.2

% Save provenance snapshot per subject
cfg.saveCfgSnapshot = true;

% NaN policy (paper-ready explicit behavior)
cfg.drop_vox_with_any_nan = true;        % exclude voxels with any NaN/Inf across time before zscore/corr
cfg.drop_tp_if_seed_nonfinite = true;    % drop timepoints where seed TS is NaN/Inf (global)

% reslice tmp folder
cfg.use_tmp_reslice_dir = true;

%% ---------- 1) Select files (MultiSelect) ----------
fprintf('Select required files (multi-select allowed):\n');

% EPI (denoised, native 4D)
[fn_EPI, fp_EPI] = uigetfile({'*.nii','4D-fMRI (*.nii)'}, ...
    'Select denoised 4D fMRI (e.g., dw_*.nii)', 'MultiSelect','on');
if isequal(fn_EPI,0), error('Aborted: no EPI selected.'); end
fn_EPI = normalize_multiselect(fn_EPI);

% LC seeds L
[fn_LCL, fp_LCL] = uigetfile({'*.nii','LC_L seeds (*.nii)'}, ...
    'Select LC_L seeds (Cancel = skip)', 'MultiSelect','on');
fn_LCL = normalize_multiselect(fn_LCL, true);

% LC seeds R
[fn_LCR, fp_LCR] = uigetfile({'*.nii','LC_R seeds (*.nii)'}, ...
    'Select LC_R seeds (Cancel = skip)', 'MultiSelect','on');
fn_LCR = normalize_multiselect(fn_LCR, true);

if isempty(fn_LCL) && isempty(fn_LCR)
    error('No LC seeds selected (L and/or R required).');
end

% ANALYSISMASK (recommended; from denoising script)
[fn_AM, fp_AM] = uigetfile({'ANALYSISMASK_*.nii;*.nii','Analysis masks (*.nii)'}, ...
    'Select optional ANALYSISMASK_*.nii (Cancel = skip)', 'MultiSelect','on');
fn_AM = normalize_multiselect(fn_AM, true);

% Optional T1 segments c1–c3 (fallback brainmask, native T1 space)
[fn_seg, fp_seg] = uigetfile({'c1*.nii;c2*.nii;c3*.nii','Segments c1–c3 (*.nii)'}, ...
    'Select optional c1–c3 (Cancel = skip)', 'MultiSelect','on');
fn_seg = normalize_multiselect(fn_seg, true);

%% ---------- 2) Build lists (Sub-ID matching) ----------
fprintf('\nBuilding lists from filenames (Sub-ID matching)...\n');

EPI_list   = make_list_from_files(fn_EPI, fp_EPI);
LC_L_list  = make_list_from_files(fn_LCL, fp_LCL);
LC_R_list  = make_list_from_files(fn_LCR, fp_LCR);
AM_list    = make_list_from_files(fn_AM, fp_AM);
SEG_list   = make_list_from_files(fn_seg, fp_seg);

if isempty(EPI_list)
    error('EPI list empty after parsing subject IDs.');
end

all_EPI = {EPI_list.sid};
all_LCL = {LC_L_list.sid};
all_LCR = {LC_R_list.sid};

lc_union = union(all_LCL, all_LCR);
common_sid = intersect(all_EPI, lc_union);

if isempty(common_sid)
    error('No subjects with EPI + (LC_L or LC_R) found.');
end

fprintf('\n Matched subjects (EPI + seed L/R):\n');
disp(common_sid(:));

%% ---------- Master table ----------
allRows = table();

%% ---------- 3) Loop over subjects ----------
for s = 1:numel(common_sid)
    sid = common_sid{s};
    fprintf('\n=============================================\n');
    fprintf('▶ Subject: %s\n', sid);
    fprintf('=============================================\n');

    EPI  = get_fname_for_sid(EPI_list, sid);
    LC_L = get_fname_for_sid_or_empty(LC_L_list, sid);
    LC_R = get_fname_for_sid_or_empty(LC_R_list, sid);

    AM = '';
    if cfg.use_analysisMask_input
        AM = get_fname_for_sid_or_empty(AM_list, sid);
    end

    segpaths = get_segpaths_for_sid(SEG_list, sid); % may be empty

    outdir = fileparts(EPI);
    [~, base] = fileparts(EPI);

    fprintf('EPI:  %s\n', EPI);
    fprintf('LC_L: %s\n', tern(~isempty(LC_L), LC_L, '(none)'));
    fprintf('LC_R: %s\n', tern(~isempty(LC_R), LC_R, '(none)'));
    fprintf('ANALYSISMASK: %s\n', tern(~isempty(AM), AM, '(none)'));

    % tmp dir for safe reslices / intermediates
    tmpdir = fullfile(outdir, sprintf('tmp_reslice_LCSeed_%s', sid));
    if cfg.use_tmp_reslice_dir
        if ~exist(tmpdir,'dir'), mkdir(tmpdir); end
    else
        tmpdir = outdir;
    end

    %% ---------- save cfg/software snapshot per subject ----------
    if isfield(cfg,'saveCfgSnapshot') && cfg.saveCfgSnapshot
        cfgFile = fullfile(outdir, sprintf('CFG_LCSeed_v2_3_1_%s.mat', base));
        try
            save(cfgFile, 'cfg', 'sid', 'EPI', 'LC_L', 'LC_R', 'AM', 'segpaths', 'matlabVersion', 'spmVersion');
        catch ME
            warning('Could not write cfg snapshot (%s): %s', sid, ME.message);
        end
    end

    
    %% ---------- Load EPI header + data ----------
% SAFE read: copy to tmp + volume-wise load
[Y4, Vepi_all] = read_4d_safe_spm(EPI, tmpdir);

nt = numel(Vepi_all);
Vepi_mid = Vepi_all(round(nt/2));
fprintf('TR fixed: %.4f s | nt=%d\n', cfg.TR, nt);

[nx,ny,nz,nt2] = size(Y4);
if nt2 ~= nt, nt = nt2; end

Y2D = reshape(Y4, [], nt)';   % [nt x Nvox]
clear Y4;

    

    % ---------------- PATCH 1 ----------------
    % valid voxels (nonzero and nonflat) — finite policy handled later (good_vox)
    Y2D_d = double(Y2D);
    sdY = std(Y2D_d, 0, 1, 'omitnan');
    nonzero = any(Y2D_d~=0, 1);
    valid = isfinite(sdY) & (sdY > 0) & nonzero;     % [1 x Nvox]

    % mean EPI (for fallback mask)
    mu_vec = mean(Y2D_d, 1, 'omitnan'); % [1 x Nvox]
    mu3 = reshape(mu_vec, [nx,ny,nz]);

    %% ---------- Build analysis mask in EPI grid ----------
    [mask_vec, mask_source] = get_analysis_mask_in_epi_grid(AM, segpaths, Vepi_mid, ...
        [nx,ny,nz], cfg, sid, outdir, tmpdir, mu3, valid);

    mask_vec = mask_vec(:)' & valid; % EPI-analysis mask ∩ valid

    % intersect with T1 brainmask ONLY for fallback masks
    doExtraT1Intersect = cfg.intersect_with_T1_brainmask;
    if contains(string(mask_source), "ANALYSISMASK_input")
        doExtraT1Intersect = false;
    end
    if doExtraT1Intersect
        t1mask_vec = get_T1_brainmask_in_epi_grid(segpaths, Vepi_mid, [nx,ny,nz], cfg, sid, tmpdir);
        if ~isempty(t1mask_vec)
            mask_vec = mask_vec & t1mask_vec(:)' & valid;
            mask_source = mask_source + "∩T1brain";
        else
            warning('No usable T1 brainmask (c1–c3) for %s -> using analysis mask only.', sid);
        end
    end

    nMask = nnz(mask_vec);
    idx_mask = find(mask_vec);

    if nMask < 1000
        warning('Analysis mask very small (%d vox). Check masks/preprocessing (%s).', nMask, sid);
    end
    fprintf(' Using analysis mask (%s): %d voxels\n', string(mask_source), nMask);

    %% ---------- Process sides ----------
    sides = {'L','R'};
    for si = 1:2
        side = sides{si};

        if strcmp(side,'L')
            seedfile = LC_L;
        else
            seedfile = LC_R;
        end

        if isempty(seedfile)
            continue;
        end

        fprintf('\n--- Side %s (%s) ---\n', side, sid);

        %% Seed to EPI grid (NN if needed) — via tmpdir to avoid collisions
        seed_epi = ensure_in_epi_grid_safe(seedfile, Vepi_mid, sid, sprintf('LC_%s', side), tmpdir);

        Vseed = spm_vol(seed_epi);
        Mseed = spm_read_vols(Vseed) > 0.5;
        nSeed_all = nnz(Mseed);

        if nSeed_all < cfg.min_seed_vox
            warning('Seed has %d voxels (<%d). Skipping (%s, %s).', nSeed_all, cfg.min_seed_vox, sid, side);
            continue;
        end

        % Seed overlap with analysis mask
        Mseed_vec = Mseed(:)';
        inside = nnz(Mseed_vec & mask_vec);
        perc_inMask = 100 * inside / max(1, nSeed_all);

        fprintf('Seed voxels total: %d | inside analysis mask: %d (%.1f%%)\n', nSeed_all, inside, perc_inMask);
        if perc_inMask < cfg.lc_inmask_warnPct
            warning('LC_%s seed has <%.0f%% inside analysis mask (%.1f%%) (%s).', side, cfg.lc_inmask_warnPct, perc_inMask, sid);
        end

        % Restrict seed to mask for time series extraction
        Mseed_use = Mseed_vec & mask_vec;
        nSeed_use = nnz(Mseed_use);
        if nSeed_use < cfg.min_seed_vox
            warning('Seed∩Mask has %d voxels. Skipping (%s, %s).', nSeed_use, sid, side);
            continue;
        end

        % QC overlay (in EPI grid)
        overlay = zeros(nx*ny*nz,1,'single');
        overlay(mask_vec)  = 0.5;
        overlay(Mseed_vec) = overlay(Mseed_vec) + 0.8;  % seed (all)
        overlay(Mseed_use) = overlay(Mseed_use) + 0.7;  % seed∩mask (brighter)
        Vq = Vseed;
        Vq.fname = fullfile(outdir, sprintf('QC_LC_%s_in_AnalysisMask_%s.nii', side, sid));
        spm_write_vol(Vq, reshape(overlay, [nx,ny,nz]));
        fprintf(' QC overlay written: %s\n', Vq.fname);

        %% ---------- Seed time series (NaN-robust) ----------
        Yseed = Y2D_d(:, Mseed_use);  % [nt x nSeed_use], double

        % Optional winsorization (NaN-safe, per voxel)
        if ~isempty(cfg.seed_winsor_prc)
            Yseed = winsorize_cols_nan_safe(Yseed, cfg.seed_winsor_prc(1), cfg.seed_winsor_prc(2));
        end

        % Seed average TS
        LC_ts_raw = mean(Yseed, 2, 'omitnan');

        % Drop timepoints with nonfinite seed (global, explicit)
        good_tp = true(nt,1);
        if cfg.drop_tp_if_seed_nonfinite
            good_tp = isfinite(LC_ts_raw);
            if nnz(~good_tp) > 0
                warning('Dropping %d timepoints with nonfinite seed TS (%s, %s).', nnz(~good_tp), sid, side);
            end
        end

        LC_ts_raw = LC_ts_raw(good_tp);
        nt_eff = numel(LC_ts_raw);

        % Seed z-score (NaN-safe)
        LC_ts = zscore_omitnan_1d(LC_ts_raw);

        % Seed voxel distribution (finite values only)
        seed_vals = Yseed(good_tp, :);
        seed_vals = seed_vals(isfinite(seed_vals));

        f = figure('Color','w','Position',[200 200 500 400]);
        histogram(seed_vals, 60);
        xlabel('BOLD signal (seed voxels)');
        ylabel('Count');
        title(sprintf('LC-%s seed voxel distribution (%s)', side, sid));
        grid on;

        png_seed_hist = fullfile(outdir, sprintf('QC_LC_%s_seed_histogram_%s.png', side, base));
        saveas(f, png_seed_hist);
        close(f);

        % voxelwise variance within seed (using good_tp)
        LC_var_vox = mean(var(Yseed(good_tp,:), 0, 1, 'omitnan'), 'omitnan');

        save(fullfile(outdir, sprintf('LC_%s_timeseries_%s.mat', side, base)), ...
            'LC_ts','LC_ts_raw','good_tp','Yseed','Mseed','mask_vec','nSeed_all','nSeed_use','perc_inMask','cfg');

        %% ---------- Band-amplitude metric (ALFF-like, explicit definition) ----------
        [LC_BandAmp_RS, freqs_h, LC_var_ts] = compute_band_amplitude_one_sided(LC_ts_raw, cfg.TR, cfg.band_amp_band);

        save(fullfile(outdir, sprintf('LC_%s_%s_%s.mat', side, cfg.band_amp_name, base)), ...
            'LC_BandAmp_RS','LC_var_ts','LC_var_vox','freqs_h','cfg');

        %% ---------- Voxelwise correlation inside mask (NaN-robust policy) ----------
        Ymask = Y2D_d(:, idx_mask);       % [nt x nMask]
        Ymask = Ymask(good_tp, :);        % [nt_eff x nMask]

        % Exclude voxels with any NaN/Inf across time
        good_vox = true(1, size(Ymask,2));
        if cfg.drop_vox_with_any_nan
            good_vox = all(isfinite(Ymask), 1);
        end

        % Exclude constant voxels (sd==0)
        sd_mask = std(Ymask, 0, 1, 'omitnan');
        good_vox = good_vox & isfinite(sd_mask) & (sd_mask > 0);

        nMask_good = nnz(good_vox);
        if nMask_good < 1000
            warning('Very few good voxels for correlation (%d). Check NaNs/mask/preproc (%s).', nMask_good, sid);
        end

        % z-score voxelwise (only good voxels)
        Ymask_good = Ymask(:, good_vox);
        Ymask_z = zscore_omitnan_cols(Ymask_good);

        if ~all(isfinite(LC_ts)) || numel(LC_ts) ~= nt_eff
            error('Seed TS invalid after preprocessing (%s, %s).', sid, side);
        end

        % correlation r = dot(zx,zy)/(n-1)
        r_good = (Ymask_z' * LC_ts) / (nt_eff - 1);

        % Clamp for Fisher z
        r_good = max(min(r_good, 0.999999), -0.999999);
        z_good = atanh(r_good);

        % Put back into mask vector
        r_mask = nan(1, nMask);
        z_mask = nan(1, nMask);
        r_mask(good_vox) = r_good;
        z_mask(good_vox) = z_good;

        r_vec = nan(1, nx*ny*nz);
        z_vec = nan(1, nx*ny*nz);
        r_vec(idx_mask) = r_mask;
        z_vec(idx_mask) = z_mask;

        t_vec = [];
        p_vec = [];
        df_eff = NaN;
        if cfg.compute_t_p
            rho = corr(LC_ts(1:end-1), LC_ts(2:end), 'Rows','complete');
            if ~isfinite(rho), rho = 0; end
            df_eff = max(10, floor(nt_eff * (1 - rho^2) / (1 + rho^2)));

            t_good = r_good .* sqrt(df_eff ./ (1 - r_good.^2));
            p_good = 2 * tcdf(-abs(t_good), df_eff);

            t_mask = nan(1, nMask); p_mask = nan(1, nMask);
            t_mask(good_vox) = t_good;
            p_mask(good_vox) = p_good;

            t_vec = nan(1, nx*ny*nz);
            p_vec = nan(1, nx*ny*nz);
            t_vec(idx_mask) = t_mask;
            p_vec(idx_mask) = p_mask;
        end

        %% ---------- Write maps ----------
        Vout = Vepi_all(1);
        Vout.dt = [16 0]; % float32 outputs

        write_map(Vout, reshape(single(r_vec), [nx,ny,nz]), fullfile(outdir, sprintf('LC_%s_r_%s.nii', side, base)));
        write_map(Vout, reshape(single(z_vec), [nx,ny,nz]), fullfile(outdir, sprintf('LC_%s_z_%s.nii', side, base)));

        f = figure('Color','w','Position',[200 200 900 300]);
        plot(LC_ts_raw,'k','LineWidth',1.2);
        xlabel('Volume (after dropping nonfinite seed TP)');
        ylabel('Seed BOLD');
        title(sprintf('LC-%s seed BOLD time series (%s)', side, sid));
        grid on;

        png_ts = fullfile(outdir, sprintf('QC_LC_%s_timeseries_%s.png', side, base));
        saveas(f, png_ts);
        close(f);

        if cfg.compute_t_p
            write_map(Vout, reshape(single(t_vec), [nx,ny,nz]), fullfile(outdir, sprintf('LC_%s_t_%s.nii', side, base)));
            write_map(Vout, reshape(single(p_vec), [nx,ny,nz]), fullfile(outdir, sprintf('LC_%s_p_%s.nii', side, base)));
        end

        %% ---------- QC figure ----------
        % ---------------- PATCH 4 ----------------
        nValidCorrVox = nnz(isfinite(r_mask));   % final valid correlation voxels

        f = figure('Color','w');
        subplot(2,1,1);
        plot(LC_ts_raw,'LineWidth',1.2);
        title(sprintf('LC_%s seed time series — n=%d vox (seed∩mask), %s', side, nSeed_use, sid));
        xlabel('Volume (after drop)'); grid on;

        subplot(2,1,2);
        z_map = reshape(z_vec, [nx,ny,nz]);
        zv = z_map(mask_vec);
        zv = zv(isfinite(zv));
        histogram(zv, cfg.hist_bins);
        xlabel('Fisher z'); ylabel('Voxel count'); grid on;
        title(sprintf('Distribution of Fisher z within analysis mask (valid corr vox: %d)', nValidCorrVox));

        qc_png = fullfile(outdir, sprintf('QC_LC_%s_seed_%s.png', side, base));
        saveas(f, qc_png);
        close(f);

        %% ---------- QC CSV per side ----------
        csvfile = fullfile(outdir, sprintf('QC_LC_%s_%s.csv', side, base));
        Subject = string(base);
        Side    = string(side);

        nTP_dropped = nnz(~good_tp);

        T = table(Subject, string(sid), Side, ...
            nSeed_all, nSeed_use, perc_inMask, ...
            nMask, string(mask_source), ...
            nValidCorrVox, nTP_dropped, nt_eff, ...
            LC_BandAmp_RS, LC_var_ts, LC_var_vox, ...
            cfg.TR, df_eff, ...
            'VariableNames', { ...
            'Subject','SubID','Side', ...
            'nSeed_total','nSeed_inMask','Seed_inMask_percent', ...
            'nMaskVox','MaskSource', ...
            'nValidCorrVox','nTimepointsDropped_seedNonfinite','nt_effective', ...
            cfg.band_amp_name,'LC_Var_TimeSeries','LC_Var_Voxels', ...
            'TR','Effective_DF_exploratory'});

        writetable(T, csvfile);
        fprintf(' QC CSV saved: %s\n', csvfile);

        allRows = [allRows; T]; %#ok<AGROW>

        %% ---------- Provenance text log ----------
        prov_txt = fullfile(outdir, sprintf('PROVENANCE_LCSeed_v2_3_1_%s_%s.txt', sid, side));
        try
            fid = fopen(prov_txt,'w');
            fprintf(fid, 'LC Seed-to-Voxel v2.3.1 provenance\n');
            fprintf(fid, 'DateTime: %s\n', datestr(now));
            fprintf(fid, 'MATLAB: %s\n', matlabVersion);
            fprintf(fid, 'SPM: %s\n', spmVersion);
            fprintf(fid, 'SubjectID: %s\n', sid);
            fprintf(fid, 'Side: %s\n', side);
            fprintf(fid, 'EPI: %s\n', EPI);
            fprintf(fid, 'SeedOriginal: %s\n', seedfile);
            fprintf(fid, 'SeedUsed(EPIgrid): %s\n', seed_epi);
            fprintf(fid, 'AnalysisMaskInput: %s\n', tern(~isempty(AM), AM, '(none)'));
            fprintf(fid, 'MaskSource: %s\n', string(mask_source));
            fprintf(fid, 'TR: %.6f\n', cfg.TR);
            fprintf(fid, 'nt_original: %d\n', nt);
            fprintf(fid, 'nt_effective: %d\n', nt_eff);
            fprintf(fid, 'nTimepointsDropped_seedNonfinite: %d\n', nTP_dropped);
            fprintf(fid, 'nMaskVox: %d\n', nMask);
            fprintf(fid, 'nValidCorrVox: %d\n', nValidCorrVox);
            fprintf(fid, 'nSeed_total: %d\n', nSeed_all);
            fprintf(fid, 'nSeed_inMask: %d\n', nSeed_use);
            fprintf(fid, 'Seed_inMask_percent: %.3f\n', perc_inMask);
            fprintf(fid, '%s_band: [%.6f %.6f]\n', cfg.band_amp_name, cfg.band_amp_band(1), cfg.band_amp_band(2));
            fprintf(fid, '%s_value: %.6f\n', cfg.band_amp_name, LC_BandAmp_RS);
            fclose(fid);
            fprintf(' Provenance log written: %s\n', prov_txt);
        catch ME
            warning('Could not write provenance log: %s', ME.message);
        end
    end

    fprintf('\nDONE — Subject %s outputs in: %s\n', sid, outdir);
end

%% ---------- Master CSV ----------
if ~isempty(allRows)
    master_csv = fullfile(fp_EPI, 'QC_LC_ALL_SUBJECTS_LR_v2_3_1.csv');
    writetable(allRows, master_csv);
    fprintf('\n Master QC CSV saved: %s\n', master_csv);
else
    fprintf('\n No entries for master CSV.\n');
end

fprintf('\n=== DONE: LC Seed-to-Voxel v2.3.1 ===\n');
end

%% ================== HELPERS ==================

function c = normalize_multiselect(x, allowCancel)
if nargin < 2, allowCancel = false; end
if isequal(x,0)
    if allowCancel, c = {}; return; end
    error('Selection canceled but this input is required.');
end
if ischar(x), c = {x}; else, c = x; end
end

function L = make_list_from_files(fn_cell, fp)
L = struct('sid',{},'fname',{});
if isempty(fn_cell), return; end
for i = 1:numel(fn_cell)
    sid = get_subid_from_name(fn_cell{i});
    if isempty(sid)
        fprintf('   Skipped (no Sub-ID): %s\n', fn_cell{i});
        continue;
    end
    L(end+1).sid = sid; %#ok<SAGROW>
    L(end).fname = fullfile(fp, fn_cell{i});
end
end

function sid = get_subid_from_name(fname)
% Robust Sub-ID parsing for patterns like:
%   ..._usub-933... , ..._ausub-933... , sub-933_... , sub_933...
% Works even if preceded/followed by '_' (word char), avoids \b pitfalls.

[~,name,~] = fileparts(fname);

% 1) match usub-/ausub-/sub- patterns with safe boundaries:
%    - left boundary: start OR non-alnum
%    - right boundary: next char is non-alnum OR end
pat = '(?i)(?:^|[^a-z0-9])(?:au)?sub[-_]?([a-z0-9]+)(?=[^a-z0-9]|$)';
tok = regexp(name, pat, 'tokens', 'once');
if ~isempty(tok)
    sid = tok{1};
    return;
end

% 2) explicit usub- (covers your cT1_usub-933... case reliably)
pat2 = '(?i)(?:^|[^a-z0-9])usub[-_]?([a-z0-9]+)(?=[^a-z0-9]|$)';
tok = regexp(name, pat2, 'tokens', 'once');
if ~isempty(tok)
    sid = tok{1};
    return;
end

% 3) fallback: first long digit group anywhere (>=2 digits)
tok = regexp(name, '(\d{2,})', 'tokens', 'once');
if ~isempty(tok)
    sid = tok{1};
else
    sid = '';
end
end
function fname = get_fname_for_sid_or_empty(list, sid)
fname = '';
if isempty(list), return; end
idx = find(strcmp({list.sid}, sid), 1, 'first');
if ~isempty(idx)
    fname = list(idx).fname;
end
end

function fname = get_fname_for_sid(list, sid)
fname = get_fname_for_sid_or_empty(list, sid);
if isempty(fname)
    error('No file found for Subject %s', sid);
end
end




function segpaths = get_segpaths_for_sid(SEG_list, sid)
segpaths = {};
if isempty(SEG_list), return; end
idx = find(strcmp({SEG_list.sid}, sid));
if isempty(idx), return; end
segpaths = {SEG_list(idx).fname};

c1 = segpaths(contains(lower(segpaths), [filesep 'c1']) | contains(lower(segpaths), 'c1'));
c2 = segpaths(contains(lower(segpaths), [filesep 'c2']) | contains(lower(segpaths), 'c2'));
c3 = segpaths(contains(lower(segpaths), [filesep 'c3']) | contains(lower(segpaths), 'c3'));
ordered = {};
if ~isempty(c1), ordered{end+1} = c1{1}; end %#ok<AGROW>
if ~isempty(c2), ordered{end+1} = c2{1}; end %#ok<AGROW>
if ~isempty(c3), ordered{end+1} = c3{1}; end %#ok<AGROW>
rest = setdiff(segpaths, ordered, 'stable');
segpaths = [ordered rest];
end

function out = tern(cond,a,b)
if cond, out=a; else, out=b; end
end

function [mask_vec, source] = get_analysis_mask_in_epi_grid(AM, segpaths, Vepi_mid, dim3, cfg, sid, outdir, tmpdir, mu3, valid1xN)
% Returns logical mask vector in EPI grid.

% 1) Preferred: supplied ANALYSISMASK_*.nii
if ~isempty(AM) && isfile(AM)
    maskfile = ensure_in_epi_grid_safe(AM, Vepi_mid, sid, 'ANALYSISMASK', tmpdir);
    Vm = spm_vol(maskfile);
    m = spm_read_vols(Vm) > 0.5;
    mask_vec = m(:);
    source = "ANALYSISMASK_input";
    return;
end

% 2) Fallback: build from c1–c3 segments if available
if cfg.make_brainmask_from_segments && numel(segpaths) >= 3
    c1 = ''; c2 = ''; c3 = '';
    for i=1:numel(segpaths)
        [~,nm,~] = fileparts(segpaths{i});
        if startsWith(lower(nm),'c1') && isempty(c1), c1 = segpaths{i}; end
        if startsWith(lower(nm),'c2') && isempty(c2), c2 = segpaths{i}; end
        if startsWith(lower(nm),'c3') && isempty(c3), c3 = segpaths{i}; end
    end
    if ~isempty(c1) && ~isempty(c2) && ~isempty(c3)
        gm = spm_read_vols(spm_vol(c1));
        wm = spm_read_vols(spm_vol(c2));
        csf= spm_read_vols(spm_vol(c3));
        bm = (gm + wm + csf) > cfg.seg_thr;

        Vbm = spm_vol(c1);
        bm_t1 = fullfile(tmpdir, sprintf('BrainMask_auto_%s_T1space.nii', sid));
        Vbm.fname = bm_t1;
        spm_write_vol(Vbm, double(bm));

        bm_epi = ensure_in_epi_grid_safe(bm_t1, Vepi_mid, sid, 'BrainMaskSegments', tmpdir);
        Vm = spm_vol(bm_epi);
        m = spm_read_vols(Vm) > 0.5;
        mask_vec = m(:);
        source = "Segments_c1c2c3";
        return;
    end
end

% 3) Last resort: mean EPI thresholding within candidates
mu = mu3(:);
validv = valid1xN(:);
cand = isfinite(mu) & (mu > 0) & validv;

if nnz(cand) < cfg.meanEPI_min_vox
    warning('MeanEPI fallback: too few candidates (%d). Using valid voxels (%s).', nnz(cand), sid);
    mask_vec = validv;
    source = "ValidVoxels_only";
else
    thr = prctile(mu(cand), cfg.meanEPI_thr_prct);
    mask_vec = cand & (mu > thr);
    source = sprintf("MeanEPI_prct%d", cfg.meanEPI_thr_prct);
end

% Optional largest connected component
if cfg.keep_largest_cc
    try
        m3 = reshape(mask_vec, dim3);
        CC = bwconncomp(m3, 6);
        if CC.NumObjects > 0
            [~,imax] = max(cellfun(@numel, CC.PixelIdxList));
            m3 = false(dim3);
            m3(CC.PixelIdxList{imax}) = true;
            mask_vec = m3(:);
            source = source + "_LCC";
        end
    catch
    end
end
end

function outmask = ensure_in_epi_grid_safe(infile, Vepi_mid, sid, label, tmpdir)
% Safe reslice to avoid r* collisions by copying to tmpdir first.
V = spm_vol(infile);
same = all(V.dim==Vepi_mid.dim) && max(abs(V.mat(:)-Vepi_mid.mat(:)))<1e-6;
outmask = infile;

if same
    return
end

fprintf(' %s ≠ EPI grid → reslice NN via tmp (%s)\n', label, sid);

[~,nm,ext] = fileparts(infile);
stamp = datestr(now,'yyyymmdd_HHMMSSFFF');
tmp_in = fullfile(tmpdir, sprintf('%s_%s_%s%s', nm, label, stamp, ext));

% ---------------- PATCH 3 ----------------
% FAIL-FAST if copy to tmp fails (prevents r* collisions on shares)
try
    copyfile(infile, tmp_in);
catch ME
    error('Could not copy to tmpdir for safe reslice (%s, %s). %s', sid, label, ME.message);
end

flags.interp=0; flags.wrap=[0 0 0]; flags.mask=0; flags.which=1; flags.mean=0;
spm_reslice(char(Vepi_mid.fname, tmp_in), flags);

outmask = spm_file(tmp_in,'prefix','r');
end

function write_map(Vout, map3, fname)
Vout.fname = fname;
spm_write_vol(Vout, map3);
fprintf(' wrote: %s\n', fname);
end

function [BandAmp, freqs_h, var_ts] = compute_band_amplitude_one_sided(ts, TR, band)
% Explicit band-amplitude metric (ALFF-like).
% Steps:
% 1) linear detrend + demean
% 2) FFT
% 3) one-sided amplitude spectrum A(f) = (2/N)*|FFT(f)| for f in (0, Fs/2)
% 4) BandAmp = sum(A(f) within [band(1), band(2)])

Fs = 1/TR;
ts = double(ts(:));

% ---------------- PATCH 2 ----------------
% Do NOT remove samples; if nonfinite exist, return NaN (policy handled upstream via good_tp)
if any(~isfinite(ts))
    BandAmp = NaN; freqs_h = []; var_ts = NaN;
    return
end

N  = numel(ts);
if N < 8
    BandAmp = NaN; freqs_h = []; var_ts = NaN;
    return
end

ts = detrend(ts, 'linear');
ts = ts - mean(ts);

fft_vals = fft(ts);
freqs = (0:N-1) * (Fs/N);

half = 2:floor(N/2); % exclude DC (0 Hz) by starting at 2
freqs_h = freqs(half);

A = (2/N) * abs(fft_vals(half)); % one-sided amplitude

low_idx = freqs_h >= band(1) & freqs_h <= band(2);
BandAmp = sum(A(low_idx));

var_ts = var(ts);
end

function Z = zscore_omitnan_cols(X)
mu = mean(X, 1, 'omitnan');
sd = std(X, 0, 1, 'omitnan');
sd(sd==0 | ~isfinite(sd)) = NaN;
Z = (X - mu) ./ sd;
end

function z = zscore_omitnan_1d(x)
x = double(x(:));
mu = mean(x, 'omitnan');
sd = std(x, 0, 'omitnan');
if ~isfinite(sd) || sd==0
    z = nan(size(x));
else
    z = (x - mu) ./ sd;
end
end

function Y = winsorize_cols_nan_safe(Y, pLo, pHi)
Y = double(Y);
[nT, nV] = size(Y);
for v = 1:nV
    col = Y(:,v);
    fin = isfinite(col);
    if nnz(fin) < max(5, round(0.2*nT))
        continue
    end
    lo = prctile(col(fin), pLo);
    hi = prctile(col(fin), pHi);
    col(fin) = max(min(col(fin), hi), lo);
    Y(:,v) = col;
end
end

function t1mask_vec = get_T1_brainmask_in_epi_grid(segpaths, Vepi_mid, dim3, cfg, sid, tmpdir)
t1mask_vec = [];

if isempty(segpaths) || numel(segpaths) < 3
    return
end

c1 = ''; c2 = ''; c3 = '';
for i=1:numel(segpaths)
    [~,nm,~] = fileparts(segpaths{i});
    if startsWith(lower(nm),'c1') && isempty(c1), c1 = segpaths{i}; end
    if startsWith(lower(nm),'c2') && isempty(c2), c2 = segpaths{i}; end
    if startsWith(lower(nm),'c3') && isempty(c3), c3 = segpaths{i}; end
end
if isempty(c1) || isempty(c2) || isempty(c3)
    return
end

V1 = spm_vol(c1); V2 = spm_vol(c2); V3 = spm_vol(c3);
GM = spm_read_vols(V1);
WM = spm_read_vols(V2);
CSF= spm_read_vols(V3);

brain_t1 = (GM + WM + CSF) > cfg.t1_brain_thr;

Vbm = V1;
bm_t1 = fullfile(tmpdir, sprintf('T1brainmask_%s_T1space.nii', sid));
Vbm.fname = bm_t1;
spm_write_vol(Vbm, double(brain_t1));

bm_epi = ensure_in_epi_grid_safe(bm_t1, Vepi_mid, sid, 'T1brainmask', tmpdir);

Vm = spm_vol(bm_epi);
m  = spm_read_vols(Vm) > 0.5;

nx = dim3(1); ny = dim3(2); nz = dim3(3);
if any(size(m) ~= [nx,ny,nz])
    warning('T1 brainmask reslice size mismatch (%s). Ignoring.', sid);
    return
end

t1mask_vec = m(:);
end
function [Y4, V4] = read_4d_safe_spm(epi4d, tmpdir)
% Robust reading: copy to tmp (avoids locks) + read volume-wise (reduces memmap stress)

if ~exist(tmpdir,'dir'), mkdir(tmpdir); end

[~,nm,ext] = fileparts(epi4d);
stamp = datestr(now,'yyyymmdd_HHMMSSFFF');
local_copy = fullfile(tmpdir, sprintf('%s_copy_%s%s', nm, stamp, ext));

try
    copyfile(epi4d, local_copy);
    epi_use = local_copy;
catch
    % if copy fails, fall back to original
    epi_use = epi4d;
end

V4 = spm_vol(epi_use);
nt = numel(V4);
dim = V4(1).dim;

% prealloc single to keep memory reasonable
Y4 = zeros([dim nt], 'single');

for t = 1:nt
    % reading one volume at a time is usually much more stable on Windows
    Y4(:,:,:,t) = single(spm_read_vols(V4(t)));
end
end
