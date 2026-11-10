function warp_drn_and_exclude_csf_subjectmatch_spm12_win()
% SPM12 (Windows GUI): Warp one MNI DRN seed mask into each subject's EPI space
% using inverse deformation fields, then exclude CSF using subject-specific c3 maps.
%
% Inputs (GUI):
%   - One MNI DRN mask (.nii/.nii.gz)
%   - Multiple reference EPIs (e.g., meanEPI*.nii)  [MultiSelect]
%   - Multiple deformation fields y_*.nii (T1->MNI) [MultiSelect]
%   - Multiple CSF maps c3*.nii (T1 space TPM)      [MultiSelect]
%   - Output directory
%
% Subject matching:
%   - Matches by token "sub-<ID>" in filenames (letters/digits allowed)
%
% Outputs per subject:
%   sub-<ID>_DRN_inEPI.nii
%   sub-<ID>_DRN_inEPI_CSFexcluded.nii
% Plus: a CSV log with voxel counts.

% -------------------- init --------------------
spm('Defaults','fMRI');
spm_jobman('initcfg');

% -------------------- choose files --------------------
[maskFile, maskPath] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
    'Select MNI DRN seed mask (single file)');
if isequal(maskFile,0), error('No MNI mask selected.'); end
mask_mni = fullfile(maskPath, maskFile);

[epiFiles, epiPath] = uigetfile({'*.nii','NIfTI (*.nii)'}, ...
    'Select reference EPI images (e.g., meanEPI) - MultiSelect', 'MultiSelect','on');
if isequal(epiFiles,0), error('No EPI references selected.'); end
if ischar(epiFiles), epiFiles = {epiFiles}; end
epi_full = cellfun(@(f) fullfile(epiPath,f), epiFiles, 'UniformOutput', false);

[yFiles, yPath] = uigetfile({'y_*.nii','SPM deformation field (y_*.nii)'}, ...
    'Select deformation fields y_*.nii (T1->MNI) - MultiSelect', 'MultiSelect','on');
if isequal(yFiles,0), error('No deformation fields selected.'); end
if ischar(yFiles), yFiles = {yFiles}; end
y_full = cellfun(@(f) fullfile(yPath,f), yFiles, 'UniformOutput', false);

[c3Files, c3Path] = uigetfile({'c3*.nii;*.nii','CSF TPM (c3*.nii)'}, ...
    'Select CSF maps c3*.nii (T1 space) - MultiSelect', 'MultiSelect','on');
if isequal(c3Files,0), error('No c3 files selected.'); end
if ischar(c3Files), c3Files = {c3Files}; end
c3_full = cellfun(@(f) fullfile(c3Path,f), c3Files, 'UniformOutput', false);

outDir = uigetdir(epiPath,'Select output directory');
if isequal(outDir,0), error('No output directory selected.'); end
if ~exist(outDir,'dir'), mkdir(outDir); end

% -------------------- parameters --------------------
thrCSF = 0.50;     % CSF probability threshold
minVox = 5;        % warn if seed has fewer voxels after exclusion
prefixTmp = 'TMP_'; % temp prefix used to find outputs robustly

% -------------------- subject matching --------------------
% ID pattern: sub-<alphanum>
get_id = @(p) regexp(p,'sub-([A-Za-z0-9]+)','tokens','once');

epi_map = containers.Map();
for i = 1:numel(epi_full)
    id = get_id(epi_full{i});
    if isempty(id), warning('No sub-<ID> in EPI: %s', epi_full{i}); continue; end
    epi_map(id{1}) = epi_full{i};
end

y_map = containers.Map();
for i = 1:numel(y_full)
    id = get_id(y_full{i});
    if isempty(id), warning('No sub-<ID> in y_: %s', y_full{i}); continue; end
    y_map(id{1}) = y_full{i};
end

c3_map = containers.Map();
for i = 1:numel(c3_full)
    id = get_id(c3_full{i});
    if isempty(id), warning('No sub-<ID> in c3: %s', c3_full{i}); continue; end
    c3_map(id{1}) = c3_full{i};
end

ids = intersect(intersect(keys(epi_map), keys(y_map)), keys(c3_map));
if isempty(ids)
    error('No matching subjects found. Ensure filenames include "sub-<ID>" for EPI, y_, and c3.');
end

fprintf('\nMatched subjects (%d):\n', numel(ids));
disp(ids(:));

% -------------------- process each subject --------------------
logRows = {'ID','EPIref','yField','c3File','seed_warped','seed_clean','vox_warped','vox_clean','note'};

for k = 1:numel(ids)
    id = ids{k};
    epi = epi_map(id);
    yField = y_map(id);
    c3 = c3_map(id);

    fprintf('\n=== sub-%s ===\n', id);
    fprintf('EPI: %s\n', epi);
    fprintf('y_:  %s\n', yField);
    fprintf('c3:  %s\n', c3);

    note = "";

    % --- sanity checks ---
    Vy = spm_vol(yField);
    if numel(Vy) < 3
        warning('yField <3 vols for sub-%s. Check deformation field: %s', id, yField);
        note = note + "yField<3vols;";
    end

    % =========================
    % (1) Warp MNI seed -> EPI
    % =========================
    t0 = now;
    clear matlabbatch

    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {yField};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space      = {epi};

    % explicit target grid
    matlabbatch{1}.spm.util.defs.comp{2}.id.space       = {epi};

    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames     = {mask_mni};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {outDir};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp     = 0;    % NN for mask
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask       = 0;    % safer for brainstem
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm       = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix     = sprintf('%sseed_sub-%s_', prefixTmp, id);

    spm_jobman('run', matlabbatch);

    % find warped seed output
    cand = dir(fullfile(outDir, sprintf('%sseed_sub-%s_*', prefixTmp, id)));
    cand = cand([cand.datenum] > t0 - 1/1440);
    if isempty(cand)
        warning('No warped seed output found for sub-%s', id);
        logRows(end+1,:) = {id, epi, yField, c3, '', '', NaN, NaN, 'no_warp_output'}; %#ok<AGROW>
        continue;
    end
    [~,ix] = max([cand.datenum]);
    seed_warp_tmp = fullfile(cand(ix).folder, cand(ix).name);

    seed_warp = fullfile(outDir, sprintf('sub-%s_DRN_inEPI.nii', id));
    if exist(seed_warp,'file'), delete(seed_warp); end
    movefile(seed_warp_tmp, seed_warp);

    % binarize warped seed (defensive)
    Vsw = spm_vol(seed_warp);
    Sw  = spm_read_vols(Vsw);
    Swb = double(Sw > 0.5);
    spm_write_vol(Vsw, Swb);  % overwrite in place as binary
    vox_warped = nnz(Swb);

    % ======================================
    % (2) Resample c3(T1) -> EPI voxel grid
    % ======================================
    t1 = now;
    clear matlabbatch
    matlabbatch{1}.spm.util.defs.comp{1}.id.space = {epi};   % target grid EPI
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {c3};  % resample c3
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {outDir};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;     % trilinear for probability map
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask   = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm   = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = sprintf('%sc3_sub-%s_', prefixTmp, id);

    spm_jobman('run', matlabbatch);

    cand = dir(fullfile(outDir, sprintf('%sc3_sub-%s_*', prefixTmp, id)));
    cand = cand([cand.datenum] > t1 - 1/1440);
    if isempty(cand)
        warning('No c3(EPI) output found for sub-%s', id);
        logRows(end+1,:) = {id, epi, yField, c3, seed_warp, '', vox_warped, NaN, 'no_c3_output'}; %#ok<AGROW>
        continue;
    end
    [~,ix] = max([cand.datenum]);
    c3_epi_tmp = fullfile(cand(ix).folder, cand(ix).name);

    c3_epi = fullfile(outDir, sprintf('sub-%s_c3_inEPI.nii', id));
    if exist(c3_epi,'file'), delete(c3_epi); end
    movefile(c3_epi_tmp, c3_epi);

    % =========================
    % (3) Exclude CSF from seed
    % =========================
    Vc = spm_vol(c3_epi);
    C  = spm_read_vols(Vc);

    csfMask = (C > thrCSF);
    seedClean = Swb & ~csfMask;

    vox_clean = nnz(seedClean);
    if vox_clean < minVox
        warning('sub-%s: Seed after CSF exclusion is very small (%d vox).', id, vox_clean);
        note = note + sprintf("seed_clean_small(%d);", vox_clean);
    end
    if vox_clean == 0
        warning('sub-%s: Seed after CSF exclusion is EMPTY.', id);
        note = note + "seed_clean_EMPTY;";
    end

    seed_clean_file = fullfile(outDir, sprintf('sub-%s_DRN_inEPI_CSFexcluded.nii', id));
    Vout = Vsw;
    Vout.fname = seed_clean_file;
    spm_write_vol(Vout, double(seedClean));

    % cleanup: optionally delete intermediate c3_epi file? (I keep it for QC)
    % optional QC:
    % spm_check_registration(char(epi, seed_clean_file));

    fprintf(' Wrote seed:  %s (vox=%d)\n', seed_warp, vox_warped);
    fprintf(' Wrote clean: %s (vox=%d)\n', seed_clean_file, vox_clean);

    logRows(end+1,:) = {id, epi, yField, c3, seed_warp, seed_clean_file, vox_warped, vox_clean, char(note)}; %#ok<AGROW>
end

% -------------------- save CSV log --------------------
csvFile = fullfile(outDir, 'DRN_seed_warp_and_CSFexclude_log.csv');
try
    writecell(logRows, csvFile);
catch
    fid = fopen(csvFile,'w');
    fprintf(fid,'ID,EPIref,yField,c3File,seed_warped,seed_clean,vox_warped,vox_clean,note\n');
    for r = 2:size(logRows,1)
        fprintf(fid,'%s,"%s","%s","%s","%s","%s",%g,%g,"%s"\n', ...
            logRows{r,1}, logRows{r,2}, logRows{r,3}, logRows{r,4}, ...
            logRows{r,5}, logRows{r,6}, logRows{r,7}, logRows{r,8}, logRows{r,9});
    end
    fclose(fid);
end
fprintf('\n Log saved: %s\n', csvFile);
fprintf('\nDone.\n');
end
