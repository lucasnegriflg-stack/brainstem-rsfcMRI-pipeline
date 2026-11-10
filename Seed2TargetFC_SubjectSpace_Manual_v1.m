function SimpleSeedTargetFC_MetricsFirst_v4()
% SimpleSeedTargetFC_MetricsFirst_v4
%
% Completely standalone seed-to-target FC extraction for denoised resting-state fMRI.
%
% GUI ORDER EXACTLY:
%   1) Select ALL EPI NIfTI files at once
%   2) Select ALL seed NIfTI files at once
%   3) Select ALL analysis-mask NIfTI files at once
%   4) Select target files from as many folders as needed:
%        - choose a target folder
%        - multi-select ALL target NIfTIs in that folder
%        - repeat with more folders
%        - cancel the folder dialog when done
%
% MATCHING:
%   All files are matched only by the subject ID found in the path/name:
%   sub-933, sub-1372, sub-01, sub-ABC123, etc.
%
% IMPORTANT:
%   This version avoids MATLAB table() constructors internally, because some
%   installations/path conflicts can trigger "Incorrect number of arguments"
%   inside table.table. CSV files are written manually.
%
% REQUIREMENT:
%   SPM must be on the MATLAB path.

clear; clc;

spm('defaults','fmri');
spm_jobman('initcfg');

SET = default_settings();

fprintf('\n=== SimpleSeedTargetFC_MetricsFirst_v4 ===\n');
fprintf('No HC/PD. No manifest. Matching only by sub-ID. No table() row construction.\n');

%% Output directory
outBase = uigetdir(pwd, 'Select output directory');
if isequal(outBase,0)
    error('Cancelled: no output directory selected.');
end
outdir = fullfile(outBase, 'results_SimpleSeedTargetFC_MetricsFirst_v4');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fprintf('\nOutput directory:\n  %s\n', outdir);

Issues = cell(0,6);

%% 1) EPIs
fprintf('\nSTEP 1/4: Select ALL EPI NIfTI files at once.\n');
epiRaw = select_multi_nifti(pwd, 'STEP 1/4: Select ALL EPI NIfTI files at once', true);
[epiFiles, epiOriginal, Issues] = prepare_nifti_files(epiRaw, outdir, 'EPI', Issues);
[EPI, Issues] = build_file_items(epiFiles, epiOriginal, 'EPI', '', Issues);
[EPI, Issues] = keep_first_per_subject(EPI, 'EPI', Issues);
fprintf('Usable EPI subjects: %d\n', numel(EPI));
if isempty(EPI)
    error('No usable EPIs after sub-ID matching.');
end

%% 2) Seeds
fprintf('\nSTEP 2/4: Select ALL seed NIfTI files at once.\n');
seedRaw = select_multi_nifti(pwd, 'STEP 2/4: Select ALL seed NIfTI files at once', true);
[seedFiles, seedOriginal, Issues] = prepare_nifti_files(seedRaw, outdir, 'seed', Issues);
[Seeds, Issues] = build_file_items(seedFiles, seedOriginal, 'seed', '', Issues);
[Seeds, Issues] = make_unique_labels_within_subject(Seeds, 'seed', Issues);
fprintf('Usable seed files: %d\n', numel(Seeds));
if isempty(Seeds)
    error('No usable seed files after sub-ID matching.');
end

%% 3) Analysis masks
fprintf('\nSTEP 3/4: Select ALL analysis-mask NIfTI files at once.\n');
fprintf('This is optional. Press Cancel if you do not want analysis masks.\n');
amRaw = select_multi_nifti(pwd, 'STEP 3/4: Select ALL analysis-mask NIfTI files at once, or Cancel to skip', false);
if isempty(amRaw)
    AM = empty_items();
    Issues = add_issue(Issues, 'INFO', 'AnalysisMaskSelection', '', 'analysis_mask', '', 'No analysis masks selected. FC will be computed without analysis masks.');
else
    [amFiles, amOriginal, Issues] = prepare_nifti_files(amRaw, outdir, 'analysis_mask', Issues);
    [AM, Issues] = build_file_items(amFiles, amOriginal, 'analysis_mask', '', Issues);
    [AM, Issues] = keep_first_per_subject(AM, 'analysis_mask', Issues);
end
fprintf('Usable analysis masks: %d\n', numel(AM));

%% 4) Targets from variable number of folders
fprintf('\nSTEP 4/4: Target selection.\n');
fprintf('Choose one target folder, then multi-select ALL target NIfTIs in it. Repeat for more folders.\n');
fprintf('Cancel the folder dialog when all target folders are done.\n');
Targets = empty_items();

while true
    targetFolder = uigetdir(pwd, 'STEP 4/4: Select a target folder, or Cancel when done');
    if isequal(targetFolder,0)
        break
    end

    folderLabel = folder_basename(targetFolder);
    targetRaw = select_multi_nifti(targetFolder, sprintf('Select ALL target NIfTIs in folder: %s', folderLabel), false);

    if isempty(targetRaw)
        Issues = add_issue(Issues, 'WARNING', 'TargetSelection', '', 'target', targetFolder, 'Target folder selected but no target NIfTI selected inside it.');
        continue
    end

    [targetFiles, targetOriginal, Issues] = prepare_nifti_files(targetRaw, outdir, 'target', Issues);
    [TargetBlock, Issues] = build_file_items(targetFiles, targetOriginal, 'target', folderLabel, Issues);
    Targets = [Targets; TargetBlock]; %#ok<AGROW>
end

[Targets, Issues] = make_unique_labels_within_subject(Targets, 'target', Issues);
fprintf('Usable target files: %d\n', numel(Targets));
if isempty(Targets)
    error('No usable target files selected.');
end

%% Match to EPI subjects and save selection reports
epiSubjects = unique_cell({EPI.SubjectID});
[SeedsMatched, Issues] = keep_only_subjects(Seeds, epiSubjects, 'seed', Issues);
[AMMatched, Issues] = keep_only_subjects(AM, epiSubjects, 'analysis_mask', Issues);
[TargetsMatched, Issues] = keep_only_subjects(Targets, epiSubjects, 'target', Issues);

if isempty(SeedsMatched)
    error('No seed files match any selected EPI subject ID.');
end
if isempty(TargetsMatched)
    error('No target files match any selected EPI subject ID.');
end

write_items_csv(fullfile(outdir, 'Matched_EPIs.csv'), EPI, true);
write_items_csv(fullfile(outdir, 'Matched_Seeds.csv'), SeedsMatched, true);
write_items_csv(fullfile(outdir, 'Matched_AnalysisMasks.csv'), AMMatched, true);
write_items_csv(fullfile(outdir, 'Matched_Targets.csv'), TargetsMatched, true);
write_cell_csv(fullfile(outdir, 'QC_matching_issues.csv'), {'Severity','Stage','SubjectID','Role','Path','Message'}, Issues);

%% Compute connectivity
ConnRows = cell(0,17);
SubjectQCRows = cell(0,9);
SeedQCRows = cell(0,9);
TargetQCRows = cell(0,9);

subjects = {EPI.SubjectID};
fprintf('\n=== Computing connectivity ===\n');
fprintf('Subjects with EPI: %d\n', numel(subjects));

for s = 1:numel(subjects)
    sid = subjects{s};
    fprintf('\n[%s] %d/%d\n', sid, s, numel(subjects));

    epiItem = EPI(s);
    seedItems = select_items_by_subject(SeedsMatched, sid);
    targetItems = select_items_by_subject(TargetsMatched, sid);
    amItems = select_items_by_subject(AMMatched, sid);

    if isempty(seedItems)
        note = 'Skipped: no matching seed files for this subject.';
        fprintf('  %s\n', note);
        SubjectQCRows(end+1,:) = {sid, epiItem.Path, '', 0, numel(targetItems), 0, NaN, 'SKIPPED', note}; %#ok<AGROW>
        continue
    end

    if isempty(targetItems)
        note = 'Skipped: no matching target files for this subject.';
        fprintf('  %s\n', note);
        SubjectQCRows(end+1,:) = {sid, epiItem.Path, '', numel(seedItems), 0, 0, NaN, 'SKIPPED', note}; %#ok<AGROW>
        continue
    end

    try
        [ConnSub, SeedQCSub, TargetQCSub, SubjectQCSub, Issues] = compute_subject_fc(sid, epiItem, seedItems, targetItems, amItems, SET, Issues);
        ConnRows = [ConnRows; ConnSub]; %#ok<AGROW>
        SeedQCRows = [SeedQCRows; SeedQCSub]; %#ok<AGROW>
        TargetQCRows = [TargetQCRows; TargetQCSub]; %#ok<AGROW>
        SubjectQCRows = [SubjectQCRows; SubjectQCSub]; %#ok<AGROW>
        fprintf('  FC rows: %d\n', size(ConnSub,1));
        if numel(SubjectQCSub) >= 9
            fprintf('  QC: %s\n', SubjectQCSub{9});
        end
    catch ME
        note = exception_message_full(ME);
        warning('[%s] %s', sid, note);
        amPath = '';
        if ~isempty(amItems)
            amPath = amItems(1).Path;
        end
        SubjectQCRows(end+1,:) = {sid, epiItem.Path, amPath, numel(seedItems), numel(targetItems), 0, NaN, 'FAILED', note}; %#ok<AGROW>
        Issues = add_issue(Issues, 'ERROR', 'ComputeSubject', sid, 'subject', epiItem.Path, note);
    end
end

%% Save outputs
connHeader = {'SubjectID','EPIFile','AnalysisMaskFile','Seed','SeedFile','Target','TargetFile','r','z','nFrames','SeedVoxOrig','SeedVoxInMask','SeedCoveragePct','TargetVoxOrig','TargetVoxInMask','TargetCoveragePct','Notes'};
subjHeader = {'SubjectID','EPIFile','AnalysisMaskFile','nSeedsMatched','nTargetsMatched','nFC','nFrames','Status','Notes'};
seedQCHeader = {'SubjectID','Seed','Path','nVoxOrig','nVoxInMask','CoveragePct','MaskValid','SignalValid','Notes'};
targetQCHeader = {'SubjectID','Target','Path','nVoxOrig','nVoxInMask','CoveragePct','MaskValid','SignalValid','Notes'};

write_cell_csv(fullfile(outdir, 'Connectivity_long.csv'), connHeader, ConnRows);
write_cell_csv(fullfile(outdir, 'QC_subjects.csv'), subjHeader, SubjectQCRows);
write_cell_csv(fullfile(outdir, 'QC_seed_masks.csv'), seedQCHeader, SeedQCRows);
write_cell_csv(fullfile(outdir, 'QC_target_masks.csv'), targetQCHeader, TargetQCRows);
write_cell_csv(fullfile(outdir, 'QC_matching_issues.csv'), {'Severity','Stage','SubjectID','Role','Path','Message'}, Issues);
write_wide_fisherz_csv(fullfile(outdir, 'Connectivity_wide_FisherZ.csv'), ConnRows);

save(fullfile(outdir, 'SimpleSeedTargetFC_MetricsFirst_v4_workspace.mat'), 'SET', 'EPI', 'SeedsMatched', 'AMMatched', 'TargetsMatched', 'ConnRows', 'SubjectQCRows', 'SeedQCRows', 'TargetQCRows', 'Issues');

fprintf('\nDONE. Results written to:\n  %s\n', outdir);
fprintf('Main CSV:\n  %s\n', fullfile(outdir, 'Connectivity_long.csv'));
fprintf('Wide Fisher-z CSV:\n  %s\n', fullfile(outdir, 'Connectivity_wide_FisherZ.csv'));
end

%% =====================================================================
%% Settings
%% =====================================================================
function SET = default_settings()
SET = struct();
SET.seed_thr = 0.5;
SET.target_thr = 0.5;
SET.analysismask_thr = 0.5;
% Metrics-first defaults: compute a metric whenever at least one voxel is available.
% QC CSVs still report voxel counts, coverage and fallbacks.
SET.min_vox_seed = 1;
SET.min_vox_target = 1;
SET.seed_cov_thr_pct = NaN;
SET.target_cov_thr_pct = NaN;
SET.fallback_without_analysismask_if_intersection_fails = true;
SET.epsClamp = 1e-6;
end

%% =====================================================================
%% Compute subject
%% =====================================================================
function [ConnRows, SeedQCRows, TargetQCRows, SubjectQCRow, Issues] = compute_subject_fc(sid, epiItem, seedItems, targetItems, amItems, SET, Issues)
ConnRows = cell(0,17);
SeedQCRows = cell(0,9);
TargetQCRows = cell(0,9);

EPIfile = epiItem.Path;
AMfile = '';
if ~isempty(amItems)
    AMfile = amItems(1).Path;
end

Vepi_all = spm_vol(EPIfile);
nt = numel(Vepi_all);
if nt < 2
    error('EPI has fewer than 2 volumes: %s', EPIfile);
end
Vref = Vepi_all(max(1, round(nt/2)));

Y = spm_read_vols(Vepi_all);
Y2D = reshape(Y, [], nt)';
clear Y;

AM = [];
amNote = 'No analysis mask.';
if ~isempty(AMfile)
    try
        AM = load_binary_mask_to_ref(AMfile, Vref, SET.analysismask_thr);
        amNote = 'Used analysis mask.';
    catch ME
        AM = [];
        amNote = ['Analysis mask failed; computed without analysis mask: ' ME.message];
        Issues = add_issue(Issues, 'WARNING', 'LoadAnalysisMask', sid, 'analysis_mask', AMfile, amNote);
    end
end

nSeeds = numel(seedItems);
nTargets = numel(targetItems);
SeedTS = nan(nt, nSeeds);
TargetTS = nan(nt, nTargets);
SeedValid = false(nSeeds, 1);
TargetValid = false(nTargets, 1);
SeedQC = repmat(empty_qc_struct(), nSeeds, 1);
TargetQC = repmat(empty_qc_struct(), nTargets, 1);

for i = 1:nSeeds
    f = seedItems(i).Path;
    [ts, q] = extract_mask_timeseries(f, Vref, Y2D, AM, SET.seed_thr, SET.min_vox_seed, SET.seed_cov_thr_pct, SET);
    zts = zscore_safe(ts);
    sigValid = all(isfinite(zts));
    if sigValid
        SeedTS(:, i) = zts;
        SeedValid(i) = true;
    end
    SeedQC(i) = q;
    SeedQCRows(end+1,:) = {sid, seedItems(i).Label, f, q.nVoxOrig, q.nVoxInMask, q.CoveragePct, q.MaskValid, sigValid, q.Notes}; %#ok<AGROW>
end

for i = 1:nTargets
    f = targetItems(i).Path;
    [ts, q] = extract_mask_timeseries(f, Vref, Y2D, AM, SET.target_thr, SET.min_vox_target, SET.target_cov_thr_pct, SET);
    zts = zscore_safe(ts);
    sigValid = all(isfinite(zts));
    if sigValid
        TargetTS(:, i) = zts;
        TargetValid(i) = true;
    end
    TargetQC(i) = q;
    TargetQCRows(end+1,:) = {sid, targetItems(i).Label, f, q.nVoxOrig, q.nVoxInMask, q.CoveragePct, q.MaskValid, sigValid, q.Notes}; %#ok<AGROW>
end

for si = 1:nSeeds
    if ~SeedValid(si)
        continue
    end
    for ti = 1:nTargets
        if ~TargetValid(ti)
            continue
        end
        r = pearson_r_complete(SeedTS(:,si), TargetTS(:,ti));
        if isempty(r) || ~isfinite(r)
            r = NaN;
            z = NaN;
        else
            r = min(max(r, -1 + SET.epsClamp), 1 - SET.epsClamp);
            z = atanh(r);
        end
        ConnRows(end+1,:) = { ...
            sid, EPIfile, AMfile, ...
            seedItems(si).Label, seedItems(si).Path, ...
            targetItems(ti).Label, targetItems(ti).Path, ...
            r, z, nt, ...
            SeedQC(si).nVoxOrig, SeedQC(si).nVoxInMask, SeedQC(si).CoveragePct, ...
            TargetQC(ti).nVoxOrig, TargetQC(ti).nVoxInMask, TargetQC(ti).CoveragePct, ...
            amNote}; %#ok<AGROW>
    end
end

if isempty(ConnRows)
    status = 'NO_VALID_FC';
    note = sprintf('No valid seed-target FC rows. Valid seeds=%d/%d, valid targets=%d/%d. %s', sum(SeedValid), nSeeds, sum(TargetValid), nTargets, amNote);
else
    status = 'OK';
    note = sprintf('Computed %d FC rows. Valid seeds=%d/%d, valid targets=%d/%d. %s', size(ConnRows,1), sum(SeedValid), nSeeds, sum(TargetValid), nTargets, amNote);
end
SubjectQCRow = {sid, EPIfile, AMfile, nSeeds, nTargets, size(ConnRows,1), nt, status, note};
end

function q = empty_qc_struct()
q = struct('nVoxOrig', NaN, 'nVoxInMask', NaN, 'CoveragePct', NaN, 'MaskValid', false, 'Notes', '');
end

%% =====================================================================
%% ROI extraction and resampling
%% =====================================================================
function [ts, q] = extract_mask_timeseries(maskFile, Vref, Y2D, AM, thr, minVox, covThrPct, SET)
q = empty_qc_struct();
ts = nan(size(Y2D,1),1);

try
    Mroi = load_binary_mask_to_ref(maskFile, Vref, thr);
catch ME
    q.Notes = ['Invalid: could not load/resample mask: ' ME.message];
    return
end

q.nVoxOrig = nnz(Mroi);
if q.nVoxOrig == 0
    q.nVoxInMask = 0;
    q.CoveragePct = NaN;
    q.Notes = 'Invalid: mask has zero voxels before analysis-mask intersection.';
    return
end

usedFallback = false;
if isempty(AM)
    M = Mroi;
    q.nVoxInMask = nnz(M);
    q.CoveragePct = 100;
    maskNote = 'No analysis mask used.';
else
    Mint = Mroi & AM;
    nInt = nnz(Mint);
    covInt = 100 * nInt / q.nVoxOrig;
    useIntersection = true;

    if nInt < minVox
        useIntersection = false;
    end
    if isfinite(covThrPct) && covInt < covThrPct
        useIntersection = false;
    end

    if useIntersection
        M = Mint;
        q.nVoxInMask = nInt;
        q.CoveragePct = covInt;
        maskNote = 'Used ROI intersected with analysis mask.';
    elseif isfield(SET, 'fallback_without_analysismask_if_intersection_fails') && SET.fallback_without_analysismask_if_intersection_fails
        M = Mroi;
        q.nVoxInMask = nnz(M);
        q.CoveragePct = 100;
        usedFallback = true;
        maskNote = sprintf('Fallback: ROI intersected with analysis mask had %d voxels, so ROI was computed WITHOUT analysis mask.', nInt);
    else
        M = Mint;
        q.nVoxInMask = nInt;
        q.CoveragePct = covInt;
        q.Notes = sprintf('Invalid: only %d voxels after analysis-mask intersection; minimum required is %d.', q.nVoxInMask, minVox);
        return
    end
end

if q.nVoxInMask < minVox
    q.Notes = sprintf('Invalid: only %d voxels available; minimum required is %d.', q.nVoxInMask, minVox);
    return
end

Yroi = Y2D(:, M(:));
ts = mean_omitnan_rows(Yroi);
if any(~isfinite(ts))
    q.Notes = ['Invalid: mean ROI time series contains non-finite values. ' maskNote];
    return
end
if safe_std(ts) == 0 || ~isfinite(safe_std(ts))
    q.Notes = ['Invalid: constant or non-finite mean ROI time series. ' maskNote];
    return
end
q.MaskValid = true;
if usedFallback
    q.Notes = ['OK_WITH_FALLBACK. ' maskNote];
else
    q.Notes = ['OK. ' maskNote];
end
end

function M = load_binary_mask_to_ref(maskFile, Vref, thr)
Vsrc_all = spm_vol(maskFile);
Vsrc = Vsrc_all(1);

same_dim = isequal(double(Vref.dim(:)), double(Vsrc.dim(:)));
same_mat = max(abs(Vref.mat(:) - Vsrc.mat(:))) < 1e-6;

if same_dim && same_mat
    M = spm_read_vols(Vsrc) > thr;
    M = logical(M);
    return
end

dim = double(Vref.dim(:)');
nVox = prod(dim);
vals = zeros(nVox,1);
chunk = 200000;
invSrcMat = inv(Vsrc.mat);

for startIdx = 1:chunk:nVox
    stopIdx = min(startIdx + chunk - 1, nVox);
    idx = startIdx:stopIdx;
    [x, y, z] = ind2sub(dim, idx(:));
    xyzRefVox = [double(x)'; double(y)'; double(z)'; ones(1,numel(idx))];
    xyzWorld = Vref.mat * xyzRefVox;
    xyzSrcVox = invSrcMat * xyzWorld;
    sampled = spm_sample_vol(Vsrc, xyzSrcVox(1,:), xyzSrcVox(2,:), xyzSrcVox(3,:), 0);
    vals(idx) = sampled(:);
end

M = reshape(vals, dim) > thr;
M = logical(M);
end

%% =====================================================================
%% Math helpers
%% =====================================================================
function z = zscore_safe(x)
x = double(x(:));
if any(~isfinite(x))
    z = nan(size(x));
    return
end
sd = safe_std(x);
if sd == 0 || ~isfinite(sd)
    z = nan(size(x));
else
    z = (x - mean(x)) ./ sd;
end
end

function r = pearson_r_complete(x, y)
x = double(x(:));
y = double(y(:));
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
if numel(x) < 2
    r = NaN;
    return
end
x = x - mean(x);
y = y - mean(y);
denom = sqrt(sum(x.^2) * sum(y.^2));
if denom == 0 || ~isfinite(denom)
    r = NaN;
else
    r = sum(x .* y) / denom;
end
end

function s = safe_std(x)
x = double(x(:));
x = x(isfinite(x));
if numel(x) < 2
    s = NaN;
else
    s = std(x);
end
end

function m = mean_omitnan_rows(Y)
finiteMask = isfinite(Y);
counts = sum(finiteMask, 2);
Y(~finiteMask) = 0;
sums = sum(Y, 2);
m = sums ./ counts;
m(counts == 0) = NaN;
end

%% =====================================================================
%% File selection and matching
%% =====================================================================
function files = select_multi_nifti(startDir, prompt, required)
if nargin < 3
    required = false;
end
if isempty(startDir) || ~isfolder(startDir)
    startDir = pwd;
end
oldDir = pwd;
cleanupObj = onCleanup(@() cd(oldDir)); %#ok<NASGU>
cd(startDir);

[fs, path] = uigetfile( ...
    {'*.nii;*.nii.gz;*.img', 'NIfTI files (*.nii, *.nii.gz, *.img)'; '*.*', 'All files (*.*)'}, ...
    prompt, ...
    'MultiSelect', 'on');

if isequal(fs,0)
    if required
        error('Required file selection cancelled: %s', prompt);
    else
        files = {};
        return
    end
end
if ischar(fs)
    fs = {fs};
end
files = cell(numel(fs),1);
for i = 1:numel(fs)
    files{i} = fullfile(path, fs{i});
end
end

function [filesOut, originalOut, Issues] = prepare_nifti_files(filesIn, outdir, role, Issues)
filesOut = {};
originalOut = {};
unzipDir = fullfile(outdir, 'unzipped_nii_gz');
if ~exist(unzipDir, 'dir')
    mkdir(unzipDir);
end

for i = 1:numel(filesIn)
    f = filesIn{i};
    if isempty(f) || ~isfile(f)
        Issues = add_issue(Issues, 'ERROR', 'PrepareFiles', '', role, f, 'Selected file does not exist. Skipped.');
        continue
    end
    fproc = f;
    if ends_with_ignorecase(f, '.nii.gz')
        try
            gunzip(f, unzipDir);
            [~, name1, ext1] = fileparts(f);
            if strcmpi(ext1, '.gz')
                [~, name2, ext2] = fileparts(name1);
                fproc = fullfile(unzipDir, [name2 ext2]);
            else
                fproc = fullfile(unzipDir, name1);
            end
            Issues = add_issue(Issues, 'INFO', 'Gunzip', '', role, f, ['Unzipped to: ' fproc]);
        catch ME
            Issues = add_issue(Issues, 'ERROR', 'Gunzip', '', role, f, ['Could not unzip .nii.gz file: ' ME.message '. Skipped.']);
            continue
        end
    end
    filesOut{end+1,1} = fproc; %#ok<AGROW>
    originalOut{end+1,1} = f; %#ok<AGROW>
end
end

function [Items, Issues] = build_file_items(files, original, role, folderLabel, Issues)
Items = empty_items();
for i = 1:numel(files)
    sid = extract_subid(files{i});
    if isempty(sid)
        Issues = add_issue(Issues, 'ERROR', 'BuildItems', '', role, files{i}, 'No sub-ID found. Skipped.');
        continue
    end
    item = empty_item();
    item.SubjectID = sid;
    item.Role = role;
    item.Path = files{i};
    item.OriginalPath = original{i};
    item.FolderLabel = folderLabel;

    if strcmpi(role, 'seed')
        item.Label = derive_seed_label(files{i}, sid);
    elseif strcmpi(role, 'target')
        item.Label = derive_target_label(files{i}, sid, folderLabel);
    else
        item.Label = role;
    end
    Items(end+1,1) = item; %#ok<AGROW>
end
end

function item = empty_item()
item = struct('SubjectID', '', 'Role', '', 'Label', '', 'FolderLabel', '', 'Path', '', 'OriginalPath', '');
end

function Items = empty_items()
Items = repmat(empty_item(), 0, 1);
end

function [ItemsOut, Issues] = keep_first_per_subject(Items, role, Issues)
ItemsOut = empty_items();
subjects = unique_cell({Items.SubjectID});
for s = 1:numel(subjects)
    sid = subjects{s};
    idx = find(strcmp({Items.SubjectID}, sid));
    if isempty(idx)
        continue
    end
    ItemsOut(end+1,1) = Items(idx(1)); %#ok<AGROW>
    if numel(idx) > 1
        for k = 2:numel(idx)
            Issues = add_issue(Issues, 'WARNING', 'DuplicatePerSubject', sid, role, Items(idx(k)).Path, sprintf('Multiple %s files for this subject. Kept first, skipped this one.', role));
        end
    end
end
end

function [Items, Issues] = make_unique_labels_within_subject(Items, role, Issues)
subjects = unique_cell({Items.SubjectID});
for s = 1:numel(subjects)
    sid = subjects{s};
    idxSub = find(strcmp({Items.SubjectID}, sid));
    labels = {Items(idxSub).Label};
    uniq = unique_cell(labels);
    for u = 1:numel(uniq)
        idxSameLocal = find(strcmp(labels, uniq{u}));
        if numel(idxSameLocal) <= 1
            continue
        end
        for k = 2:numel(idxSameLocal)
            rowIndex = idxSub(idxSameLocal(k));
            oldLabel = Items(rowIndex).Label;
            newLabel = sprintf('%s__dup%d', oldLabel, k);
            Items(rowIndex).Label = newLabel;
            Issues = add_issue(Issues, 'WARNING', 'DuplicateROILabel', sid, role, Items(rowIndex).Path, sprintf('Duplicate %s label within subject. Renamed %s to %s.', role, oldLabel, newLabel));
        end
    end
end
end

function [ItemsOut, Issues] = keep_only_subjects(Items, subjects, role, Issues)
ItemsOut = empty_items();
for i = 1:numel(Items)
    if any(strcmp(subjects, Items(i).SubjectID))
        ItemsOut(end+1,1) = Items(i); %#ok<AGROW>
    else
        Issues = add_issue(Issues, 'WARNING', 'NoMatchingEPI', Items(i).SubjectID, role, Items(i).Path, 'File has a subject ID that does not match any selected EPI. Skipped.');
    end
end
end

function ItemsSub = select_items_by_subject(Items, sid)
if isempty(Items)
    ItemsSub = empty_items();
    return
end
idx = strcmp({Items.SubjectID}, sid);
ItemsSub = Items(idx);
end

function sid = extract_subid(pathstr)
m = regexp(pathstr, '(sub-[A-Za-z0-9]+)', 'match', 'once');
if isempty(m)
    sid = '';
else
    sid = m;
end
end

function label = derive_seed_label(pathstr, sid)
base = filename_without_nii_ext(pathstr);
label = strip_roi_label(base, sid);
if isempty(label)
    label = base;
end
label = make_safe_label(label);
end

function label = derive_target_label(pathstr, sid, folderLabel)
base = filename_without_nii_ext(pathstr);
cleanBase = strip_roi_label(base, sid);
folderLabel = make_safe_label(folderLabel);

if isempty(cleanBase) || is_generic_roi_name(cleanBase)
    label = folderLabel;
elseif contains_ignorecase(cleanBase, folderLabel) || contains_ignorecase(folderLabel, cleanBase)
    label = folderLabel;
else
    label = [folderLabel '__' make_safe_label(cleanBase)];
end
label = make_safe_label(label);
end

function tf = is_generic_roi_name(x)
x = lower(make_safe_label(x));
generic = {'roi','mask','target','native','space_native','space_T1w','space_epi'};
tf = any(strcmp(x, generic));
end

function label = strip_roi_label(base, sid)
label = base;
label = regexprep(label, ['(?i)' regexptranslate('escape', sid)], '');
label = regexprep(label, '(?i)^sub-[A-Za-z0-9]+[_\-]*', '');
label = regexprep(label, '(?i)^(seed|roi|mask|native|space-native|space_T1w|space-MNI|space_epi|target)[_\-]*', '');
label = regexprep(label, '(?i)[_\-]*(seed|roi|mask|native|target)$', '');
label = regexprep(label, '^[_\-\s]+|[_\-\s]+$', '');
end

function base = filename_without_nii_ext(pathstr)
[~, name, ext] = fileparts(pathstr);
if strcmpi(ext, '.gz')
    [~, name2, ext2] = fileparts(name);
    if strcmpi(ext2, '.nii')
        base = name2;
    else
        base = name;
    end
else
    base = name;
end
end

function label = make_safe_label(label)
label = char(label);
label = regexprep(label, '[^A-Za-z0-9]+', '_');
label = regexprep(label, '^_+|_+$', '');
if isempty(label)
    label = 'ROI';
end
end

function name = folder_basename(pathstr)
[~, name] = fileparts(pathstr);
if isempty(name)
    name = 'TargetFolder';
end
name = make_safe_label(name);
end

function c = unique_cell(c)
if isempty(c)
    c = {};
    return
end
[~, ia] = unique(c, 'stable');
c = c(sort(ia));
% Restore true stable order manually because sort(ia) is increasing original positions.
end

function tf = contains_ignorecase(a,b)
tf = ~isempty(strfind(lower(a), lower(b))); %#ok<STREMP>
end

function tf = ends_with_ignorecase(s, suffix)
s = lower(s);
suffix = lower(suffix);
tf = numel(s) >= numel(suffix) && strcmp(s(end-numel(suffix)+1:end), suffix);
end

%% =====================================================================
%% CSV writers
%% =====================================================================
function write_items_csv(outfile, Items, includeOriginal)
if nargin < 3
    includeOriginal = true;
end
if includeOriginal
    header = {'SubjectID','Role','Label','FolderLabel','Path','OriginalPath'};
else
    header = {'SubjectID','Role','Label','FolderLabel','Path'};
end
rows = cell(numel(Items), numel(header));
for i = 1:numel(Items)
    if includeOriginal
        rows(i,:) = {Items(i).SubjectID, Items(i).Role, Items(i).Label, Items(i).FolderLabel, Items(i).Path, Items(i).OriginalPath};
    else
        rows(i,:) = {Items(i).SubjectID, Items(i).Role, Items(i).Label, Items(i).FolderLabel, Items(i).Path};
    end
end
write_cell_csv(outfile, header, rows);
end

function write_cell_csv(outfile, header, rows)
fid = fopen(outfile, 'w');
if fid < 0
    error('Could not open CSV for writing: %s', outfile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

write_csv_line(fid, header);
for i = 1:size(rows,1)
    write_csv_line(fid, rows(i,:));
end
end

function write_csv_line(fid, row)
for j = 1:numel(row)
    if j > 1
        fprintf(fid, ',');
    end
    fprintf(fid, '%s', csv_escape(row{j}));
end
fprintf(fid, '\n');
end

function s = csv_escape(x)
if isnumeric(x)
    if isempty(x) || (isscalar(x) && isnan(x))
        s = '';
    elseif isscalar(x)
        s = sprintf('%.12g', x);
    else
        s = mat2str(x);
    end
elseif islogical(x)
    if x
        s = 'true';
    else
        s = 'false';
    end
elseif isstring(x)
    if ismissing(x) || strlength(x) == 0
        s = '';
    else
        s = char(x);
    end
elseif ischar(x)
    s = x;
else
    try
        s = char(string(x));
    catch
        s = '';
    end
end
s = strrep(s, '"', '""');
if any(s == ',') || any(s == '"') || any(s == sprintf('\n')) || any(s == sprintf('\r'))
    s = ['"' s '"'];
end
end

function write_wide_fisherz_csv(outfile, ConnRows)
% ConnRows columns:
% 1 SubjectID, 4 Seed, 6 Target, 9 z
if isempty(ConnRows)
    write_cell_csv(outfile, {'SubjectID'}, cell(0,1));
    return
end
subjects = unique_cell(ConnRows(:,1)');
features = cell(size(ConnRows,1),1);
for i = 1:size(ConnRows,1)
    features{i} = make_safe_label(['FCz_' ConnRows{i,4} '__' ConnRows{i,6}]);
end
features = unique_cell(features');

header = [{'SubjectID'}, features(:)'];
rows = cell(numel(subjects), numel(header));
for s = 1:numel(subjects)
    rows{s,1} = subjects{s};
    for f = 1:numel(features)
        rows{s,f+1} = NaN;
    end
end

for i = 1:size(ConnRows,1)
    sid = ConnRows{i,1};
    feat = make_safe_label(['FCz_' ConnRows{i,4} '__' ConnRows{i,6}]);
    z = ConnRows{i,9};
    ri = find(strcmp(subjects, sid), 1, 'first');
    ci = find(strcmp(features, feat), 1, 'first') + 1;
    rows{ri,ci} = z;
end
write_cell_csv(outfile, header, rows);
end

%% =====================================================================
%% Issues and errors
%% =====================================================================
function Issues = add_issue(Issues, severity, stage, sid, role, path, message)
Issues(end+1,:) = {severity, stage, sid, role, path, message};
end

function msg = exception_message_full(ME)
msg = ['FAILED: ' ME.message];
if ~isempty(ME.stack)
    parts = cell(numel(ME.stack),1);
    for i = 1:numel(ME.stack)
        parts{i} = sprintf('%s line %d', ME.stack(i).name, ME.stack(i).line);
    end
    msg = [msg ' | Stack: ' strjoin(parts, ' <- ')];
end
end
