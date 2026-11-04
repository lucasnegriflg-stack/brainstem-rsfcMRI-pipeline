function CreateBilateralMasks_MatchBySubID()
% CreateBilateralMasks_MatchBySubID.m
% ============================================================
% PURPOSE
%   Create bilateral NIfTI masks from left and right ROI masks.
%
% WORKFLOW
%   1) Select ALL LEFT mask NIfTI files at once.
%   2) Select ALL RIGHT mask NIfTI files at once.
%   3) Select output folder.
%   4) Match left/right masks by subject ID extracted from file paths.
%   5) For each subject:
%        - if left and right exist: create union mask left OR right
%        - if only one side exists: still create bilateral output from that side
%        - write warnings to QC CSV
%
% REQUIREMENTS
%   - MATLAB
%   - SPM on MATLAB path
%   - NIfTI files as .nii, not .nii.gz
%
% OUTPUTS
%   - sub-XXX_<label>_bilateral.nii
%   - BilateralMask_QC.csv
%   - BilateralMask_Issues.csv
%   - Selected_LeftMasks.csv
%   - Selected_RightMasks.csv
% ============================================================

clear; clc;

if exist('spm','file') ~= 2
    error('SPM is not on the MATLAB path. Add SPM first, then run this script again.');
end

spm('defaults','fmri');
spm_jobman('initcfg');

fprintf('\n=== Create bilateral masks by subject ID ===\n');
fprintf('Subject IDs are extracted as sub-XX / sub-XXX / sub-XXXX / etc.\n');
fprintf('If one side is missing, the available side is still written as bilateral with a warning.\n\n');

% ---------------- SETTINGS ----------------
answer = inputdlg( ...
    {'Output mask label, e.g. Amygdala or DRN:', 'Binarization threshold:'}, ...
    'Bilateral mask settings', [1 80], {'BilateralROI', '0.5'});
if isempty(answer)
    error('Cancelled.');
end

maskLabel = strtrim(answer{1});
if isempty(maskLabel)
    maskLabel = 'BilateralROI';
end
maskLabel = char(matlab.lang.makeValidName(maskLabel));

thr = str2double(answer{2});
if ~isfinite(thr)
    error('Invalid binarization threshold.');
end

% ---------------- SELECT FILES ----------------
leftFiles = select_multi_files(pwd, ...
    'STEP 1/3: Select ALL LEFT mask NIfTI files (.nii) using multi-selection', ...
    true, '*.nii');

rightFiles = select_multi_files(pwd, ...
    'STEP 2/3: Select ALL RIGHT mask NIfTI files (.nii) using multi-selection', ...
    true, '*.nii');

outdir = uigetdir(pwd, 'STEP 3/3: Select output folder for bilateral masks');
if isequal(outdir,0)
    error('Cancelled: no output folder selected.');
end

if ~exist(outdir,'dir')
    mkdir(outdir);
end

fprintf('\nOutput folder:\n  %s\n', outdir);
fprintf('Left masks selected:  %d\n', numel(leftFiles));
fprintf('Right masks selected: %d\n', numel(rightFiles));

% ---------------- GROUP BY SUBJECT ----------------
[leftMap, leftSelected, issues] = group_files_by_subject(leftFiles, 'left', {});
[rightMap, rightSelected, issues] = group_files_by_subject(rightFiles, 'right', issues);

writecell([{'SubjectID','Side','Path'}; leftSelected], fullfile(outdir, 'Selected_LeftMasks.csv'));
writecell([{'SubjectID','Side','Path'}; rightSelected], fullfile(outdir, 'Selected_RightMasks.csv'));

allSubjects = unique([keys(leftMap), keys(rightMap)]);
allSubjects = sort_subject_ids(allSubjects);

fprintf('Matched unique subjects: %d\n', numel(allSubjects));

% ---------------- PROCESS SUBJECTS ----------------
qcRows = {};

for i = 1:numel(allSubjects)
    sid = allSubjects{i};

    if isKey(leftMap, sid)
        Lfiles = leftMap(sid);
    else
        Lfiles = {};
    end

    if isKey(rightMap, sid)
        Rfiles = rightMap(sid);
    else
        Rfiles = {};
    end

    fprintf('\n[%s] %d/%d | left files: %d | right files: %d\n', ...
        sid, i, numel(allSubjects), numel(Lfiles), numel(Rfiles));

    subjectNotes = {};
    status = 'OK';

    if isempty(Lfiles)
        status = 'WARNING';
        subjectNotes{end+1} = 'Missing left mask; output uses right side only.'; %#ok<AGROW>
        issues = add_issue(issues, sid, 'WARNING', 'Missing left mask; using right side only.', '');
    end
    if isempty(Rfiles)
        status = 'WARNING';
        subjectNotes{end+1} = 'Missing right mask; output uses left side only.'; %#ok<AGROW>
        issues = add_issue(issues, sid, 'WARNING', 'Missing right mask; using left side only.', '');
    end
    if numel(Lfiles) > 1
        subjectNotes{end+1} = sprintf('Multiple left masks combined: %d.', numel(Lfiles)); %#ok<AGROW>
        issues = add_issue(issues, sid, 'INFO', sprintf('Multiple left masks combined: %d.', numel(Lfiles)), strjoin(Lfiles, ' | '));
    end
    if numel(Rfiles) > 1
        subjectNotes{end+1} = sprintf('Multiple right masks combined: %d.', numel(Rfiles)); %#ok<AGROW>
        issues = add_issue(issues, sid, 'INFO', sprintf('Multiple right masks combined: %d.', numel(Rfiles)), strjoin(Rfiles, ' | '));
    end

    try
        % Reference grid preference:
        %   - if left exists: use first readable left mask as reference
        %   - otherwise: use first readable right mask as reference
        Vref = [];

        if ~isempty(Lfiles)
            [Mleft, Vref, nVoxLeft, notesLeft, newIssues] = combine_masks_to_ref(Lfiles, Vref, thr, sid, 'left');
            issues = [issues; newIssues]; %#ok<AGROW>
            subjectNotes = [subjectNotes, notesLeft]; %#ok<AGROW>
        else
            Mleft = [];
            nVoxLeft = 0;
        end

        if ~isempty(Rfiles)
            [Mright, Vref, nVoxRight, notesRight, newIssues] = combine_masks_to_ref(Rfiles, Vref, thr, sid, 'right');
            issues = [issues; newIssues]; %#ok<AGROW>
            subjectNotes = [subjectNotes, notesRight]; %#ok<AGROW>
        else
            Mright = [];
            nVoxRight = 0;
        end

        if isempty(Vref)
            status = 'FAILED';
            subjectNotes{end+1} = 'No readable left or right mask.'; %#ok<AGROW>
            issues = add_issue(issues, sid, 'ERROR', 'No readable left or right mask.', '');
            qcRows(end+1,:) = {sid, '', numel(Lfiles), numel(Rfiles), nVoxLeft, nVoxRight, 0, status, strjoin(subjectNotes, ' ')}; %#ok<AGROW>
            fprintf('  FAILED: no readable masks.\n');
            continue
        end

        if isempty(Mleft)
            Mleft = false(Vref.dim);
        end
        if isempty(Mright)
            Mright = false(Vref.dim);
        end

        Mbilat = Mleft | Mright;
        nVoxBilateral = nnz(Mbilat);

        if nVoxBilateral < 1
            status = 'FAILED';
            subjectNotes{end+1} = 'Bilateral mask has zero voxels after thresholding.'; %#ok<AGROW>
            issues = add_issue(issues, sid, 'ERROR', 'Bilateral mask has zero voxels after thresholding.', '');
            qcRows(end+1,:) = {sid, '', numel(Lfiles), numel(Rfiles), nVoxLeft, nVoxRight, nVoxBilateral, status, strjoin(subjectNotes, ' ')}; %#ok<AGROW>
            fprintf('  FAILED: bilateral mask has zero voxels.\n');
            continue
        end

        outfile = fullfile(outdir, sprintf('%s_%s_bilateral.nii', sid, maskLabel));
        write_binary_nifti(uint8(Mbilat), Vref, outfile);

        qcRows(end+1,:) = {sid, outfile, numel(Lfiles), numel(Rfiles), nVoxLeft, nVoxRight, nVoxBilateral, status, strjoin(subjectNotes, ' ')}; %#ok<AGROW>

        fprintf('  Written: %s\n', outfile);
        fprintf('  Voxels: left=%d | right=%d | bilateral=%d | status=%s\n', ...
            nVoxLeft, nVoxRight, nVoxBilateral, status);

    catch ME
        status = 'FAILED';
        msg = sprintf('%s', ME.message);
        subjectNotes{end+1} = msg; %#ok<AGROW>
        issues = add_issue(issues, sid, 'ERROR', msg, '');
        qcRows(end+1,:) = {sid, '', numel(Lfiles), numel(Rfiles), NaN, NaN, NaN, status, strjoin(subjectNotes, ' ')}; %#ok<AGROW>
        warning('[%s] FAILED: %s', sid, msg);
    end
end

% ---------------- SAVE QC ----------------
qcHeader = {'SubjectID','OutputFile','nLeftFiles','nRightFiles','nVoxLeft','nVoxRight','nVoxBilateral','Status','Notes'};
writecell([qcHeader; qcRows], fullfile(outdir, 'BilateralMask_QC.csv'));

issueHeader = {'SubjectID','Severity','Issue','Files'};
writecell([issueHeader; issues], fullfile(outdir, 'BilateralMask_Issues.csv'));

fprintf('\nDONE.\n');
fprintf('Bilateral masks and QC files written to:\n  %s\n', outdir);
fprintf('Main QC file:\n  %s\n', fullfile(outdir, 'BilateralMask_QC.csv'));

end

%% =====================================================================
%% FILE SELECTION
%% =====================================================================
function files = select_multi_files(startDir, prompt, required, pattern)
if nargin < 4 || isempty(pattern)
    pattern = '*.*';
end
if isempty(startDir) || ~isfolder(startDir)
    startDir = pwd;
end

filter = fullfile(startDir, pattern);
[fs, path] = uigetfile(filter, prompt, 'MultiSelect', 'on');

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
for k = 1:numel(fs)
    files{k} = fullfile(path, fs{k});
end
end

%% =====================================================================
%% SUBJECT MATCHING
%% =====================================================================
function [M, selectedRows, issues] = group_files_by_subject(files, sideLabel, issues)
M = containers.Map('KeyType','char','ValueType','any');
selectedRows = {};

for i = 1:numel(files)
    f = files{i};
    sid = extract_subid(f);

    if isempty(sid)
        issues = add_issue(issues, '', 'ERROR', sprintf('Could not extract subject ID from %s mask.', sideLabel), f);
        continue
    end

    if isKey(M, sid)
        tmp = M(sid);
        tmp{end+1} = f; %#ok<AGROW>
        M(sid) = tmp;
    else
        M(sid) = {f};
    end

    selectedRows(end+1,:) = {sid, sideLabel, f}; %#ok<AGROW>
end
end

function sid = extract_subid(pathstr)
% Extracts sub-... until first separator-like character after the ID.
% Examples matched:
%   sub-933
%   sub-1372
%   sub-001
%   sub-ABC123, if present
m = regexp(char(pathstr), '(sub-[A-Za-z0-9]+)', 'match', 'once');
if isempty(m)
    sid = '';
else
    sid = m;
end
end

function sorted = sort_subject_ids(ids)
if isempty(ids)
    sorted = {};
    return
end

ids = unique(ids);
nums = nan(size(ids));
for i = 1:numel(ids)
    m = regexp(ids{i}, 'sub-(\d+)', 'tokens', 'once');
    if ~isempty(m)
        nums(i) = str2double(m{1});
    end
end

if all(isfinite(nums))
    [~, idx] = sort(nums);
else
    [~, idx] = sort(ids);
end
sorted = ids(idx);
end

function issues = add_issue(issues, sid, severity, issue, files)
if nargin < 5
    files = '';
end
issues(end+1,:) = {sid, severity, issue, files};
end

%% =====================================================================
%% MASK COMBINATION
%% =====================================================================
function [M, Vref, nVox, notes, issues] = combine_masks_to_ref(files, Vref, thr, sid, sideLabel)
issues = {};
notes = {};
M = [];
nVox = 0;

if isempty(files)
    return
end

% If no reference is provided, find the first readable mask and use its grid.
if isempty(Vref)
    Vref = [];
    for i = 1:numel(files)
        try
            VtmpAll = spm_vol(files{i});
            Vref = VtmpAll(1);
            break
        catch ME
            issues = add_issue(issues, sid, 'ERROR', sprintf('Could not read %s mask for reference: %s', sideLabel, ME.message), files{i});
        end
    end

    if isempty(Vref)
        issues = add_issue(issues, sid, 'ERROR', sprintf('No readable %s mask found.', sideLabel), strjoin(files, ' | '));
        return
    end
end

M = false(Vref.dim);

for i = 1:numel(files)
    f = files{i};
    try
        [Mi, wasResampled] = load_binary_mask_to_ref(f, Vref, thr);
        M = M | Mi;

        if wasResampled
            txt = sprintf('%s mask resampled to reference grid: %s', sideLabel, f);
            notes{end+1} = txt; %#ok<AGROW>
            issues = add_issue(issues, sid, 'INFO', txt, f);
        end
    catch ME
        issues = add_issue(issues, sid, 'ERROR', sprintf('Could not load %s mask: %s', sideLabel, ME.message), f);
    end
end

nVox = nnz(M);
if nVox < 1
    issues = add_issue(issues, sid, 'WARNING', sprintf('%s mask has zero voxels after thresholding.', sideLabel), strjoin(files, ' | '));
end
end

function [M, wasResampled] = load_binary_mask_to_ref(maskFile, Vref, thr)
VsrcAll = spm_vol(maskFile);
Vsrc = VsrcAll(1);

sameGrid = isequal(Vref.dim, Vsrc.dim) && max(abs(Vref.mat(:) - Vsrc.mat(:))) < 1e-6;
wasResampled = ~sameGrid;

if sameGrid
    M = spm_read_vols(Vsrc) > thr;
    M = logical(M);
    return
end

% Nearest-neighbor sampling into reference grid.
dim = Vref.dim;
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

function write_binary_nifti(M, Vref, outfile)
Vout = Vref;
Vout.fname = outfile;
Vout.dt = [spm_type('uint8') 0];
Vout.pinfo = [1; 0; 0];
Vout.descrip = 'Bilateral union mask created by CreateBilateralMasks_MatchBySubID';
spm_write_vol(Vout, uint8(M));
end
