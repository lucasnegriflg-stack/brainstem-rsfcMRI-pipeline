%% simple_reslice_sn_masks_to_epi_by_subid.m
% Simple version:
%   - Select EPI files with MultiSelect
%   - Select SN mask files with MultiSelect
%   - Match by subject number after "sub"
%       Sub-123, sub123, sub_123, sub-001 all work
%   - No coregistration is estimated
%   - Masks are only resliced to the EPI voxel grid
%   - Uses nearest-neighbour interpolation
%
% Requirement:
%   SPM12 must be on the MATLAB path.

clear; clc;

%% --- SPM setup ---
if exist('spm', 'file') ~= 2
    spmDir = uigetdir(pwd, 'Select SPM12 folder');
    if isequal(spmDir, 0)
        error('SPM12 not selected.');
    end
    addpath(spmDir);
end

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

%% --- Select files ---
[epiNames, epiPath] = uigetfile( ...
    {'*.nii;*.nii.gz', 'NIfTI files (*.nii, *.nii.gz)'}, ...
    'Select ALL EPI reference files', ...
    'MultiSelect', 'on');

if isequal(epiNames, 0)
    error('No EPI files selected.');
end

[maskNames, maskPath] = uigetfile( ...
    {'*.nii;*.nii.gz', 'NIfTI files (*.nii, *.nii.gz)'}, ...
    'Select ALL SN mask files', ...
    'MultiSelect', 'on');

if isequal(maskNames, 0)
    error('No mask files selected.');
end

outDir = uigetdir(pwd, 'Select output folder');

if isequal(outDir, 0)
    error('No output folder selected.');
end

epiNames  = makeCell(epiNames);
maskNames = makeCell(maskNames);

epiFiles = cell(size(epiNames));
for i = 1:numel(epiNames)
    epiFiles{i} = fullfile(epiPath, epiNames{i});
end

maskFiles = cell(size(maskNames));
for i = 1:numel(maskNames)
    maskFiles{i} = fullfile(maskPath, maskNames{i});
end

%% --- Extract subject numbers ---
epiSub = cell(numel(epiFiles), 1);
maskSub = cell(numel(maskFiles), 1);

fprintf('\nEPI subject IDs:\n');
for i = 1:numel(epiFiles)
    [epiSub{i}, ok] = getSubNumber(epiFiles{i});
    if ok
        fprintf('  %s  ->  sub-%s\n', epiNames{i}, epiSub{i});
    else
        fprintf('  %s  ->  NO SUBJECT ID FOUND\n', epiNames{i});
    end
end

fprintf('\nMask subject IDs:\n');
for i = 1:numel(maskFiles)
    [maskSub{i}, ok] = getSubNumber(maskFiles{i});
    if ok
        fprintf('  %s  ->  sub-%s\n', maskNames{i}, maskSub{i});
    else
        fprintf('  %s  ->  NO SUBJECT ID FOUND\n', maskNames{i});
    end
end

%% --- Reslice matching masks to matching EPI grid ---
fprintf('\nStarting reslicing...\n');

nDone = 0;
nSkipped = 0;

for e = 1:numel(epiFiles)

    sid = epiSub{e};

    if isempty(sid)
        warning('Skipping EPI without subject ID: %s', epiFiles{e});
        nSkipped = nSkipped + 1;
        continue;
    end

    matchingMasks = find(strcmp(maskSub, sid));

    if isempty(matchingMasks)
        warning('No masks found for sub-%s', sid);
        nSkipped = nSkipped + 1;
        continue;
    end

    subjOut = fullfile(outDir, ['sub-' sid]);
    if ~exist(subjOut, 'dir')
        mkdir(subjOut);
    end

    % Copy/unzip EPI into output folder so SPM can read it safely.
    epiCopy = copyOrGunzipForSpm(epiFiles{e}, subjOut, ['REF_sub-' sid]);

    fprintf('\n============================================================\n');
    fprintf('Subject sub-%s\n', sid);
    fprintf('EPI: %s\n', epiFiles{e});
    fprintf('Matching masks: %d\n', numel(matchingMasks));
    fprintf('Output folder: %s\n', subjOut);
    fprintf('============================================================\n');

    for k = 1:numel(matchingMasks)

        m = matchingMasks(k);

        fprintf('\nReslicing mask:\n%s\n', maskFiles{m});

        % Copy/unzip mask into subject output folder.
        maskCopy = copyOrGunzipForSpm(maskFiles{m}, subjOut, 'MASK');

        % Reslice copied mask to EPI grid.
        flags = struct;
        flags.interp = 0;          % nearest neighbour: correct for masks
        flags.wrap   = [0 0 0];
        flags.mask   = 0;
        flags.which  = 1;          % reslice only second image, not reference
        flags.mean   = 0;
        flags.prefix = 'EPIgrid_';

        P = char([epiCopy ',1'], [maskCopy ',1']);
        spm_reslice(P, flags);

        [p, n, ext] = fileparts(maskCopy);
        reslicedMask = fullfile(p, ['EPIgrid_' n ext]);

        if ~exist(reslicedMask, 'file')
            warning('Expected resliced mask not found: %s', reslicedMask);
            nSkipped = nSkipped + 1;
            continue;
        end

        % Binarize final mask.
        V = spm_vol(reslicedMask);
        Y = spm_read_vols(V);
        Y(~isfinite(Y)) = 0;
        Ybin = uint8(Y > 0.5);

        [~, originalMaskName] = fileparts(stripNiiGz(maskFiles{m}));
        finalMask = fullfile(subjOut, ['sub-' sid '_' originalMaskName '_EPIgrid_bin.nii']);

        Vout = V;
        Vout.fname = finalMask;
        Vout.dt = [spm_type('uint8') 0];
        Vout.pinfo = [1; 0; 0];

        spm_write_vol(Vout, Ybin);

        fprintf('Saved final mask:\n%s\n', finalMask);

        nDone = nDone + 1;
    end
end

fprintf('\n============================================================\n');
fprintf('Done.\n');
fprintf('Final masks created: %d\n', nDone);
fprintf('Skipped items: %d\n', nSkipped);
fprintf('Output root folder:\n%s\n', outDir);
fprintf('============================================================\n');

%% ========================================================================
%% Local helper functions
%% ========================================================================

function c = makeCell(x)
    if ischar(x)
        c = {x};
    else
        c = x;
    end
end

function [sid, ok] = getSubNumber(filePath)
    % Robustly extracts the number after "sub".
    %
    % Examples:
    %   Sub-123_something.nii  -> 123
    %   sub123_mask.nii       -> 123
    %   sub_123_epi.nii       -> 123
    %   SUB 001_mask.nii      -> 1
    %
    % This intentionally ignores case and separator style.

    [~, name, ext] = fileparts(filePath);

    if strcmpi(ext, '.gz')
        [~, name2, ext2] = fileparts(name);
        base = [name2 ext2 ext];
    else
        base = [name ext];
    end

    base = lower(base);

    tok = regexp(base, 'sub[-_ ]*0*([0-9]+)', 'tokens', 'once');

    if isempty(tok)
        sid = '';
        ok = false;
    else
        sid = num2str(str2double(tok{1}));
        ok = true;
    end
end

function outFile = copyOrGunzipForSpm(inFile, outDir, prefix)
    % Copies .nii or unzips .nii.gz to outDir.
    % Returns an uncompressed .nii path for SPM.

    [~, name, ext] = fileparts(inFile);
    prefix = regexprep(prefix, '[^A-Za-z0-9_\-]', '_');

    if strcmpi(ext, '.gz')
        gunzip(inFile, outDir);

        % For something.nii.gz, fileparts returns name = something.nii.
        unzipped = fullfile(outDir, name);

        [~, n2, e2] = fileparts(unzipped);
        outFile = fullfile(outDir, [prefix '_' n2 e2]);

        if exist(outFile, 'file')
            delete(outFile);
        end

        movefile(unzipped, outFile);

    elseif strcmpi(ext, '.nii')
        outFile = fullfile(outDir, [prefix '_' name ext]);

        if exist(outFile, 'file')
            delete(outFile);
        end

        copyfile(inFile, outFile);

    else
        error('Only .nii and .nii.gz are supported in this simple script: %s', inFile);
    end
end

function cleanPath = stripNiiGz(filePath)
    % Returns path without .nii or .nii.gz extension for output naming.

    [p, name, ext] = fileparts(filePath);

    if strcmpi(ext, '.gz')
        [~, name2, ext2] = fileparts(name);
        if strcmpi(ext2, '.nii')
            cleanPath = fullfile(p, name2);
        else
            cleanPath = fullfile(p, name);
        end
    elseif strcmpi(ext, '.nii')
        cleanPath = fullfile(p, name);
    else
        cleanPath = fullfile(p, name);
    end
end
