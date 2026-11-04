% ============================================================
% Extract ALL ROI labels from a 3D label atlas (SPM required)
% Windows file selection included
% ============================================================

clc; clear;

% --- Atlas auswählen ---
[atlasFile, atlasPath] = uigetfile('*.nii', 'Select 3D label atlas NIfTI (e.g., atlas.nii)');
if isequal(atlasFile,0)
    error('No atlas selected.');
end
atlasFull = fullfile(atlasPath, atlasFile);

% --- Output-Ordner auswählen ---
outDir = uigetdir(pwd, 'Select output folder for ROI masks');
if isequal(outDir,0)
    error('No output folder selected.');
end

% --- Atlas laden ---
V = spm_vol(atlasFull);
if numel(V) > 1
    error('This file appears to be 4D. Script expects a 3D label atlas.');
end

Y = spm_read_vols(V);

% --- Alle Labels finden (0 und NaN entfernen) ---
vals = Y(:);
vals = vals(~isnan(vals));
labels = unique(vals);
labels = labels(labels ~= 0);          % 0 = background raus
labels = labels(labels == round(labels)); % nur Integer-Labels (Label-Atlas)

if isempty(labels)
    error('No non-zero integer labels found. Is this a label atlas?');
end

% --- Optional: Nachfrage ---
msg = sprintf('Found %d labels (non-zero integer). Extract ALL?', numel(labels));
choice = questdlg(msg, 'Extract all labels', 'Yes', 'No', 'Yes');
if ~strcmp(choice,'Yes')
    disp('Cancelled by user.');
    return;
end

fprintf('\nFound %d ROI labels. Extracting...\n', numel(labels));

% --- ROIs erzeugen ---
for k = 1:numel(labels)
    id = labels(k);

    mask = (Y == id);
    if ~any(mask(:))
        fprintf('Skipping label %d (empty)\n', id);
        continue;
    end

    Vo = V;
    Vo.fname = fullfile(outDir, sprintf('ROI_%04d.nii', id)); % 4-stellig, sauber sortierbar
    Vo.dt = [spm_type('uint8') 0];

    spm_write_vol(Vo, uint8(mask));

    if mod(k,50)==0
        fprintf('... %d/%d done\n', k, numel(labels));
    end
end

fprintf('\nDone. Saved %d ROI masks to:\n%s\n', numel(labels), outDir);
