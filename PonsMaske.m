 % ============================================================
% nm_native_pons_precise.m  (SPM12)
% ------------------------------------------------------------
% Präzise Pipeline:
%   MNI-Ponsmaske  --(iy_T1 via Pull, space=T1)-->  T1-Maske
%   NM  --(rigid->T1)-->  T1     ;  T1-Maske --(Pull, space=NM)--> NM-Maske
%   => exakte Geometrie (Affine, Voxel, BB) durch 'space' Zielbild
% ============================================================

clear; clc;
spm('defaults','fmri');
spm_jobman('initcfg');

%% --- Eingaben wählen ---
fprintf('\nWähle die MNI-Pons-Maske (z.B. MNI_PonsMask.nii)\n');
[mask_file, mask_path] = uigetfile({'*.nii;*.nii.gz','NIfTI'}, 'MNI-Maske auswählen');
if isequal(mask_file,0), return; end
mask_mni = fullfile(mask_path, mask_file);

fprintf('\nWähle alle T1-Dateien (Multi-Select)\n');
[t1_files, t1_path] = uigetfile({'*.nii;*.nii.gz','T1 NIfTI'}, 'T1 auswählen','MultiSelect','on');
if isequal(t1_files,0), return; end
if ischar(t1_files), t1_files = {t1_files}; end

fprintf('\nWähle alle NM-Dateien (Multi-Select)\n');
[nm_files, nm_path] = uigetfile({'*.nii;*.nii.gz','NM NIfTI'}, 'NM auswählen','MultiSelect','on');
if isequal(nm_files,0), return; end
if ischar(nm_files), nm_files = {nm_files}; end

out_dir = uigetdir(t1_path, 'Ausgabeordner wählen');
if isequal(out_dir,0), return; end
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% --- ID-Matching: sub-<Zahl> ---
get_id = @(f) regexp(f,'sub-(\d+)','tokens','once');
pairs = {};
for i = 1:numel(t1_files)
    id_t = get_id(t1_files{i});
    if isempty(id_t), warning('Keine ID in T1: %s', t1_files{i}); continue; end
    for j = 1:numel(nm_files)
        id_n = get_id(nm_files{j});
        if ~isempty(id_n) && strcmp(id_t{1}, id_n{1})
            pairs(end+1,:) = {id_t{1}, fullfile(t1_path,t1_files{i}), fullfile(nm_path,nm_files{j})};
        end
    end
end
if isempty(pairs), error('Keine passenden Paare (erwarte "sub-<ID>" in Namen).'); end
fprintf('\nGefundene Subjekte:\n'); disp(pairs(:,1));

%% --- Helfer: gunzip falls nötig ---
unzip_if_needed = @(p) ( ...
    (endsWith(lower(p),'.gz') && exist(p,'file') && ~exist(p(1:end-3),'file') && gunzip(p, fileparts(p))) ...
);

%% --- Hauptschleife ---
results = {'ID','MeanSignal','SD','VoxelCount'};
for k = 1:size(pairs,1)
    id = pairs{k,1}; T1 = pairs{k,2}; NM = pairs{k,3};
    fprintf('\n=== Subjekt %s ===\n', id);

    % ggf. entpacken
    if endsWith(lower(T1),'.gz'), unzip_if_needed(T1); T1 = T1(1:end-3); end
    if endsWith(lower(NM),'.gz'), unzip_if_needed(NM); NM = NM(1:end-3); end
    if endsWith(lower(mask_mni),'.gz')
        unzip_if_needed(mask_mni);
        mask_mni = mask_mni(1:end-3);
    end

    [t1_dir, t1_name, ~] = fileparts(T1);

    %% 1) Normalisierung T1->MNI (y_*)
    yT1 = fullfile(t1_dir, ['y_', t1_name, '.nii']);
    if ~exist(yT1,'file')
        fprintf('   Schätze Deformation (y_%s.nii)\n', t1_name);
        clear matlabbatch
        matlabbatch{1}.spm.spatial.normalise.est.subj.vol = {T1};
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.affreg = 'mni';
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.fwhm = 0;
        matlabbatch{1}.spm.spatial.normalise.est.eoptions.samp = 3;
        spm_jobman('run', matlabbatch);
    end

    %% 2) Inverse Deformation (iy_* oder y_iy_*)
    iyT1 = fullfile(t1_dir, ['iy_', t1_name, '.nii']);
    if ~exist(iyT1,'file')
        fprintf('   Erzeuge inverse Deformation (iy_%s.nii)\n', t1_name);
        clear matlabbatch
        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {yT1};
        matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {T1};
        matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = ['iy_', t1_name];
        matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {t1_dir};
        spm_jobman('run', matlabbatch);
    end
    if ~exist(iyT1,'file')
        cand = [dir(fullfile(t1_dir,'iy_*.nii')); dir(fullfile(t1_dir,'y_iy_*.nii')); ...
                dir(fullfile(pwd,'iy_*.nii'));     dir(fullfile(pwd,'y_iy_*.nii'))];
        if isempty(cand), error('Keine inverse Deformation gefunden.'); end
        [~,ix] = max([cand.datenum]); found = cand(ix);
        fprintf('   Gefunden: %s\n', found.name);
        iyT1 = fullfile(found.folder, found.name);
        std_name = fullfile(t1_dir, ['iy_', t1_name, '.nii']);
        if ~strcmp(iyT1, std_name), copyfile(iyT1, std_name); iyT1 = std_name; end
    end

    %% 3) MNI-Maske -> T1 **präzise** (Pull mit space = T1)
    fprintf('   Ziehe Maske (MNI->T1) mit exact space=T1 ...\n');
    t1_mask_out = fullfile(out_dir, sprintf('sub-%s_PonsMask_inT1.nii', id));
    clear matlabbatch
    % comp: nur inverse Deformation
    matlabbatch{1}.spm.util.defs.comp{1}.def = {iyT1};
    % pull: resample der MNI-Maske in EXAKTEN T1-Space
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {mask_mni};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {out_dir};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;    % nearest neighbour (labels sicher)
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'T1_';
    % entscheidend: Zielraum exakt T1
    matlabbatch{1}.spm.util.defs.comp{2}.id.space = {T1};
    spm_jobman('run', matlabbatch);

    % Finde erzeugte Datei (T1_-Präfix) & benenne um
    tmp = dir(fullfile(out_dir,'T1_*'));
    [~,ix] = max([tmp.datenum]); t1_mask_src = fullfile(out_dir, tmp(ix).name);
    if exist(t1_mask_out,'file'), delete(t1_mask_out); end
    movefile(t1_mask_src, t1_mask_out);

    %% 4) Coreg: NM -> T1 (rigid)
    fprintf('   Coregistriere NM -> T1 (rigid)\n');
    clear matlabbatch
    matlabbatch{1}.spm.spatial.coreg.estimate.ref = {T1};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {NM};
    matlabbatch{1}.spm.spatial.coreg.estimate.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = ...
        [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
    spm_jobman('run', matlabbatch);

    %% 5) T1-Maske -> NM **präzise** (Pull mit space = NM)
    fprintf('   Ziehe T1-Maske in NM-Space (exact space=NM) ...\n');
    nm_mask_out = fullfile(out_dir, sprintf('sub-%s_PonsMask_inNM.nii', id));
    clear matlabbatch
    % keine Deformation nötig (nur affiner Coreg von NM->T1 liegt jetzt in Headern),
    % wir resamplen die T1-Maske EXAKT in den NM-Space:
    matlabbatch{1}.spm.util.defs.comp{1}.id.space = {NM};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {t1_mask_out};
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {out_dir};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;  % labels sicher
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'NM_';
    spm_jobman('run', matlabbatch);

    tmp = dir(fullfile(out_dir,'NM_*'));
    [~,ix] = max([tmp.datenum]); nm_mask_src = fullfile(out_dir, tmp(ix).name);
    if exist(nm_mask_out,'file'), delete(nm_mask_out); end
    movefile(nm_mask_src, nm_mask_out);

    %% 6) Signal im NM (Maske binär/Label-tolerant)
    Vm = spm_vol(nm_mask_out); M = spm_read_vols(Vm);
    % falls deine Maske binär ist: M>0.5; falls Label: M>0 reicht ebenfalls (nearest!)
    mask = M > 0;
    Vn = spm_vol(NM); X = spm_read_vols(Vn);
    vals = X(mask);
    if isempty(vals)
        warning('Leere Maske im NM-Raum: %s', nm_mask_out);
        mu = NaN; sd = NaN; nv = 0;
    else
        mu = mean(vals(:),'omitnan'); sd = std(vals(:),'omitnan'); nv = nnz(mask);
    end
    fprintf('📊 Pons (NM-native): mean=%.6g  SD=%.6g  vox=%d\n', mu, sd, nv);
    results(end+1,:) = {id, mu, sd, nv};

    % Optional visueller Check:
    % spm_check_registration(char(NM, nm_mask_out));
end

%% CSV speichern
out_csv = fullfile(out_dir, 'NM_native_Pons_Signals.csv');
try
    writecell(results, out_csv);
    fprintf('\n✅ CSV gespeichert: %s\n', out_csv);
catch
    % Fallback für ältere MATLABs
    fid = fopen(out_csv,'w');
    fprintf(fid, 'ID,MeanSignal,SD,VoxelCount\n');
    for r = 2:size(results,1)
        fprintf(fid, '%s,%.6g,%.6g,%d\n', results{r,1}, results{r,2}, results{r,3}, results{r,4});
    end
    fclose(fid);
    fprintf('\n✅ CSV gespeichert (Fallback): %s\n', out_csv);
end

disp('✅ Fertig: Präzise Ponsmasken im NM-Raum & Signalwerte berechnet.');
