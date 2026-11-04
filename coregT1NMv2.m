%% ==========================================================
%% LC Coregistration (Option B, präzise) - BATCH
%% NM -> T1  (Maske sicher, NN + Transform-Übertragung)
%% ==========================================================
clear; clc;
spm('defaults','fmri'); spm_jobman('initcfg');
fprintf('\n🧠 START: LC Coregistration NM→T1 (Option B, präzise) [BATCH]\n\n');

%% ---------- 1) Dateien wählen (MultiSelect) ----------
fprintf('Bitte wähle die benötigten Dateien (Mehrfachauswahl möglich):\n');

% --- T1 ---
[fn_T1,fp_T1] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
                          'Wähle T1 (Referenz) - Mehrfachauswahl', ...
                          'MultiSelect','on');
if isequal(fn_T1,0), error('Abbruch (keine T1 gewählt).'); end
if ischar(fn_T1), fn_T1 = {fn_T1}; end   % auf Cell vereinheitlichen

% --- NM ---
[fn_NM,fp_NM] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
                          'Wähle Neuromelanin (Source) - Mehrfachauswahl', ...
                          'MultiSelect','on');
if isequal(fn_NM,0), error('Abbruch (kein NM gewählt).'); end
if ischar(fn_NM), fn_NM = {fn_NM}; end

% --- LC-Maske ---
[fn_LC,fp_LC] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
                          'Wähle LC-Masken (im NM-Raum, binär) - Mehrfachauswahl', ...
                          'MultiSelect','on');
if isequal(fn_LC,0), error('Abbruch (keine LC-Maske gewählt).'); end
if ischar(fn_LC), fn_LC = {fn_LC}; end

%% ---------- 2) Listen mit Subject-IDs bauen ----------
fprintf('\n📂 Erzeuge Listen mit Subject-IDs aus Dateinamen...\n');

T1_list = struct('sid',{},'fname',{});
for i = 1:numel(fn_T1)
    sid = get_subid_from_name(fn_T1{i});
    if isempty(sid)
        fprintf('  ⚠️  T1-Datei ohne erkennbaren Sub-Tag übersprungen: %s\n', fn_T1{i});
        continue;
    end
    T1_list(end+1).sid   = sid; %#ok<SAGROW>
    T1_list(end).fname   = fullfile(fp_T1, fn_T1{i});
end

NM_list = struct('sid',{},'fname',{});
for i = 1:numel(fn_NM)
    sid = get_subid_from_name(fn_NM{i});
    if isempty(sid)
        fprintf('  ⚠️  NM-Datei ohne erkennbaren Sub-Tag übersprungen: %s\n', fn_NM{i});
        continue;
    end
    NM_list(end+1).sid   = sid; %#ok<SAGROW>
    NM_list(end).fname   = fullfile(fp_NM, fn_NM{i});
end

LC_list = struct('sid',{},'fname',{});
for i = 1:numel(fn_LC)
    sid = get_subid_from_name(fn_LC{i});
    if isempty(sid)
        fprintf('  ⚠️  LC-Datei ohne erkennbaren Sub-Tag übersprungen: %s\n', fn_LC{i});
        continue;
    end
    LC_list(end+1).sid   = sid; %#ok<SAGROW>
    LC_list(end).fname   = fullfile(fp_LC, fn_LC{i});
end

if isempty(T1_list) || isempty(NM_list) || isempty(LC_list)
    error('Nach dem Filtern anhand der Subject-IDs sind eine oder mehrere Listen leer.');
end

%% ---------- 3) Subjects matchen ----------
all_T1_sid = {T1_list.sid};
all_NM_sid = {NM_list.sid};
all_LC_sid = {LC_list.sid};

common_sid = intersect(intersect(all_T1_sid, all_NM_sid), all_LC_sid);

if isempty(common_sid)
    error('Keine gemeinsamen Subject-IDs in T1, NM und LC-Masken gefunden.');
end

fprintf('\n✅ Gefundene gemeinsame Subjects:\n');
disp(common_sid(:));

%% ---------- 4) Loop über alle Subjects (Batch, ohne OK-Klicks) ----------
for s = 1:numel(common_sid)
    sid = common_sid{s};
    fprintf('\n=============================================\n');
    fprintf('▶ Verarbeite Subject: %s\n', sid);
    fprintf('=============================================\n');

    % --- passende Dateien für dieses Subject holen ---
    T1 = get_fname_for_sid(T1_list, sid);
    NM = get_fname_for_sid(NM_list, sid);
    LCmask = get_fname_for_sid(LC_list, sid);

    fprintf('T1:     %s\n', T1);
    fprintf('NM:     %s\n', NM);
    fprintf('LCmask: %s\n', LCmask);

    %% ---------- PRE-QC: Header prüfen (wie im Original) ----------
    V_NM = spm_vol(NM);
    V_LC = spm_vol(LCmask);

    same_dim = all(V_NM.dim == V_LC.dim);
    same_mat = all(abs(V_NM.mat(:) - V_LC.mat(:)) < 1e-6);
    if same_dim && same_mat
        fprintf('\n✅ Maske & NM haben identischen Header.\n');
    else
        fprintf('\n⚠️  Header unterschiedlich – kopiere NM-Header auf Maske (keine Interpolation)...\n');
        Y = spm_read_vols(V_LC);
        V_LC.mat   = V_NM.mat;
        V_LC.dim   = V_NM.dim;
        V_LC.pinfo = [1;0;0];
        V_LC.fname = spm_file(LCmask,'prefix',['geomfix_' sid '_']);
        spm_write_vol(V_LC,Y);
        LCmask = V_LC.fname;
        fprintf('   → Neue Maske: %s\n', LCmask);
    end

    % Optional: QC-Ansicht, aber OHNE msgbox/uiwait → blockiert nicht
    try
        spm_check_registration(char(T1,NM,LCmask));
        spm_orthviews('Caption',1,['T1 ' sid]);
        spm_orthviews('Caption',2,['NM ' sid]);
        spm_orthviews('Caption',3,['LC (NM) ' sid]);
    catch
        fprintf('   (QC-View konnte nicht geöffnet werden, fahre trotzdem fort.)\n');
    end

    %% ---------- Coregistration (NM -> T1, estimate only) ----------
    fprintf('\n=== Coreg (estimate only) NM -> T1 für %s ===\n', sid);

    % Alte NM-Matrix sichern
    V_NM_before = spm_vol(NM);
    M_before = V_NM_before.mat;

    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref    = {T1};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {NM};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2 1];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = ...
       [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
    spm_jobman('run',matlabbatch);
    fprintf('✅ Coreg estimate abgeschlossen (NM-Header aktualisiert) für %s.\n', sid);

    %% ---------- Transform auf LC-Maske übertragen ----------
    V_NM_after = spm_vol(NM);
    M_after = V_NM_after.mat;
    T = M_after / M_before;  % 4x4 starre Transform (NM→T1)

    fprintf('Übertrage NM→T1-Transform auf Maske (%s)...\n', sid);
    V_LC = spm_vol(LCmask);
    Y = spm_read_vols(V_LC);
    V_LC.mat   = T * V_LC.mat;      % ⬅️ entscheidender Schritt
    V_LC.fname = spm_file(LCmask,'prefix',['xform_' sid '_']);
    spm_write_vol(V_LC,Y);
    LCmask_T1hdr = V_LC.fname;
    fprintf('✅ Neue Maske mit Transform geschrieben: %s\n', LCmask_T1hdr);

    %% ---------- Reslice Maske (NN) ins T1-Gitter ----------
    fprintf('\n=== Reslice Maske (NN) in T1 für %s ===\n', sid);
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.write.ref    = {T1};
    matlabbatch{1}.spm.spatial.coreg.write.source = {LCmask_T1hdr};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; % Nearest Neighbour
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap   = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask   = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'cT1_';
    spm_jobman('run',matlabbatch);

    LC_inT1 = spm_file(LCmask_T1hdr,'prefix','cT1_');
    fprintf('✅ Fertig (%s): %s\n', sid, LC_inT1);

    %% ---------- QC nachher (ohne msgbox-Stop) ----------
    try
        spm_check_registration(char(T1,NM,LC_inT1));
        spm_orthviews('Caption',1,['T1 ' sid]);
        spm_orthviews('Caption',2,['NM ' sid]);
        spm_orthviews('Caption',3,['LC_inT1 ' sid]);
        title(['QC nach Coregistration - ' sid]);
    catch
        fprintf('   (QC-View nachher konnte nicht geöffnet werden.)\n');
    end

    Vmask = spm_vol(LC_inT1);
    M = spm_read_vols(Vmask) > 0.5;
    [idx_x,idx_y,idx_z] = ind2sub(size(M), find(M));
    if isempty(idx_x)
        warning('Maske enthält keine Voxel > 0.5 im T1-Raster! (%s)', sid);
    else
        com_mm = Vmask.mat * [mean(idx_x); mean(idx_y); mean(idx_z); 1];
        fprintf('\n[POST-QC %s] LC Center-of-Mass (mm): [%.2f %.2f %.2f]\n', ...
            sid, com_mm(1),com_mm(2),com_mm(3));
    end

    V_T1 = spm_vol(T1);
    same_pix = all(abs(diag(V_T1.mat(1:3,1:3)) - diag(Vmask.mat(1:3,1:3))) < 1e-6);
    if same_pix
        fprintf('[POST-QC %s] Voxelraster identisch mit T1: JA ✅\n', sid);
    else
        fprintf('[POST-QC %s] Voxelraster identisch mit T1: NEIN ⚠️\n', sid);
    end
    fprintf('[POST-QC %s] LC Voxels (binär): %d\n', sid, nnz(M));

    fprintf('\n🏁 DONE für %s.\n', sid);
end

fprintf('\n🎉 Alle Subjects abgearbeitet.\n');

%% ================== Helper-Funktionen ==================
function sid = get_subid_from_name(fname)
    % Subject-ID aus Dateinamen mit Pattern:
    % Sub-XXX, subXX, Sub_XX, Sub-0001, ...
    [~,name,~] = fileparts(fname);
    tok = regexp(name, '(?i)sub[-_]?(\d+)', 'tokens');
    if isempty(tok)
        sid = '';
    else
        sid = tok{1}{1};
    end
end

function fname = get_fname_for_sid(list, sid)
    % Findet den ersten Eintrag in list mit passender sid
    idx = find(strcmp({list.sid}, sid));
    if isempty(idx)
        error('Keine Datei für Subject %s gefunden.', sid);
    elseif numel(idx) > 1
        warning('Mehrere Dateien für Subject %s gefunden, nehme die erste.', sid);
    end
    fname = list(idx(1)).fname;
end
