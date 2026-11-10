%% ==========================================================
%% LC COM → EPI-Seeds (L/R getrennt, Liquor-Cropping)
%% Batch per Klick (4x MultiSelect): erst alle LC_L, dann alle LC_R, dann alle c3, dann alle EPI
%% + Output-Ordner frei wählbar
%% + Seed-Radius frei einstellbar
%% + QC-Logging (PASS/FAIL) ohne Skript-Abbruch
%%
%% Erwartete Inputs (je Subject-ID konsistent):
%%   - LC links  (T1-Raum): ..._L_...LC.nii
%%   - LC rechts (T1-Raum): ..._R_...LC.nii
%%   - CSF/c3 (T1 oder EPI-grid): c3...sub-<id>...nii  (oder enthält "c3")
%%   - EPI 4D (T1-aligned, nicht normalisiert): ...epi.nii
%%
%% Output:
%%   - Seeds (nur wenn PASS für die jeweilige Seite):
%%       LC_L_seed<r>mm_Sub<id>_<EPIname>.nii
%%       LC_R_seed<r>mm_Sub<id>_<EPIname>.nii
%%   - QC Log:
%%       LC_seed_QC_log.csv  (im Output-Ordner)
%% ==========================================================

clear; clc;
spm('defaults','fmri'); spm_jobman('initcfg');
fprintf('\n START: LC COM (L/R) → Seeds im EPI-Raum (mit Liquor-Cropping)\n\n');

%% ---------------- Parameter UI ----------------
dlg = inputdlg( ...
    {'Seed-Radius in mm (z.B. 2 oder 3):', 'CSF/c3 Threshold (0..1, z.B. 0.7):'}, ...
    'Parameter', [1 45], {'2','0.8'});

if isempty(dlg), error('Abbruch (keine Parameter gewählt).'); end

radius_mm = str2double(dlg{1});
thr       = str2double(dlg{2});

if ~isfinite(radius_mm) || radius_mm <= 0
    error('Ungültiger Seed-Radius. Bitte > 0 eingeben.');
end
if ~isfinite(thr) || thr < 0 || thr > 1
    error('Ungültiger CSF-Threshold. Bitte Wert zwischen 0 und 1 eingeben.');
end

outdir = uigetdir(pwd, 'Wähle Output-Ordner (Seeds + QC-Log)');
if isequal(outdir,0), error('Abbruch (kein Output-Ordner gewählt).'); end

fprintf('Parameter: radius=%.2f mm | CSF thr=%.2f | Output=%s\n\n', radius_mm, thr, outdir);

%% ---------------- 1) Dateien wählen (4x MultiSelect) ----------------
fprintf('Bitte wähle nun die Dateien in 4 Schritten (MultiSelect):\n');

LC_L_files = local_uigetfiles('Wähle ALLE LINKEN LC-Masken (LC_L) im T1-Raum', {'*.nii;*.nii.gz','NIfTI'});
LC_R_files = local_uigetfiles('Wähle ALLE RECHTEN LC-Masken (LC_R) im T1-Raum', {'*.nii;*.nii.gz','NIfTI'});
C3_files   = local_uigetfiles('Wähle ALLE c3/CSF-Masken (T1 oder EPI-grid)',       {'*.nii;*.nii.gz','NIfTI'});
EPI_files  = local_uigetfiles('Wähle ALLE 4D-EPI-Dateien (T1-aligned)',           {'*.nii;*.nii.gz','NIfTI'});

%% ---------------- 2) Subject-IDs extrahieren + Matchen ----------------
get_subid = @(s) local_get_subid(s);

LC_L_map = local_build_map(LC_L_files, get_subid, 'LC_L');
LC_R_map = local_build_map(LC_R_files, get_subid, 'LC_R');
C3_map   = local_build_map(C3_files,   get_subid, 'c3/CSF');
EPI_map  = local_build_map(EPI_files,  get_subid, 'EPI');

% Welche Subjects können verarbeitet werden? -> Schnittmenge aller 4 Kategorien
subs = intersect(intersect(keys(LC_L_map), keys(LC_R_map)), intersect(keys(C3_map), keys(EPI_map)));
subs = sort(string(subs));

if isempty(subs)
    fprintf('\n❌ Keine Subject-ID wurde in ALLEN 4 Listen gefunden.\n');
    fprintf('Tipp: Prüfe Dateinamen (Sub933 / sub-933 / usub-933 / xform_933) und ob pro Subject je 1 Datei je Kategorie gewählt ist.\n');
    return;
end

fprintf(' Gefundene Subjects mit vollständigem Set: %d\n', numel(subs));
disp(subs');

%% ---------------- 3) QC Log vorbereiten ----------------
logFile = fullfile(outdir, 'LC_seed_QC_log.csv');
writeHeader = ~exist(logFile, 'file');

fid = fopen(logFile, 'a');
if fid < 0
    error('Konnte QC-Logfile nicht öffnen: %s', logFile);
end
if writeHeader
    fprintf(fid, ['timestamp,subid,epi,lc_l,lc_r,csf_in,csf_epi,', ...
                  'radius_mm,csf_thr,', ...
                  'L_before,L_after,R_before,R_after,', ...
                  'fail_L,fail_R,status,note\n']);
end
fclose(fid);

%% ---------------- 4) Loop über Subjects ----------------
fprintf('\n=== Verarbeitung startet ===\n');

for si = 1:numel(subs)
    SUBID = subs(si);
    fprintf('\n---------------------------------------------\n');
    fprintf('Subject %s\n', SUBID);
    fprintf('---------------------------------------------\n');

    LC_L_T1 = LC_L_map(char(SUBID));
    LC_R_T1 = LC_R_map(char(SUBID));
    LIQ_IN  = C3_map(char(SUBID));
    EPI_T1  = EPI_map(char(SUBID));

    fprintf('LC_L: %s\nLC_R: %s\nc3:   %s\nEPI:  %s\n', LC_L_T1, LC_R_T1, LIQ_IN, EPI_T1);

    % ---------- COM berechnen ----------
    try
        com_L_mm = local_calc_com_mm(LC_L_T1);
        com_R_mm = local_calc_com_mm(LC_R_T1);
    catch ME
        local_log_qc(logFile, SUBID, EPI_T1, LC_L_T1, LC_R_T1, LIQ_IN, "", radius_mm, thr, ...
            NaN, NaN, NaN, NaN, 1, 1, "FAIL", "COM calc failed: "+string(ME.message));
        warning('COM-Berechnung fehlgeschlagen (%s). Subject wird übersprungen.', ME.message);
        continue;
    end

    % ---------- EPI Header ----------
    try
        Vepi_all = spm_vol(EPI_T1);
        Vref     = Vepi_all(1);
        M_epi    = Vref.mat;
        dims     = Vref.dim;
    catch ME
        local_log_qc(logFile, SUBID, EPI_T1, LC_L_T1, LC_R_T1, LIQ_IN, "", radius_mm, thr, ...
            NaN, NaN, NaN, NaN, 1, 1, "FAIL", "EPI read failed: "+string(ME.message));
        warning('EPI-Header konnte nicht gelesen werden (%s). Subject wird übersprungen.', ME.message);
        continue;
    end

    % ---------- CSF ins EPI-grid reslicen falls nötig ----------
    LIQ_EPI = LIQ_IN;
    try
        Vliq = spm_vol(LIQ_IN);
        sameGrid = isequal(Vliq.dim, dims) && max(abs(Vliq.mat(:) - M_epi(:))) < 1e-4;

        if ~sameGrid
            flags = struct('interp',1,'wrap',[0 0 0],'mask',0,'which',1,'mean',0);
            spm_reslice(char(Vref.fname, LIQ_IN), flags);
            LIQ_EPI = spm_file(LIQ_IN,'prefix','r'); % r*
        end

        Vliq = spm_vol(LIQ_EPI);
        Yliq = spm_read_vols(Vliq);
        liq_bin = Yliq > thr;

    catch ME
        local_log_qc(logFile, SUBID, EPI_T1, LC_L_T1, LC_R_T1, LIQ_IN, LIQ_EPI, radius_mm, thr, ...
            NaN, NaN, NaN, NaN, 1, 1, "FAIL", "CSF read/reslice failed: "+string(ME.message));
        warning('CSF Reslice/Read fehlgeschlagen (%s). Subject wird übersprungen.', ME.message);
        continue;
    end

    % ---------- Koordinatengitter ----------
    [Y,X,Z]  = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
    vox_grid = [Y(:) X(:) Z(:) ones(numel(Y),1)]';
    mm_grid  = M_epi * vox_grid;
    mm_grid  = mm_grid(1:3,:);

    % ---------- 2mm/Custom Kugeln ----------
    dL = sqrt(sum( (mm_grid - com_L_mm).^2, 1));
    mask_L_sphere = reshape(dL <= radius_mm, dims);

    dR = sqrt(sum( (mm_grid - com_R_mm).^2, 1));
    mask_R_sphere = reshape(dR <= radius_mm, dims);

    nL_before = nnz(mask_L_sphere);
    nR_before = nnz(mask_R_sphere);

    % ---------- Cropping ----------
    mask_L = mask_L_sphere & ~liq_bin;
    mask_R = mask_R_sphere & ~liq_bin;

    nL_after = nnz(mask_L);
    nR_after = nnz(mask_R);

    qc_fail_L = (nL_after == 0);
    qc_fail_R = (nR_after == 0);

    status = "PASS";
    note   = "ok";
    if qc_fail_L || qc_fail_R
        status = "FAIL";
        note = "seed fully removed by CSF-cropping";
    end

    fprintf('Voxel: L %d→%d | R %d→%d | status=%s\n', nL_before, nL_after, nR_before, nR_after, status);

    % ---------- Schreiben (nur pro Seite wenn nicht fail) ----------
    seed_L_fname = "";
    seed_R_fname = "";

    [~, epiname, ~] = fileparts(EPI_T1);

    if ~qc_fail_L
        V_L = Vref;
        V_L.fname   = fullfile(outdir, sprintf('LC_L_seed%gmm_Sub%s_%s.nii', radius_mm, char(SUBID), epiname));
        V_L.descrip = sprintf('LC left seed (%.2fmm sphere, COM-based, CSF-cropped thr=%.2f)', radius_mm, thr);
        spm_write_vol(V_L, double(mask_L));
        seed_L_fname = V_L.fname;
        fprintf('   wrote L seed: %s\n', seed_L_fname);
    else
        fprintf('   L seed not written (0 vox after cropping)\n');
    end

    if ~qc_fail_R
        V_R = Vref;
        V_R.fname   = fullfile(outdir, sprintf('LC_R_seed%gmm_Sub%s_%s.nii', radius_mm, char(SUBID), epiname));
        V_R.descrip = sprintf('LC right seed (%.2fmm sphere, COM-based, CSF-cropped thr=%.2f)', radius_mm, thr);
        spm_write_vol(V_R, double(mask_R));
        seed_R_fname = V_R.fname;
        fprintf('   wrote R seed: %s\n', seed_R_fname);
    else
        fprintf('   R seed not written (0 vox after cropping)\n');
    end

    % ---------- QC Log (immer) ----------
    local_log_qc(logFile, SUBID, EPI_T1, LC_L_T1, LC_R_T1, LIQ_IN, LIQ_EPI, radius_mm, thr, ...
        nL_before, nL_after, nR_before, nR_after, qc_fail_L, qc_fail_R, status, note);

    % ---------- Quick-QC Overlay (optional; zeigt nur existierende Seeds) ----------
    try
        toShow = {Vref.fname};
        if seed_L_fname ~= "" && exist(seed_L_fname,'file'), toShow{end+1} = seed_L_fname; end %#ok<AGROW>
        if seed_R_fname ~= "" && exist(seed_R_fname,'file'), toShow{end+1} = seed_R_fname; end %#ok<AGROW>

        spm_check_registration(char(toShow{:}));
        spm_orthviews('Caption',1,'EPI (T1-aligned)');
        capIdx = 2;
        if seed_L_fname ~= "", spm_orthviews('Caption',capIdx,'LC_L seed'); capIdx = capIdx + 1; end
        if seed_R_fname ~= "", spm_orthviews('Caption',capIdx,'LC_R seed'); end

        if status == "FAIL"
            uiwait(msgbox(sprintf('QC-FAIL geloggt.\nL_after=%d, R_after=%d\nSub=%s', nL_after, nR_after, SUBID), ...
                'LC Seed QC (FAIL)'));
        else
            uiwait(msgbox(sprintf('QC-PASS.\nSub=%s', SUBID), 'LC Seed QC (PASS)'));
        end
    catch
        fprintf('  (QC-Ansicht konnte nicht geöffnet werden – ggf. SPM-GUI nicht aktiv.)\n');
    end
end

fprintf('\n🏁 DONE. QC-Log: %s\n', logFile);

%% ==========================================================
%% -------------------- Helper functions --------------------
%% ==========================================================

function files = local_uigetfiles(titleStr, filterSpec)
    [fn, fp] = uigetfile(filterSpec, titleStr, 'MultiSelect','on');
    if isequal(fn,0)
        error('Abbruch (keine Dateien gewählt) bei: %s', titleStr);
    end
    if ischar(fn), fn = {fn}; end
    files = fullfile(fp, fn);
end

function m = local_build_map(files, get_subid, label)
    m = containers.Map('KeyType','char','ValueType','char');
    for i=1:numel(files)
        sid = get_subid(files{i});
        if sid == ""
            warning('%s: Datei ohne erkennbaren SubID-Token ignoriert: %s', label, files{i});
            continue;
        end
        key = char(sid);
        if isKey(m, key)
            warning('%s: Mehrfacher SubID=%s. Verwende LETZTE Datei: %s', label, key, files{i});
        end
        m(key) = files{i};
    end
    if isempty(keys(m))
        error('%s: Keine Dateien mit erkennbarer Subject-ID.', label);
    end
end

function com_mm = local_calc_com_mm(fname)
    V = spm_vol(fname);
    Y = spm_read_vols(V);
    M = Y > 0.5;
    idx = find(M);
    if isempty(idx)
        error('Maske enthält keine Voxel >0.5: %s', fname);
    end
    [ix,iy,iz] = ind2sub(size(M), idx);
    vox = [ix iy iz ones(numel(ix),1)]';
    mm  = V.mat * vox;
    mm  = mm(1:3,:);
    com_mm = mean(mm,2);
end

function subid = local_get_subid(fname)
    s = lower(fname);

    % 1) sub-933 / sub933
    tok = regexp(s, 'sub[-_ ]?(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        subid = string(tok{1}); return;
    end

    % 2) usub-933
    tok = regexp(s, 'usub[-_ ]?(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        subid = string(tok{1}); return;
    end

    % 3) xform_933
    tok = regexp(s, 'xform[-_ ]?(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        subid = string(tok{1}); return;
    end

    subid = "";
end

function local_log_qc(logFile, subid, epi, lc_l, lc_r, csf_in, csf_epi, radius_mm, thr, ...
                      nL_before, nL_after, nR_before, nR_after, fail_L, fail_R, status, note)
    fid = fopen(logFile, 'a');
    if fid < 0
        warning('Konnte QC-Logfile nicht öffnen: %s', logFile);
        return;
    end
    ts = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

    % NaN-safe printing for voxel counts
    if isnan(nL_before), nL_before = -1; end
    if isnan(nL_after),  nL_after  = -1; end
    if isnan(nR_before), nR_before = -1; end
    if isnan(nR_after),  nR_after  = -1; end

    fprintf(fid, '%s,%s,"%s","%s","%s","%s","%s",%.2f,%.2f,%d,%d,%d,%d,%d,%d,%s,"%s"\n', ...
        ts, char(subid), epi, lc_l, lc_r, csf_in, csf_epi, ...
        radius_mm, thr, ...
        nL_before, nL_after, nR_before, nR_after, ...
        logical(fail_L), logical(fail_R), char(status), char(note));
    fclose(fid);
end
