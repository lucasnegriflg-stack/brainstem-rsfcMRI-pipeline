%% ==========================================================
%% SAFE & CLEAN (BATCH): EPI (4D, preprocessed) → T1 Raum
%% erzeugt meanEPI aus Realign+Unwarp+STC EPI
%% Batch über mehrere Subjects mit Sub-ID-Matching
%% ==========================================================
clear; clc;
spm('defaults','fmri'); spm_jobman('initcfg');
fprintf('\n🧠 START: 4D EPI (realigned) → T1 Raum (mit meanEPI) [BATCH]\n\n');

%% ---------- 1) Dateien wählen (MultiSelect) ----------
fprintf('Bitte wähle die Dateien (Mehrfachauswahl möglich):\n');

% --- T1 ---
[fn_T1,fp_T1] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
    'Wähle T1-Strukturbilder (Referenz) - Mehrfachauswahl', ...
    'MultiSelect','on');
if isequal(fn_T1,0), error('Abbruch (keine T1 gewählt).'); end
if ischar(fn_T1), fn_T1 = {fn_T1}; end   % auf Cell vereinheitlichen

% --- 4D-EPI ---
[fn_EPI,fp_EPI] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
    'Wähle 4D-EPI (nach Realign/Unwarp/STC aus CONN) - Mehrfachauswahl', ...
    'MultiSelect','on');
if isequal(fn_EPI,0), error('Abbruch (kein EPI gewählt).'); end
if ischar(fn_EPI), fn_EPI = {fn_EPI}; end

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

EPI_list = struct('sid',{},'fname',{});
for i = 1:numel(fn_EPI)
    sid = get_subid_from_name(fn_EPI{i});
    if isempty(sid)
        fprintf('  ⚠️  EPI-Datei ohne erkennbaren Sub-Tag übersprungen: %s\n', fn_EPI{i});
        continue;
    end
    EPI_list(end+1).sid   = sid; %#ok<SAGROW>
    EPI_list(end).fname   = fullfile(fp_EPI, fn_EPI{i});
end

if isempty(T1_list) || isempty(EPI_list)
    error('Nach dem Filtern anhand der Subject-IDs ist eine der Listen leer.');
end

%% ---------- 3) Subjects matchen ----------
all_T1_sid  = {T1_list.sid};
all_EPI_sid = {EPI_list.sid};

common_sid = intersect(all_T1_sid, all_EPI_sid);

if isempty(common_sid)
    error('Keine gemeinsamen Subject-IDs in T1 und EPI gefunden.');
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
    T1    = get_fname_for_sid(T1_list,  sid);
    EPI4D = get_fname_for_sid(EPI_list, sid);

    fprintf('T1:    %s\n', T1);
    fprintf('EPI4D: %s\n', EPI4D);

    %% ---------- 2) meanEPI erstellen ----------
    fprintf('\n=== Erzeuge meanEPI aus 4D-EPI (%s) ===\n', sid);
    V_epi = spm_vol(EPI4D);
    nVols = numel(V_epi);

    fprintf('  %d Volumes gefunden. Berechne Mittelwert...\n', nVols);
    Ysum = [];
    for i = 1:nVols
        Y = spm_read_vols(V_epi(i));
        if i==1
            Ysum = zeros(size(Y), 'single');
        end
        Ysum = Ysum + single(Y);
    end
    Ymean = Ysum / nVols;

    Vmean = V_epi(1);
    Vmean.fname = spm_file(EPI4D,'prefix','mean_');
    spm_write_vol(Vmean, Ymean);
    meanEPI = Vmean.fname;
    fprintf('✅ meanEPI geschrieben: %s\n', meanEPI);

    %% ---------- 3) QC vorher (optional, aber ohne msgbox) ----------
    try
        spm_check_registration(char(T1, meanEPI));
        spm_orthviews('Caption',1,['T1 ' sid]);
        spm_orthviews('Caption',2,['meanEPI ' sid]);
    catch
        fprintf('   (QC-View vor Coreg konnte nicht geöffnet werden, fahre fort.)\n');
    end

    %% ---------- 4) Coregistration: meanEPI -> T1 (estimate only) ----------
    fprintf('\n=== Coreg (estimate only) meanEPI -> T1 (%s) ===\n', sid);
    V_before = spm_vol(meanEPI);
    M_before = V_before.mat;

    matlabbatch = [];
    matlabbatch{1}.spm.spatial.coreg.estimate.ref    = {T1};
    matlabbatch{1}.spm.spatial.coreg.estimate.source = {meanEPI};
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep      = [4 2 1];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol      = ...
     [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm     = [7 7];
    spm_jobman('run', matlabbatch);

    V_after = spm_vol(meanEPI);
    M_after = V_after.mat;
    T = M_after / M_before;   % rigid Transform (meanEPI_before → T1)
    fprintf('✅ Coreg-Schätzung ok. (Rigid Transform berechnet für %s)\n', sid);

    %% ---------- 5) Apply Transform: 4D-EPI → T1 Raum (linear) ----------
    fprintf('\n=== Wende Transform auf 4D-EPI an (%s) ===\n', sid);

    V_epi = spm_vol(EPI4D);
    nVols = numel(V_epi);
    fprintf('  Anzahl Volumes: %d\n', nVols);

    V_new = V_epi;  % vorinitialisieren
    for i = 1:nVols
        V_new(i)      = V_epi(i);
        V_new(i).mat  = T * V_epi(i).mat;
        V_new(i).fname = spm_file(V_epi(i).fname,'prefix','cT1_');
    end

    spm_progress_bar('Init', nVols, ['Header-Transform EPI ' sid], '');
    for i = 1:nVols
        Y = spm_read_vols(V_epi(i));
        spm_write_vol(V_new(i), Y);
        spm_progress_bar('Set', i);
    end
    spm_progress_bar('Clear');

    EPI_inT1 = spm_file(EPI4D,'prefix','cT1_');
    fprintf('✅ Fertig (%s): %s (4D-EPI jetzt im T1-Raum)\n', sid, EPI_inT1);

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
        sid = tok{1}{1};   % nur die Nummer, z.B. '01', '001', '123'
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
