function warp_drn_and_exclude_csf_subjectmatch_spm12_win()
% SPM12 (Windows GUI)
%
% Warp MULTIPLE MNI seed masks into each subject's EPI space
% using inverse deformation fields, then exclude CSF using
% subject-specific c3 maps.
%
% Seed selection:
%   Seeds are selected ONE AT A TIME.
%   After each selected seed, the script asks whether another
%   seed should be added.
%
% Inputs (GUI):
%   - Multiple MNI seed masks (.nii/.nii.gz), selected sequentially
%   - Multiple reference EPIs (e.g., meanEPI*.nii) [MultiSelect]
%   - Multiple deformation fields y_*.nii (T1->MNI) [MultiSelect]
%   - Multiple CSF maps c3*.nii (T1 space TPM) [MultiSelect]
%   - Output directory
%
% Subject matching:
%   - Matches by token "sub-<ID>" in filenames
%   - Letters/digits allowed
%
% Example outputs:
%   sub-001_DRN_inEPI.nii
%   sub-001_DRN_inEPI_CSFexcluded.nii
%
%   sub-001_LC_inEPI.nii
%   sub-001_LC_inEPI_CSFexcluded.nii
%
% Plus:
%   CSV log containing one row per subject x seed combination.


% ============================================================
% INIT
% ============================================================

spm('Defaults','fMRI');
spm_jobman('initcfg');


% ============================================================
% SELECT SEEDS ONE BY ONE
% ============================================================

seed_full   = {};
seed_labels = {};

fprintf('\n====================================================\n');
fprintf('SEED SELECTION\n');
fprintf('====================================================\n');

while true

    [seedFile, seedPath] = uigetfile( ...
        {'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
        'Select ONE MNI seed mask');

    % --------------------------------------------------------
    % User pressed Cancel
    % --------------------------------------------------------
    if isequal(seedFile,0)

        if isempty(seed_full)
            error('No seed mask selected.');
        else
            % If seeds were already selected, Cancel simply
            % finishes the seed-selection stage.
            break;
        end
    end


    % --------------------------------------------------------
    % Full path
    % --------------------------------------------------------
    thisSeed = fullfile(seedPath, seedFile);


    % --------------------------------------------------------
    % Generate seed name from filename
    %
    % Examples:
    %   DRN.nii        -> DRN
    %   LC_seed.nii.gz -> LC_seed
    % --------------------------------------------------------

    [~, seedLabel, ext] = fileparts(seedFile);

    % .nii.gz needs fileparts twice
    if strcmpi(ext,'.gz')
        [~, seedLabel, ext2] = fileparts(seedLabel);

        if ~strcmpi(ext2,'.nii')
            warning('Unexpected compressed file extension: %s', seedFile);
        end
    end


    % --------------------------------------------------------
    % Make filename-safe seed label
    % --------------------------------------------------------

    seedLabel = regexprep(seedLabel, '[^A-Za-z0-9_-]', '_');


    % --------------------------------------------------------
    % Avoid duplicate seed names
    % --------------------------------------------------------

    if any(strcmpi(seed_labels, seedLabel))

        warning(['A seed named "%s" has already been selected.\n' ...
                 'This seed will NOT be added again.'], seedLabel);

    else

        seed_full{end+1}   = thisSeed;   %#ok<AGROW>
        seed_labels{end+1} = seedLabel;  %#ok<AGROW>

        fprintf('Added seed: %s\n', seedLabel);
        fprintf('  %s\n', thisSeed);

    end


    % --------------------------------------------------------
    % Ask whether another seed should be selected
    % --------------------------------------------------------

    answer = questdlg( ...
        'Do you want to select another seed?', ...
        'Add another seed?', ...
        'Yes', ...
        'No', ...
        'Yes');


    if isempty(answer) || strcmpi(answer,'No')
        break;
    end

end


if isempty(seed_full)
    error('No seeds selected.');
end


fprintf('\nSelected seeds (%d):\n', numel(seed_full));

for s = 1:numel(seed_full)
    fprintf('  %d: %s\n', s, seed_labels{s});
end



% ============================================================
% SELECT EPI FILES
% ============================================================

[epiFiles, epiPath] = uigetfile( ...
    {'*.nii','NIfTI (*.nii)'}, ...
    'Select reference EPI images (e.g., meanEPI) - MultiSelect', ...
    'MultiSelect','on');

if isequal(epiFiles,0)
    error('No EPI references selected.');
end

if ischar(epiFiles)
    epiFiles = {epiFiles};
end

epi_full = cellfun( ...
    @(f) fullfile(epiPath,f), ...
    epiFiles, ...
    'UniformOutput', false);



% ============================================================
% SELECT DEFORMATION FIELDS
% ============================================================

[yFiles, yPath] = uigetfile( ...
    {'y_*.nii','SPM deformation field (y_*.nii)'}, ...
    'Select deformation fields y_*.nii (T1->MNI) - MultiSelect', ...
    'MultiSelect','on');

if isequal(yFiles,0)
    error('No deformation fields selected.');
end

if ischar(yFiles)
    yFiles = {yFiles};
end

y_full = cellfun( ...
    @(f) fullfile(yPath,f), ...
    yFiles, ...
    'UniformOutput', false);



% ============================================================
% SELECT CSF MAPS
% ============================================================

[c3Files, c3Path] = uigetfile( ...
    {'c3*.nii;*.nii','CSF TPM (c3*.nii)'}, ...
    'Select CSF maps c3*.nii (T1 space) - MultiSelect', ...
    'MultiSelect','on');

if isequal(c3Files,0)
    error('No c3 files selected.');
end

if ischar(c3Files)
    c3Files = {c3Files};
end

c3_full = cellfun( ...
    @(f) fullfile(c3Path,f), ...
    c3Files, ...
    'UniformOutput', false);



% ============================================================
% OUTPUT DIRECTORY
% ============================================================

outDir = uigetdir(epiPath, 'Select output directory');

if isequal(outDir,0)
    error('No output directory selected.');
end

if ~exist(outDir,'dir')
    mkdir(outDir);
end



% ============================================================
% PARAMETERS
% ============================================================

thrCSF = 0.50;      % CSF probability threshold

minVox = 5;         % Warn if seed becomes very small

prefixTmp = 'TMP_'; % Temporary SPM prefix



% ============================================================
% SUBJECT MATCHING
% ============================================================

% Subject ID pattern:
%
% sub-001
% sub-A01
% sub-test123
%
% Only letters/digits after sub- are included.

get_id = @(p) regexp( ...
    p, ...
    'sub-([A-Za-z0-9]+)', ...
    'tokens', ...
    'once');



% ------------------------------------------------------------
% EPI map
% ------------------------------------------------------------

epi_map = containers.Map();

for i = 1:numel(epi_full)

    id = get_id(epi_full{i});

    if isempty(id)
        warning('No sub-<ID> in EPI: %s', epi_full{i});
        continue;
    end

    epi_map(id{1}) = epi_full{i};
end



% ------------------------------------------------------------
% Deformation map
% ------------------------------------------------------------

y_map = containers.Map();

for i = 1:numel(y_full)

    id = get_id(y_full{i});

    if isempty(id)
        warning('No sub-<ID> in y_: %s', y_full{i});
        continue;
    end

    y_map(id{1}) = y_full{i};
end



% ------------------------------------------------------------
% CSF map
% ------------------------------------------------------------

c3_map = containers.Map();

for i = 1:numel(c3_full)

    id = get_id(c3_full{i});

    if isempty(id)
        warning('No sub-<ID> in c3: %s', c3_full{i});
        continue;
    end

    c3_map(id{1}) = c3_full{i};
end



% ------------------------------------------------------------
% Subjects available in ALL three file types
% ------------------------------------------------------------

ids = intersect( ...
    intersect(keys(epi_map), keys(y_map)), ...
    keys(c3_map));


if isempty(ids)

    error(['No matching subjects found. ' ...
           'Ensure filenames include "sub-<ID>" ' ...
           'for EPI, y_, and c3.']);

end


fprintf('\n====================================================\n');
fprintf('MATCHED SUBJECTS: %d\n', numel(ids));
fprintf('====================================================\n');

disp(ids(:));


fprintf('\nNumber of seeds: %d\n', numel(seed_full));
fprintf('Total seed transformations: %d\n', ...
    numel(ids) * numel(seed_full));



% ============================================================
% LOG
% ============================================================

logRows = { ...
    'ID', ...
    'SeedName', ...
    'SeedMNI', ...
    'EPIref', ...
    'yField', ...
    'c3File', ...
    'seed_warped', ...
    'seed_clean', ...
    'vox_warped', ...
    'vox_clean', ...
    'note'};



% ============================================================
% PROCESS SUBJECTS
% ============================================================

for k = 1:numel(ids)

    id     = ids{k};
    epi    = epi_map(id);
    yField = y_map(id);
    c3     = c3_map(id);


    fprintf('\n\n');
    fprintf('====================================================\n');
    fprintf('SUBJECT: sub-%s\n', id);
    fprintf('====================================================\n');

    fprintf('EPI: %s\n', epi);
    fprintf('y_:  %s\n', yField);
    fprintf('c3:  %s\n', c3);


    subjectNote = "";


    % ========================================================
    % SANITY CHECK: deformation field
    % ========================================================

    Vy = spm_vol(yField);

    if numel(Vy) < 3

        warning( ...
            'yField <3 vols for sub-%s. Check deformation field: %s', ...
            id, yField);

        subjectNote = subjectNote + "yField<3vols;";

    end



    % ========================================================
    % STEP 1:
    % RESAMPLE c3(T1) -> EPI GRID
    %
    % Only needs to be done ONCE per subject.
    % All seeds use the same c3-in-EPI map.
    % ========================================================

    fprintf('\nResampling c3 -> EPI grid...\n');


    % Unique temporary prefix
    c3Prefix = sprintf( ...
        '%sc3_sub-%s_', ...
        prefixTmp, id);


    % --------------------------------------------------------
    % Delete old temporary files with same prefix
    % --------------------------------------------------------

    oldTmp = dir(fullfile(outDir, [c3Prefix '*']));

    for ii = 1:numel(oldTmp)

        tmpFile = fullfile( ...
            oldTmp(ii).folder, ...
            oldTmp(ii).name);

        if ~oldTmp(ii).isdir
            delete(tmpFile);
        end

    end


    % --------------------------------------------------------
    % SPM resampling
    % --------------------------------------------------------

    clear matlabbatch


    matlabbatch{1}.spm.util.defs.comp{1}.id.space = {epi};


    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {c3};


    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = ...
        {outDir};


    % Trilinear interpolation for probability map
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 1;


    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;


    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = ...
        [0 0 0];


    matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = ...
        c3Prefix;


    spm_jobman('run', matlabbatch);



    % --------------------------------------------------------
    % Find generated c3 file
    % --------------------------------------------------------

    cand = dir(fullfile(outDir, [c3Prefix '*']));


    if isempty(cand)

        warning( ...
            'No c3(EPI) output found for sub-%s', ...
            id);

        % No CSF map means none of this subject's seeds
        % can be processed correctly.
        for s = 1:numel(seed_full)

            logRows(end+1,:) = { ...
                id, ...
                seed_labels{s}, ...
                seed_full{s}, ...
                epi, ...
                yField, ...
                c3, ...
                '', ...
                '', ...
                NaN, ...
                NaN, ...
                'no_c3_output'}; %#ok<AGROW>

        end

        continue;

    end


    % If several candidates exist, use newest one
    [~,ix] = max([cand.datenum]);


    c3_epi_tmp = fullfile( ...
        cand(ix).folder, ...
        cand(ix).name);


    % --------------------------------------------------------
    % Final c3 filename
    % --------------------------------------------------------

    c3_epi = fullfile( ...
        outDir, ...
        sprintf('sub-%s_c3_inEPI.nii', id));


    if exist(c3_epi,'file')
        delete(c3_epi);
    end


    movefile(c3_epi_tmp, c3_epi);


    % --------------------------------------------------------
    % Read CSF probability map once
    % --------------------------------------------------------

    Vc = spm_vol(c3_epi);

    C = spm_read_vols(Vc);


    % Binary CSF mask
    csfMask = (C > thrCSF);


    fprintf('c3 resampled successfully:\n');
    fprintf('  %s\n', c3_epi);



    % ========================================================
    % STEP 2:
    % PROCESS EVERY SEED FOR THIS SUBJECT
    % ========================================================

    for s = 1:numel(seed_full)

        mask_mni  = seed_full{s};
        seedLabel = seed_labels{s};


        fprintf('\n----------------------------------------------------\n');
        fprintf('Seed %d/%d: %s\n', ...
            s, numel(seed_full), seedLabel);
        fprintf('----------------------------------------------------\n');

        fprintf('MNI seed:\n');
        fprintf('  %s\n', mask_mni);


        note = subjectNote;



        % ====================================================
        % WARP MNI SEED -> SUBJECT EPI
        % ====================================================

        clear matlabbatch


        % ----------------------------------------------------
        % Inverse deformation:
        % MNI -> native subject
        % ----------------------------------------------------

        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = ...
            {yField};


        matlabbatch{1}.spm.util.defs.comp{1}.inv.space = ...
            {epi};


        % ----------------------------------------------------
        % Explicit target voxel grid = EPI
        % ----------------------------------------------------

        matlabbatch{1}.spm.util.defs.comp{2}.id.space = ...
            {epi};


        % ----------------------------------------------------
        % Seed input
        % ----------------------------------------------------

        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = ...
            {mask_mni};


        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = ...
            {outDir};


        % Nearest neighbour for binary seed
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;


        % Do not mask outside deformation field
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 0;


        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = ...
            [0 0 0];


        % ----------------------------------------------------
        % Unique temporary prefix:
        %
        % TMP_seed_sub-001_DRN_
        % TMP_seed_sub-001_LC_
        % ----------------------------------------------------

        seedPrefix = sprintf( ...
            '%sseed_sub-%s_%s_', ...
            prefixTmp, ...
            id, ...
            seedLabel);


        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = ...
            seedPrefix;



        % ----------------------------------------------------
        % Delete old temporary files for same subject/seed
        % ----------------------------------------------------

        oldTmp = dir(fullfile(outDir, [seedPrefix '*']));

        for ii = 1:numel(oldTmp)

            tmpFile = fullfile( ...
                oldTmp(ii).folder, ...
                oldTmp(ii).name);

            if ~oldTmp(ii).isdir
                delete(tmpFile);
            end

        end



        % ----------------------------------------------------
        % Run SPM
        % ----------------------------------------------------

        spm_jobman('run', matlabbatch);



        % ====================================================
        % FIND WARPED SEED
        % ====================================================

        cand = dir(fullfile( ...
            outDir, ...
            [seedPrefix '*']));


        if isempty(cand)

            warning( ...
                'No warped output for sub-%s / seed %s', ...
                id, seedLabel);


            logRows(end+1,:) = { ...
                id, ...
                seedLabel, ...
                mask_mni, ...
                epi, ...
                yField, ...
                c3, ...
                '', ...
                '', ...
                NaN, ...
                NaN, ...
                'no_warp_output'}; %#ok<AGROW>


            continue;

        end


        % Use newest matching file
        [~,ix] = max([cand.datenum]);


        seed_warp_tmp = fullfile( ...
            cand(ix).folder, ...
            cand(ix).name);



        % ====================================================
        % FINAL WARPED-SEED NAME
        %
        % Example:
        % sub-001_DRN_inEPI.nii
        % ====================================================

        seed_warp = fullfile( ...
            outDir, ...
            sprintf( ...
                'sub-%s_%s_inEPI.nii', ...
                id, ...
                seedLabel));


        if exist(seed_warp,'file')
            delete(seed_warp);
        end


        movefile(seed_warp_tmp, seed_warp);



        % ====================================================
        % BINARIZE WARPED SEED
        % ====================================================

        Vsw = spm_vol(seed_warp);

        Sw = spm_read_vols(Vsw);


        Swb = double(Sw > 0.5);


        % Defensive overwrite as binary image
        spm_write_vol(Vsw, Swb);


        vox_warped = nnz(Swb);



        % ====================================================
        % CHECK DIMENSIONS
        % ====================================================

        if ~isequal(size(Swb), size(csfMask))

            warning( ...
                ['Dimension mismatch for sub-%s / seed %s.\n' ...
                 'Seed size: [%s]\n' ...
                 'c3 size:   [%s]'], ...
                id, ...
                seedLabel, ...
                num2str(size(Swb)), ...
                num2str(size(csfMask)));


            logRows(end+1,:) = { ...
                id, ...
                seedLabel, ...
                mask_mni, ...
                epi, ...
                yField, ...
                c3, ...
                seed_warp, ...
                '', ...
                vox_warped, ...
                NaN, ...
                'dimension_mismatch'}; %#ok<AGROW>


            continue;

        end



        % ====================================================
        % EXCLUDE CSF
        % ====================================================

        seedClean = logical(Swb) & ~csfMask;


        vox_clean = nnz(seedClean);



        % ====================================================
        % WARNINGS
        % ====================================================

        if vox_clean < minVox

            warning( ...
                'sub-%s / %s: Seed after CSF exclusion is very small (%d vox).', ...
                id, ...
                seedLabel, ...
                vox_clean);


            note = note + ...
                sprintf("seed_clean_small(%d);", vox_clean);

        end


        if vox_clean == 0

            warning( ...
                'sub-%s / %s: Seed after CSF exclusion is EMPTY.', ...
                id, ...
                seedLabel);


            note = note + ...
                "seed_clean_EMPTY;";

        end



        % ====================================================
        % WRITE CLEAN SEED
        %
        % Example:
        % sub-001_DRN_inEPI_CSFexcluded.nii
        % ====================================================

        seed_clean_file = fullfile( ...
            outDir, ...
            sprintf( ...
                'sub-%s_%s_inEPI_CSFexcluded.nii', ...
                id, ...
                seedLabel));


        Vout = Vsw;

        Vout.fname = seed_clean_file;


        spm_write_vol( ...
            Vout, ...
            double(seedClean));



        % ====================================================
        % CONSOLE OUTPUT
        % ====================================================

        fprintf('\nSUCCESS:\n');

        fprintf( ...
            'Warped seed:   %s\n', ...
            seed_warp);

        fprintf( ...
            '  voxels = %d\n', ...
            vox_warped);


        fprintf( ...
            'CSF excluded:  %s\n', ...
            seed_clean_file);

        fprintf( ...
            '  voxels = %d\n', ...
            vox_clean);



        % ====================================================
        % LOG ROW
        % ====================================================

        logRows(end+1,:) = { ...
            id, ...
            seedLabel, ...
            mask_mni, ...
            epi, ...
            yField, ...
            c3, ...
            seed_warp, ...
            seed_clean_file, ...
            vox_warped, ...
            vox_clean, ...
            char(note)}; %#ok<AGROW>


        % ----------------------------------------------------
        % Optional visual QC
        % ----------------------------------------------------
        %
        % Uncomment if desired:
        %
        % spm_check_registration( ...
        %     char(epi, seed_clean_file));
        %

    end

end



% ============================================================
% SAVE CSV LOG
% ============================================================

csvFile = fullfile( ...
    outDir, ...
    'Seed_warp_and_CSFexclude_log.csv');


try

    writecell(logRows, csvFile);

catch

    % --------------------------------------------------------
    % Fallback for older MATLAB versions
    % --------------------------------------------------------

    fid = fopen(csvFile,'w');


    fprintf(fid, ...
        ['ID,SeedName,SeedMNI,EPIref,yField,c3File,' ...
         'seed_warped,seed_clean,vox_warped,vox_clean,note\n']);


    for r = 2:size(logRows,1)

        fprintf(fid, ...
            '%s,"%s","%s","%s","%s","%s","%s","%s",%g,%g,"%s"\n', ...
            logRows{r,1}, ...
            logRows{r,2}, ...
            logRows{r,3}, ...
            logRows{r,4}, ...
            logRows{r,5}, ...
            logRows{r,6}, ...
            logRows{r,7}, ...
            logRows{r,8}, ...
            logRows{r,9}, ...
            logRows{r,10}, ...
            logRows{r,11});

    end


    fclose(fid);

end



% ============================================================
% FINISHED
% ============================================================

fprintf('\n\n');
fprintf('====================================================\n');
fprintf('DONE\n');
fprintf('====================================================\n');

fprintf('\nSubjects processed: %d\n', numel(ids));
fprintf('Seeds selected:     %d\n', numel(seed_full));

fprintf('\nSeeds:\n');

for s = 1:numel(seed_labels)
    fprintf('  - %s\n', seed_labels{s});
end

fprintf('\nLog saved:\n');
fprintf('  %s\n\n', csvFile);

end