%% LC-Volumen, Peak/Mean und Pons-Kontrast Berechnung (Nummern-Matching)
clc; clear;

%% --- 1. LC-Masken auswählen ---
[files_masks, path_masks] = uigetfile({'*.nii;*.nii.gz','NIfTI files (*.nii, *.nii.gz)'}, ...
   'Select LC mask files', 'MultiSelect', 'on');
if isequal(files_masks,0), disp('Keine Maskendatei ausgewählt. Abbruch.'); return; end
if ~iscell(files_masks), files_masks = {files_masks}; end
n = length(files_masks);

%% --- 2. Option: Originalbilder auswählen ---
answer = questdlg('Sollen Peak/Mean aus Originalbildern berechnet werden?', ...
   'Intensity source','Ja','Nein','Nein');
use_orig_images = strcmp(answer,'Ja');
if use_orig_images
   [files_imgs, path_imgs] = uigetfile({'*.nii;*.nii.gz','NIfTI files (*.nii, *.nii.gz)'}, ...
       'Select original images (any number)', 'MultiSelect', 'on');
   if isequal(files_imgs,0)
       disp('Keine Originalbilder ausgewählt. Peak/Mean werden übersprungen.');
       use_orig_images = false;
   else
       if ~iscell(files_imgs), files_imgs = {files_imgs}; end
   end
end

%% --- 3. Option: Pons-Masken auswählen ---
answer_pons = questdlg('Sollen Pons-Masken für Kontrast verwendet werden?', ...
   'Pons reference','Ja','Nein','Nein');
use_pons_masks = strcmp(answer_pons,'Ja');
if use_pons_masks
   [files_pons, path_pons] = uigetfile({'*.nii;*.nii.gz','NIfTI files (*.nii, *.nii.gz)'}, ...
       'Select Pons mask(es)', 'MultiSelect', 'on');
   if isequal(files_pons,0)
       warning('Keine Pons-Masken ausgewählt. Kontrast wird übersprungen.');
       use_pons_masks = false;
   else
       if ~iscell(files_pons), files_pons = {files_pons}; end
   end
end

%% --- 4. Ergebnisse vorbereiten ---
Subj = cell(n,1);
Volume_mm3 = zeros(n,1);
VoxelCount = zeros(n,1);
PeakIntensity = nan(n,1);
MeanIntensity = nan(n,1);
LC_Pons_Contrast = nan(n,1);
MaskPath = cell(n,1);
ImgPath = cell(n,1);

%% --- 5. Schleife über Masken ---
for i = 1:n
   mask_file = fullfile(path_masks, files_masks{i});
   MaskPath{i} = mask_file;
   
   % Maske einlesen
   mask = niftiread(mask_file);
   info_mask = niftiinfo(mask_file);
   bw = mask>0.5;
   voxel_size_mm3 = prod(info_mask.PixelDimensions);
   
   % Volumen
   num_vox = nnz(bw);
   Volume_mm3(i) = num_vox*voxel_size_mm3;
   VoxelCount(i) = num_vox;
   
   % Subjekt-Nummer extrahieren
   numMask = regexp(files_masks{i}, '(?i)sub-(\d{2,5})', 'tokens', 'once');
   if isempty(numMask)
       [~, nme, ~] = fileparts(files_masks{i});
       Subj{i} = nme;
       warning('Keine Subjekt-Nummer in %s gefunden.', files_masks{i});
   else
       Subj{i} = ['sub-' numMask{1}];
   end
   
   if num_vox==0
       warning('Maske %s enthält keine Voxels.', files_masks{i});
   end
   
   %% --- Peak/Mean aus Originalbild ---
   if use_orig_images
       idx_img = [];
       numMaskStr = numMask{1};
       for j = 1:length(files_imgs)
           numImg = regexp(files_imgs{j}, '(?i)sub-(\d{2,5})', 'tokens', 'once');
           if ~isempty(numImg) && strcmp(numMaskStr, numImg{1})
               idx_img = j;
               break;
           end
       end
       if isempty(idx_img)
           warning('Kein Originalbild für Subjekt %s gefunden.', Subj{i});
       else
           img_file = fullfile(path_imgs, files_imgs{idx_img});
           ImgPath{i} = img_file;
           try
               img = niftiread(img_file);
               if isequal(size(img), size(bw))
                   vals = double(img(bw));
                   if ~isempty(vals)
                       PeakIntensity(i) = max(vals);
                       MeanIntensity(i) = mean(vals);
                   end
               else
                   warning('Größenunterschied zwischen Maske und Originalbild für %s.', Subj{i});
               end
           catch ME
               warning('Fehler beim Einlesen Originalbild %s: %s', files_imgs{idx_img}, ME.message);
           end
       end
   end
   
   %% --- LC-Pons-Kontrast ---
   if use_pons_masks && ~isempty(ImgPath{i})
       % Falls nur eine Pons-Maske → für alle Subjekte
       if length(files_pons)==1
           pons_mask = niftiread(fullfile(path_pons, files_pons{1}));
       else
           % Matching nach Nummer
           idx_pons = [];
           for j = 1:length(files_pons)
               numPons = regexp(files_pons{j}, '(?i)sub-(\d{2,5})', 'tokens', 'once');
               if ~isempty(numPons) && strcmp(numMaskStr, numPons{1})
                   idx_pons = j;
                   break;
               end
           end
           if isempty(idx_pons)
               warning('Keine Pons-Maske für Subjekt %s gefunden.', Subj{i});
               continue;
           end
           pons_mask = niftiread(fullfile(path_pons, files_pons{idx_pons}));
       end
       img = niftiread(ImgPath{i});
       if isequal(size(img), size(pons_mask)) && ~isnan(MeanIntensity(i))
           vals_pons = double(img(pons_mask>0.5));
           if ~isempty(vals_pons)
               LC_Pons_Contrast(i) = MeanIntensity(i)/mean(vals_pons);
           end
       else
           warning('Größenunterschied zwischen Pons-Maske und Bild für %s.', Subj{i});
       end
   end
end

%% --- 6. Tabelle erstellen ---
T = table(Subj, VoxelCount, Volume_mm3, PeakIntensity, MeanIntensity, LC_Pons_Contrast, MaskPath, ImgPath, ...
   'VariableNames', {'Subject','VoxelCount','LC_Volume_mm3','LC_Peak','LC_Mean','LC_Pons_Contrast','MaskPath','ImagePath'});
disp(T);

%% --- 7. CSV speichern ---
save_csv = questdlg('CSV-Datei speichern?', 'Export', 'Ja', 'Nein', 'Ja');
if strcmp(save_csv,'Ja')
   csv_file = fullfile(path_masks, 'LC_Results.csv');
   writetable(T, csv_file);
   fprintf('Tabelle gespeichert: %s\n', csv_file);
end

%% --- 8. Kurze Summary ---
fprintf('\nSummary (n = %d):\n', n);
validVol = ~isnan(Volume_mm3) & Volume_mm3>0;
if any(validVol)
   fprintf('Volume (mm^3): mean=%.2f, sd=%.2f, median=%.2f\n', ...
       mean(Volume_mm3(validVol)), std(Volume_mm3(validVol)), median(Volume_mm3(validVol)));
end
validMean = ~isnan(MeanIntensity);
if any(validMean)
   fprintf('LC Mean Intensity: mean=%.3f, sd=%.3f\n', mean(MeanIntensity(validMean)), std(MeanIntensity(validMean)));
end
validContrast = ~isnan(LC_Pons_Contrast);
if any(validContrast)
   fprintf('LC/Pons Contrast: mean=%.3f, sd=%.3f\n', mean(LC_Pons_Contrast(validContrast)), std(LC_Pons_Contrast(validContrast)));
end
fprintf('Fertig.\n');
 
