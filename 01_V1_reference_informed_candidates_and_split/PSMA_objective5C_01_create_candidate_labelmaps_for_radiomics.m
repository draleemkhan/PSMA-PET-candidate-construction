%% PSMA_objective5C_01_create_candidate_labelmaps_for_radiomics.m
% Objective 5C - Step 01
% Create per-case candidate labelmaps for radiomics extraction.
%
% This script regenerates Objective 4-style candidate objects:
%   TRUE candidates  = tumor-mask connected components >= minVol_mL
%   FALSE candidates = PET hot regions SUV >= SUVthr outside 1-voxel dilated tumor mask,
%                      connected components >= minVol_mL
%
% Instead of saving one NIfTI mask per candidate, it saves ONE candidate labelmap
% per case. Each candidate has a unique integer label value within that case.
% This is much more efficient for PyRadiomics.
%
% Outputs:
%   1) candidate_labelmaps\<case_id>_Objective5C_candidate_labelmap.nii
%   2) PSMA_objective5C_candidate_labelmap_index.xlsx/csv
%   3) PyRadiomics batch CSV files for PET and CT
%   4) processing log and error table
%
% Requirements:
%   - MATLAB Image Processing Toolbox
%   - corrected D-drive manifest:
%       D:\PSMA-PET-CT PROJECT\PSMA_AI_Project\01_manifest\PSMA_final_manifest_D_DRIVE.xlsx
%
% Author: PSMA PET/CT project workflow
% -------------------------------------------------------------------------

clear; clc;

%% ---------------- USER SETTINGS ----------------
projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));

manifestPath = fullfile(projectRoot, "01_manifest", "PSMA_final_manifest_D_DRIVE.xlsx");

% Existing Objective 4 candidate table is optional but useful for split labels.
objective4CandidateTable = fullfile(projectRoot, "06_results", ...
    "Objective4_TRUE_vs_FALSE_candidate_dataset", ...
    "PSMA_objective4_candidate_features_TRUE_vs_FALSE.xlsx");

outRoot = fullfile(projectRoot, "06_results", "Objective5C_Radiomics");
labelmapDir = fullfile(outRoot, "candidate_masks");
featureDir  = fullfile(outRoot, "radiomics_features");
logDir      = fullfile(outRoot, "logs");

% Candidate-generation settings used for Objective 4 radiomics candidate table.
SUVthr = 3.0;
minVol_mL = 1.0;
dilationRadiusVox = 1;   % one-voxel tumor-mask dilation
connectivity = 26;

% Run mode. Use 'pilot' first; then change to 'full'.
runMode = 'full';       % 'pilot' or 'full'
pilotNcases = 10;

% Overwrite existing per-case labelmaps?
overwriteExisting = true;

%% ---------------- CREATE FOLDERS ----------------
if ~isfolder(outRoot), mkdir(outRoot); end
if ~isfolder(labelmapDir), mkdir(labelmapDir); end
if ~isfolder(featureDir), mkdir(featureDir); end
if ~isfolder(logDir), mkdir(logDir); end

logFile = fullfile(logDir, "Objective5C_candidate_labelmap_creation_log.txt");
if isfile(logFile), delete(logFile); end
diary(logFile);

fprintf("\n============================================================\n");
fprintf("Objective 5C Step 01: Candidate labelmap creation for radiomics\n");
fprintf("Project root: %s\n", projectRoot);
fprintf("Manifest:     %s\n", manifestPath);
fprintf("Output root:  %s\n", outRoot);
fprintf("SUV threshold: %.2f | min volume: %.2f mL | connectivity: %d\n", SUVthr, minVol_mL, connectivity);
fprintf("Run mode: %s\n", runMode);
fprintf("============================================================\n\n");

%% ---------------- LOAD MANIFEST ----------------
assert(isfile(manifestPath), "Manifest not found: %s", manifestPath);
M = readtable(manifestPath, 'VariableNamingRule', 'preserve');
M = normalizeVarNames(M);

requiredManifestVars = ["case_id", "CT_path", "PET_path", "Label_path"];
assert(all(ismember(requiredManifestVars, string(M.Properties.VariableNames))), ...
    "Manifest must contain case_id, CT_path, PET_path, and Label_path columns.");

% Ensure path variables are string.
M.case_id = string(M.case_id);
M.CT_path = string(M.CT_path);
M.PET_path = string(M.PET_path);
M.Label_path = string(M.Label_path);

% Standardize optional metadata variable names.
if ~ismember("SubjectID", string(M.Properties.VariableNames)) && ismember("subjectid", lower(string(M.Properties.VariableNames)))
    % already handled poorly by preserve; leave later fallback
end

%% ---------------- OPTIONAL: LOAD OBJECTIVE 4 SPLIT INFO ----------------
splitMap = table();
if isfile(objective4CandidateTable)
    fprintf("Loading Objective 4 candidate table for split labels...\n");
    C4 = readtable(objective4CandidateTable, 'VariableNamingRule', 'preserve');
    C4 = normalizeVarNames(C4);
    if all(ismember(["case_id", "split_set"], string(C4.Properties.VariableNames)))
        C4.case_id = string(C4.case_id);
        C4.split_set = string(C4.split_set);
        splitMap = unique(C4(:, {"case_id", "split_set"}));
        fprintf("Loaded split labels for %d unique cases.\n", height(splitMap));
    else
        fprintf("Objective 4 table found, but split_set column not found. Continuing without split labels.\n");
    end
else
    fprintf("Objective 4 candidate table not found. Continuing without split labels.\n");
end

%% ---------------- CASE SELECTION ----------------
% Use cases with valid files only.
validCase = isfile(M.CT_path) & isfile(M.PET_path) & isfile(M.Label_path);
fprintf("Valid CT/PET/Label cases in manifest: %d / %d\n", sum(validCase), height(M));
M = M(validCase, :);

switch lower(runMode)
    case 'pilot'
        M = M(1:min(pilotNcases, height(M)), :);
    case 'full'
        % keep all valid cases
    otherwise
        error("runMode must be 'pilot' or 'full'.");
end

fprintf("Cases to process in this run: %d\n\n", height(M));

%% ---------------- PROCESS CASES ----------------
indexRows = cell(0, 26);
indexVarNames = { ...
    'case_id', 'SubjectID', 'StudyDate', 'candidate_labelmap_path', ...
    'candidate_label_value', 'candidate_global_id', 'candidate_source', ...
    'candidate_label_binary', 'candidate_label_name', 'candidate_number_in_case', ...
    'candidate_voxels', 'candidate_volume_mL', ...
    'PET_path', 'CT_path', 'Label_path', ...
    'SUVthr', 'minVol_mL', 'voxelVol_mL', ...
    'centroid_x_voxel', 'centroid_y_voxel', 'centroid_z_voxel', ...
    'scanner', 'radionuclide', 'contrast', 'split_set', 'notes'};

errorRows = cell(0, 4);
errorVarNames = {'case_id', 'stage', 'message', 'time'};

globalCounter = 0;

for i = 1:height(M)
    caseID = string(M.case_id(i));
    fprintf("[%d/%d] %s\n", i, height(M), caseID);

    try
        PETpath = string(M.PET_path(i));
        CTpath  = string(M.CT_path(i));
        LABpath = string(M.Label_path(i));

        info = niftiinfo(PETpath);
        PET = double(niftiread(PETpath));
        LAB = niftiread(LABpath) > 0;

        if ~isequal(size(PET), size(LAB))
            error("PET and label size mismatch: PET %s vs Label %s", mat2str(size(PET)), mat2str(size(LAB)));
        end

        voxelSize = getVoxelSize_mm(info);
        voxelVol_mL = prod(voxelSize) / 1000;
        minVox = max(1, ceil(minVol_mL / voxelVol_mL));

        candidateLabelMap = zeros(size(PET), 'uint16');
        localCounter = 0;

        % ---------- TRUE tumor candidates ----------
        CCtrue = bwconncomp(LAB, connectivity);
        for k = 1:CCtrue.NumObjects
            voxIdx = CCtrue.PixelIdxList{k};
            nVox = numel(voxIdx);
            vol_mL = nVox * voxelVol_mL;
            if vol_mL < minVol_mL
                continue;
            end

            localCounter = localCounter + 1;
            globalCounter = globalCounter + 1;
            candidateLabelMap(voxIdx) = uint16(localCounter);
            [cx, cy, cz] = centroidFromIdx(voxIdx, size(PET));

            indexRows(end+1, :) = makeIndexRow(M, i, labelmapPathForCase(labelmapDir, caseID), ...
                localCounter, globalCounter, "true_tumor_mask_component", 1, "true_tumor", ...
                localCounter, nVox, vol_mL, PETpath, CTpath, LABpath, SUVthr, minVol_mL, voxelVol_mL, ...
                cx, cy, cz, getSplit(splitMap, caseID), ""); %#ok<SAGROW>
        end

        % ---------- FALSE hot-region candidates ----------
        hotMask = PET >= SUVthr;
        if dilationRadiusVox > 0
            se = true(1 + 2*dilationRadiusVox, 1 + 2*dilationRadiusVox, 1 + 2*dilationRadiusVox);
            dilatedTumor = imdilate(LAB, se);
        else
            dilatedTumor = LAB;
        end
        falseMask = hotMask & ~dilatedTumor;
        falseMask = bwareaopen(falseMask, minVox, connectivity);

        CCfalse = bwconncomp(falseMask, connectivity);
        for k = 1:CCfalse.NumObjects
            voxIdx = CCfalse.PixelIdxList{k};
            nVox = numel(voxIdx);
            vol_mL = nVox * voxelVol_mL;
            if vol_mL < minVol_mL
                continue;
            end

            localCounter = localCounter + 1;
            globalCounter = globalCounter + 1;

            if localCounter > intmax('uint16')
                error("Too many candidates in one case for uint16 labelmap.");
            end

            candidateLabelMap(voxIdx) = uint16(localCounter);
            [cx, cy, cz] = centroidFromIdx(voxIdx, size(PET));

            indexRows(end+1, :) = makeIndexRow(M, i, labelmapPathForCase(labelmapDir, caseID), ...
                localCounter, globalCounter, "false_hot_region_outside_dilated_tumor", 0, "false_hot_region", ...
                localCounter, nVox, vol_mL, PETpath, CTpath, LABpath, SUVthr, minVol_mL, voxelVol_mL, ...
                cx, cy, cz, getSplit(splitMap, caseID), ""); %#ok<SAGROW>
        end

        % Save candidate labelmap.
        outLabelmap = labelmapPathForCase(labelmapDir, caseID);
        if overwriteExisting || ~isfile(outLabelmap)
            outInfo = info;
            outInfo.Datatype = 'uint16';
            outInfo.BitsPerPixel = 16;
            niftiwrite(candidateLabelMap, outLabelmap, outInfo, 'Compressed', false);
        end

        fprintf("   saved labelmap with %d candidates: %s\n", localCounter, outLabelmap);

    catch ME
        fprintf(2, "   ERROR: %s\n", ME.message);
        errorRows(end+1, :) = {caseID, "case_processing", string(ME.message), string(datetime('now'))}; %#ok<SAGROW>
        continue;
    end
end

%% ---------------- SAVE INDEX TABLES ----------------
if isempty(indexRows)
    Tindex = cell2table(cell(0, numel(indexVarNames)), 'VariableNames', indexVarNames);
else
    Tindex = cell2table(indexRows, 'VariableNames', indexVarNames);
end

indexXlsx = fullfile(outRoot, "PSMA_objective5C_candidate_labelmap_index.xlsx");
indexCsv  = fullfile(outRoot, "PSMA_objective5C_candidate_labelmap_index.csv");
writetable(Tindex, indexXlsx);
writetable(Tindex, indexCsv);

fprintf("\nSaved candidate index:\n  %s\n  %s\n", indexXlsx, indexCsv);

if isempty(errorRows)
    Terror = cell2table(cell(0, numel(errorVarNames)), 'VariableNames', errorVarNames);
else
    Terror = cell2table(errorRows, 'VariableNames', errorVarNames);
end
errorXlsx = fullfile(logDir, "PSMA_objective5C_candidate_labelmap_errors.xlsx");
writetable(Terror, errorXlsx);
fprintf("Saved error table:\n  %s\n", errorXlsx);

%% ---------------- CREATE PYRADIOMICS BATCH FILES ----------------
% PyRadiomics supports batch CSV with columns Image, Mask, Label.
% We create separate PET and CT batch files because radiomic intensity features
% should be extracted separately from PET and CT images.

if height(Tindex) > 0
    PETbatch = table();
    PETbatch.Image = string(Tindex.PET_path);
    PETbatch.Mask  = string(Tindex.candidate_labelmap_path);
    PETbatch.Label = double(Tindex.candidate_label_value);
    PETbatch.case_id = string(Tindex.case_id);
    PETbatch.candidate_global_id = double(Tindex.candidate_global_id);
    PETbatch.candidate_label_binary = double(Tindex.candidate_label_binary);
    PETbatch.candidate_source = string(Tindex.candidate_source);
    PETbatch.split_set = string(Tindex.split_set);

    CTbatch = table();
    CTbatch.Image = string(Tindex.CT_path);
    CTbatch.Mask  = string(Tindex.candidate_labelmap_path);
    CTbatch.Label = double(Tindex.candidate_label_value);
    CTbatch.case_id = string(Tindex.case_id);
    CTbatch.candidate_global_id = double(Tindex.candidate_global_id);
    CTbatch.candidate_label_binary = double(Tindex.candidate_label_binary);
    CTbatch.candidate_source = string(Tindex.candidate_source);
    CTbatch.split_set = string(Tindex.split_set);
else
    PETbatch = table();
    CTbatch = table();
end

petBatchCsv = fullfile(featureDir, "pyradiomics_batch_PET_candidates.csv");
ctBatchCsv  = fullfile(featureDir, "pyradiomics_batch_CT_candidates.csv");
writetable(PETbatch, petBatchCsv);
writetable(CTbatch, ctBatchCsv);

fprintf("Saved PyRadiomics batch CSV files:\n  %s\n  %s\n", petBatchCsv, ctBatchCsv);

%% ---------------- SUMMARY ----------------
fprintf("\n============================================================\n");
fprintf("DONE Objective 5C Step 01.\n");
fprintf("Cases processed:        %d\n", height(M));
fprintf("Candidate rows created: %d\n", height(Tindex));
fprintf("Errors:                 %d\n", height(Terror));
fprintf("Output folder:          %s\n", outRoot);
fprintf("\nNext step: run PyRadiomics using the batch CSV files.\n");
fprintf("============================================================\n");

diary off;

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function T = normalizeVarNames(T)
    % Converts common path column variants to standard names while preserving others.
    vn = string(T.Properties.VariableNames);

    replacements = {
        'CaseID', 'case_id'
        'caseID', 'case_id'
        'case', 'case_id'
        'CTpath', 'CT_path'
        'ct_path', 'CT_path'
        'PETpath', 'PET_path'
        'pet_path', 'PET_path'
        'SUV_path', 'PET_path'
        'Labelpath', 'Label_path'
        'label_path', 'Label_path'
        'Mask_path', 'Label_path'
        'split', 'split_set'
        'Split', 'split_set'
    };

    for r = 1:size(replacements,1)
        old = string(replacements{r,1});
        new = string(replacements{r,2});
        idx = find(vn == old, 1);
        if ~isempty(idx)
            vn(idx) = new;
        end
    end

    % Prevent duplicate names after replacement.
    T.Properties.VariableNames = matlab.lang.makeUniqueStrings(cellstr(vn));
end

function voxelSize = getVoxelSize_mm(info)
    if isfield(info, 'PixelDimensions') && numel(info.PixelDimensions) >= 3
        voxelSize = double(info.PixelDimensions(1:3));
    else
        voxelSize = [1 1 1];
    end
end

function [cx, cy, cz] = centroidFromIdx(voxIdx, imgSize)
    [x, y, z] = ind2sub(imgSize, voxIdx);
    cx = mean(x);
    cy = mean(y);
    cz = mean(z);
end

function outPath = labelmapPathForCase(labelmapDir, caseID)
    safeCase = regexprep(char(caseID), '[^a-zA-Z0-9_\-]', '_');
    outPath = fullfile(labelmapDir, sprintf('%s_Objective5C_candidate_labelmap.nii', safeCase));
end

function split = getSplit(splitMap, caseID)
    split = "unknown";
    if isempty(splitMap) || ~all(ismember(["case_id", "split_set"], string(splitMap.Properties.VariableNames)))
        return;
    end
    idx = find(string(splitMap.case_id) == string(caseID), 1);
    if ~isempty(idx)
        split = string(splitMap.split_set(idx));
    end
end

function row = makeIndexRow(M, i, labelmapPath, labelValue, globalID, source, labelBinary, labelName, numInCase, nVox, vol_mL, PETpath, CTpath, LABpath, SUVthr, minVol_mL, voxelVol_mL, cx, cy, cz, splitSet, notes)
    caseID = string(M.case_id(i));

    SubjectID = getOptionalValue(M, i, ["SubjectID", "subjectID", "subjectid"], "");
    StudyDate = getOptionalValue(M, i, ["StudyDate", "studyDate", "studydate"], "");
    scanner = getOptionalValue(M, i, ["scanner", "manufacturer_model_name", "ManufacturerModelName"], "unknown");
    radionuclide = getOptionalValue(M, i, ["radionuclide", "pet_radionuclide", "PET_radionuclide"], "unknown");
    contrast = getOptionalValue(M, i, ["contrast", "ct_contrast_agent", "CT_contrast_agent"], "unknown");

    row = { ...
        caseID, string(SubjectID), string(StudyDate), string(labelmapPath), ...
        double(labelValue), double(globalID), string(source), ...
        double(labelBinary), string(labelName), double(numInCase), ...
        double(nVox), double(vol_mL), ...
        string(PETpath), string(CTpath), string(LABpath), ...
        double(SUVthr), double(minVol_mL), double(voxelVol_mL), ...
        double(cx), double(cy), double(cz), ...
        string(scanner), string(radionuclide), string(contrast), string(splitSet), string(notes)};
end

function val = getOptionalValue(T, i, possibleNames, defaultVal)
    vn = string(T.Properties.VariableNames);
    val = defaultVal;
    for k = 1:numel(possibleNames)
        idx = find(vn == possibleNames(k), 1);
        if ~isempty(idx)
            try
                raw = T{i, idx};
                if iscell(raw)
                    raw = raw{1};
                end
                val = string(raw);
                return;
            catch
                val = defaultVal;
                return;
            end
        end
    end
end
