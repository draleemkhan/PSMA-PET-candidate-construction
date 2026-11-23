clear; clc; close all;

%% ============================================================
%  OBJECTIVE 4A
%  CREATE TRUE-TUMOR VS FALSE-POSITIVE PET HOT-REGION CANDIDATE TABLE
%
%  TRUE candidates:
%    tumor mask connected components >= 1 mL
%
%  FALSE candidates:
%    PET hot regions with SUV >= threshold
%    outside dilated tumor mask
%    volume >= 1 mL
%
%  Confirmed orientation:
%    x_low  = right
%    x_mid  = midline
%    x_high = left
%
%    y_low  = posterior
%    y_mid  = central
%    y_high = anterior
%
%    z_low_index_third    = superior
%    z_middle_index_third = middle
%    z_high_index_third   = inferior
%  ============================================================

projectDir = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));

manifestFile = fullfile(projectDir, ...
    "01_manifest", ...
    "PSMA_final_manifest.xlsx");

resultDir = fullfile(projectDir, ...
    "06_results", ...
    "Objective4_TRUE_vs_FALSE_candidate_dataset");

if ~exist(resultDir, "dir")
    mkdir(resultDir);
end

outCandidateXlsx = fullfile(resultDir, ...
    "PSMA_objective4_candidate_features_TRUE_vs_FALSE.xlsx");

outCandidateCsv = fullfile(resultDir, ...
    "PSMA_objective4_candidate_features_TRUE_vs_FALSE.csv");

outSummaryXlsx = fullfile(resultDir, ...
    "PSMA_objective4_candidate_generation_summary.xlsx");

%% ============================================================
%  SETTINGS
%  ============================================================

SUV_THRESHOLD = 3.0;              % threshold for PET hot regions
MIN_CANDIDATE_VOLUME_ML = 1.0;    % minimum candidate size
DILATE_TUMOR_MASK = true;         % remove immediate tumor boundary from false candidates

fprintf("Objective 4A: TRUE vs FALSE candidate table generation\n");
fprintf("SUV threshold for false hot candidates: %.2f\n", SUV_THRESHOLD);
fprintf("Minimum candidate volume: %.2f mL\n", MIN_CANDIDATE_VOLUME_ML);

%% ============================================================
%  READ MANIFEST
%  ============================================================

M = readtable(manifestFile);

fprintf("Total manifest cases: %d\n", height(M));

%% ============================================================
%  STORAGE
%  ============================================================

allRows = {};
errorRows = {};

caseSummaryRows = {};

globalCandidateID = 0;

%% ============================================================
%  MAIN LOOP
%  ============================================================

for i = 1:height(M)

    try

        case_id = getTableString(M.case_id, i);

        fprintf("\nProcessing case %d/%d: %s\n", i, height(M), case_id);

        ctPath    = getTableString(M.CT_path, i);
        petPath   = getTableString(M.PET_path, i);
        labelPath = getTableString(M.Label_path, i);

        SubjectID = getTableString(M.SubjectID, i);
        StudyDate = getTableString(M.StudyDate, i);

        age = M.age(i);

        scanner = getTableString(M.manufacturer_model_name, i);
        radionuclide = getTableString(M.pet_radionuclide, i);
        contrast = getTableString(M.ct_contrast_agent, i);

        %% Read images
        CT  = double(niftiread(ctPath));
        PET = double(niftiread(petPath));
        LAB = niftiread(labelPath) > 0;

        info = niftiinfo(labelPath);

        imgSize = info.ImageSize;
        voxelSize = info.PixelDimensions(1:3);
        voxelVolume_mL = prod(voxelSize) / 1000;

        if ~isequal(size(CT), size(PET), size(LAB))
            error("CT, PET, and label images have mismatched dimensions.");
        end

        %% ====================================================
        %  A. TRUE TUMOR CANDIDATES FROM TUMOR MASK
        %  ====================================================

        CCtrue = bwconncomp(LAB, 26);

        nTrueKept = 0;
        nTrueAll = CCtrue.NumObjects;

        for c = 1:CCtrue.NumObjects

            idx = CCtrue.PixelIdxList{c};

            candidateVolume_mL = numel(idx) * voxelVolume_mL;

            if candidateVolume_mL < MIN_CANDIDATE_VOLUME_ML
                continue;
            end

            globalCandidateID = globalCandidateID + 1;
            nTrueKept = nTrueKept + 1;

            row = makeCandidateRow( ...
                globalCandidateID, ...
                case_id, ...
                SubjectID, ...
                StudyDate, ...
                "true_tumor_mask_component", ...
                1, ...
                "true_tumor", ...
                nTrueKept, ...
                idx, ...
                PET, ...
                CT, ...
                imgSize, ...
                voxelSize, ...
                voxelVolume_mL, ...
                age, ...
                scanner, ...
                radionuclide, ...
                contrast);

            allRows(end+1,:) = row;

        end

        %% ====================================================
        %  B. FALSE-POSITIVE HOT PET CANDIDATES
        %  ====================================================

        if DILATE_TUMOR_MASK
            tumorExclusionMask = imdilate(LAB, true(3,3,3));
        else
            tumorExclusionMask = LAB;
        end

        hotMask = isfinite(PET) & PET >= SUV_THRESHOLD;

        % Remove tumor and immediate tumor boundary
        falseHotMask = hotMask & ~tumorExclusionMask;

        CCfalse = bwconncomp(falseHotMask, 26);

        nFalseKept = 0;
        nFalseAll = CCfalse.NumObjects;

        for c = 1:CCfalse.NumObjects

            idx = CCfalse.PixelIdxList{c};

            candidateVolume_mL = numel(idx) * voxelVolume_mL;

            if candidateVolume_mL < MIN_CANDIDATE_VOLUME_ML
                continue;
            end

            globalCandidateID = globalCandidateID + 1;
            nFalseKept = nFalseKept + 1;

            row = makeCandidateRow( ...
                globalCandidateID, ...
                case_id, ...
                SubjectID, ...
                StudyDate, ...
                "false_non_tumor_hot_region", ...
                0, ...
                "false_hot_region", ...
                nFalseKept, ...
                idx, ...
                PET, ...
                CT, ...
                imgSize, ...
                voxelSize, ...
                voxelVolume_mL, ...
                age, ...
                scanner, ...
                radionuclide, ...
                contrast);

            allRows(end+1,:) = row;

        end

        fprintf("True tumor candidates kept:       %d\n", nTrueKept);
        fprintf("False hot-region candidates kept: %d\n", nFalseKept);

        caseSummaryRows(end+1,:) = {
            case_id, ...
            SubjectID, ...
            StudyDate, ...
            age, ...
            scanner, ...
            radionuclide, ...
            contrast, ...
            nTrueAll, ...
            nTrueKept, ...
            nFalseAll, ...
            nFalseKept, ...
            nTrueKept + nFalseKept
            };

    catch ME

        warning("Error in case %d: %s", i, ME.message);

        errorRows(end+1,:) = {
            i, ...
            getTableString(M.case_id, i), ...
            string(ME.message)
            };

    end

end

%% ============================================================
%  CREATE CANDIDATE TABLE
%  ============================================================

varNames = {
    'candidate_global_id'
    'case_id'
    'SubjectID'
    'StudyDate'
    'candidate_source'
    'candidate_label'
    'candidate_label_name'
    'candidate_number_within_source'
    'candidate_voxels'
    'candidate_volume_mL'
    'SUVmax'
    'SUVmean'
    'SUVmedian'
    'SUVstd'
    'CTmean'
    'CTmedian'
    'CTstd'
    'centroid_x_voxel'
    'centroid_y_voxel'
    'centroid_z_voxel'
    'centroid_x_norm'
    'centroid_y_norm'
    'centroid_z_norm'
    'centroid_x_mm'
    'centroid_y_mm'
    'centroid_z_mm'
    'bbox_x_mm'
    'bbox_y_mm'
    'bbox_z_mm'
    'distance_from_center_mm'
    'lateral_zone'
    'AP_zone'
    'craniocaudal_zone'
    'anatomical_spatial_zone'
    'age'
    'scanner'
    'radionuclide'
    'contrast'
    };

Tcand = cell2table(allRows, 'VariableNames', varNames);

%% ============================================================
%  SAVE CANDIDATE TABLE
%  ============================================================

fprintf("\nSaving candidate feature table...\n");

writetable(Tcand, outCandidateXlsx);
writetable(Tcand, outCandidateCsv);

%% ============================================================
%  SUMMARY TABLES
%  ============================================================

totalCandidates = height(Tcand);
trueCandidates = sum(Tcand.candidate_label == 1);
falseCandidates = sum(Tcand.candidate_label == 0);

SummaryMetric = [
    "Total candidates"
    "True tumor candidates"
    "False hot-region candidates"
    "SUV threshold"
    "Minimum candidate volume mL"
    "Tumor mask dilation used"
    ];

SummaryValue = [
    totalCandidates
    trueCandidates
    falseCandidates
    SUV_THRESHOLD
    MIN_CANDIDATE_VOLUME_ML
    double(DILATE_TUMOR_MASK)
    ];

Tsummary = table(SummaryMetric, SummaryValue);

TcaseSummary = cell2table(caseSummaryRows, ...
    'VariableNames', {
    'case_id'
    'SubjectID'
    'StudyDate'
    'age'
    'scanner'
    'radionuclide'
    'contrast'
    'n_true_components_all'
    'n_true_candidates_kept'
    'n_false_hot_components_all'
    'n_false_candidates_kept'
    'n_total_candidates_kept'
    });

if isempty(errorRows)
    Terrors = table([], strings(0,1), strings(0,1), ...
        'VariableNames', {'case_index','case_id','error_message'});
else
    Terrors = cell2table(errorRows, ...
        'VariableNames', {
        'case_index'
        'case_id'
        'error_message'
        });
end

writetable(Tsummary, outSummaryXlsx, "Sheet", "Overall_summary");
writetable(TcaseSummary, outSummaryXlsx, "Sheet", "Case_summary");
writetable(Terrors, outSummaryXlsx, "Sheet", "Errors");

%% ============================================================
%  DISPLAY CLASS BALANCE
%  ============================================================

fprintf("\nObjective 4A candidate generation complete.\n");
fprintf("Total candidates: %d\n", totalCandidates);
fprintf("True tumor candidates: %d\n", trueCandidates);
fprintf("False hot-region candidates: %d\n", falseCandidates);

fprintf("\nSaved candidate table:\n%s\n", outCandidateXlsx);
fprintf("\nSaved summary:\n%s\n", outSummaryXlsx);

%% ============================================================
%  LOCAL FUNCTIONS
%  ============================================================

function s = getTableString(col, i)

    if iscell(col)
        s = string(col{i});
    else
        s = string(col(i));
    end

end

function row = makeCandidateRow( ...
    candidate_global_id, ...
    case_id, ...
    SubjectID, ...
    StudyDate, ...
    candidate_source, ...
    candidate_label, ...
    candidate_label_name, ...
    candidate_number_within_source, ...
    idx, ...
    PET, ...
    CT, ...
    imgSize, ...
    voxelSize, ...
    voxelVolume_mL, ...
    age, ...
    scanner, ...
    radionuclide, ...
    contrast)

    petVals = PET(idx);
    ctVals  = CT(idx);

    candidate_voxels = numel(idx);
    candidate_volume_mL = candidate_voxels * voxelVolume_mL;

    SUVmax = max(petVals);
    SUVmean = mean(petVals, "omitnan");
    SUVmedian = median(petVals, "omitnan");
    SUVstd = std(petVals, "omitnan");

    CTmean = mean(ctVals, "omitnan");
    CTmedian = median(ctVals, "omitnan");
    CTstd = std(ctVals, "omitnan");

    [x, y, z] = ind2sub(imgSize, idx);

    centroid_x_voxel = mean(x);
    centroid_y_voxel = mean(y);
    centroid_z_voxel = mean(z);

    centroid_x_norm = centroid_x_voxel / imgSize(1);
    centroid_y_norm = centroid_y_voxel / imgSize(2);
    centroid_z_norm = centroid_z_voxel / imgSize(3);

    centroid_x_mm = centroid_x_voxel * voxelSize(1);
    centroid_y_mm = centroid_y_voxel * voxelSize(2);
    centroid_z_mm = centroid_z_voxel * voxelSize(3);

    bbox_x_mm = (max(x) - min(x) + 1) * voxelSize(1);
    bbox_y_mm = (max(y) - min(y) + 1) * voxelSize(2);
    bbox_z_mm = (max(z) - min(z) + 1) * voxelSize(3);

    center_x = imgSize(1) / 2;
    center_y = imgSize(2) / 2;
    center_z = imgSize(3) / 2;

    dx_mm = (centroid_x_voxel - center_x) * voxelSize(1);
    dy_mm = (centroid_y_voxel - center_y) * voxelSize(2);
    dz_mm = (centroid_z_voxel - center_z) * voxelSize(3);

    distance_from_center_mm = sqrt(dx_mm^2 + dy_mm^2 + dz_mm^2);

    %% Anatomical localization from confirmed orientation

    % x: right -> left
    if centroid_x_norm < 0.45
        lateral_zone = "right";
    elseif centroid_x_norm <= 0.55
        lateral_zone = "midline";
    else
        lateral_zone = "left";
    end

    % y: posterior -> anterior
    if centroid_y_norm < 0.45
        AP_zone = "posterior";
    elseif centroid_y_norm <= 0.55
        AP_zone = "central";
    else
        AP_zone = "anterior";
    end

    % z: superior -> inferior
    if centroid_z_norm < 1/3
        craniocaudal_zone = "superior";
    elseif centroid_z_norm < 2/3
        craniocaudal_zone = "middle";
    else
        craniocaudal_zone = "inferior";
    end

    anatomical_spatial_zone = lateral_zone + "_" + AP_zone + "_" + craniocaudal_zone;

    row = {
        candidate_global_id, ...
        case_id, ...
        SubjectID, ...
        StudyDate, ...
        candidate_source, ...
        candidate_label, ...
        candidate_label_name, ...
        candidate_number_within_source, ...
        candidate_voxels, ...
        candidate_volume_mL, ...
        SUVmax, ...
        SUVmean, ...
        SUVmedian, ...
        SUVstd, ...
        CTmean, ...
        CTmedian, ...
        CTstd, ...
        centroid_x_voxel, ...
        centroid_y_voxel, ...
        centroid_z_voxel, ...
        centroid_x_norm, ...
        centroid_y_norm, ...
        centroid_z_norm, ...
        centroid_x_mm, ...
        centroid_y_mm, ...
        centroid_z_mm, ...
        bbox_x_mm, ...
        bbox_y_mm, ...
        bbox_z_mm, ...
        distance_from_center_mm, ...
        lateral_zone, ...
        AP_zone, ...
        craniocaudal_zone, ...
        anatomical_spatial_zone, ...
        age, ...
        scanner, ...
        radionuclide, ...
        contrast
        };

end