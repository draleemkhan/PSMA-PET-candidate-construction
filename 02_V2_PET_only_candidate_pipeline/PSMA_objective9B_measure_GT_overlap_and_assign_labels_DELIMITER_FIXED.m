function PSMA_objective9B_measure_GT_overlap_and_assign_labels()
% PSMA Objective 9B
% Ground-truth overlap measurement for PET-derived candidates generated in 9A.
%
% IMPORTANT SCIENTIFIC DESIGN
%   - Candidate geometry comes entirely from PET-only Objective 9A.
%   - Ground-truth masks are introduced only now.
%   - This script measures overlap and creates transparent candidate labels.
%   - The primary exported binary label is:
%         candidate_label_any_overlap = 1 when candidate overlaps >=1 GT voxel
%         candidate_label_any_overlap = 0 when candidate overlaps 0 GT voxels
%
%   Because a single-voxel overlap may be insufficient for a final locked
%   learning rule, the script also exports:
%       candidate_overlap_fraction
%       overlapping_volume_mL
%       intersected_GT_lesion_count
%       maximum_GT_lesion_coverage
%       candidate_label_overlap_10pct
%       candidate_label_overlap_25pct
%       candidate_label_overlap_50pct
%       overlap_category = no_overlap / limited_overlap / substantial_overlap
%
% The final supervised labeling rule should be selected using TRAINING data
% only, then locked before validation, test, and external analysis.
%
% Inputs
%   Objective 9A candidate label maps and candidate index tables
%   Internal ground-truth masks from existing physical splits
%   External DEEP-PSMA TTB masks for locked train_0001 ... train_0095
%
% Outputs
%   Candidate-level overlap table
%   GT lesion-level coverage table
%   Case-level summary
%   Split-level summary
%   QC workbook, MAT workspace, and text report
%
% Resumability
%   Candidate and case results are checkpointed after every completed case.

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 9B: GT overlap measurement and candidate labels\n');
fprintf('Candidate geometry remains PET-derived from Objective 9A\n');
fprintf('============================================================\n\n');

%% ========================================================================
% USER SETTINGS
% ========================================================================

sourceProjectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));
newProjectRoot = ...
    string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));

internalSplitRoot = fullfile(sourceProjectRoot,"00_data_split_physical");
externalDatasetRoot = psma_root("PSMA_EXTERNAL_DATA");

objective9ARoot = fullfile(newProjectRoot,"06_results", ...
    "Objective9_PET_Derived_Candidate_Pipeline", ...
    "Objective9A_PET_Only_Candidate_Generation");

outRoot = fullfile(newProjectRoot,"06_results", ...
    "Objective9_PET_Derived_Candidate_Pipeline", ...
    "Objective9B_GT_Overlap_and_Label_Assignment");

if ~isfolder(outRoot), mkdir(outRoot); end

CONNECTIVITY = 26;
overwriteSuccessfulCases = false;
saveExcelOutputs = true;

runStamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

%% ========================================================================
% DATASET CONFIGURATIONS
% ========================================================================

configs = {
    struct('dataset_name',"PSMA-PET-CT-Lesions", ...
           'cohort',"internal",'split',"train", ...
           'mode',"internal_split")
    struct('dataset_name',"PSMA-PET-CT-Lesions", ...
           'cohort',"internal",'split',"validation", ...
           'mode',"internal_split")
    struct('dataset_name',"PSMA-PET-CT-Lesions", ...
           'cohort',"internal",'split',"test", ...
           'mode',"internal_split")
    struct('dataset_name',"DEEP-PSMA", ...
           'cohort',"external",'split',"external", ...
           'mode',"external_locked95")
    };

AllCandidateOverlap = table();
AllLesionCoverage = table();
AllCaseSummary = table();
AllSplitSummary = table();

%% ========================================================================
% PROCESS EACH COHORT/SPLIT
% ========================================================================

for d = 1:numel(configs)

    cfg = configs{d};

    fprintf('\n============================================================\n');
    fprintf('DATASET: %s | COHORT: %s | SPLIT: %s\n', ...
        cfg.dataset_name,cfg.cohort,cfg.split);
    fprintf('============================================================\n');

    inputDir = fullfile(objective9ARoot,cfg.cohort,cfg.split);
    if ~isfolder(inputDir)
        error('Objective 9A folder not found:\n%s',inputDir);
    end

    requiredCandidateVars = [ ...
        "candidate_global_id","case_id","candidate_label_value", ...
        "candidate_volume_mL","candidate_voxels","PET_path", ...
        "candidate_labelmap_path","voxel_volume_mL"];

    % Prefer the finalized Objective 9A table, but fall back to the
    % resumable checkpoint if the final file is absent or malformed.
    candidateIndexPath = locate_valid_9A_candidate_table( ...
        inputDir,requiredCandidateVars);

    C = read_csv_strict(candidateIndexPath);

    C.Properties.VariableNames = matlab.lang.makeValidName( ...
        C.Properties.VariableNames);

    must_have(C,requiredCandidateVars,"Objective 9A candidate table");

    C.candidate_global_id = string(C.candidate_global_id);
    C.case_id = string(C.case_id);
    C.candidate_label_value = double(C.candidate_label_value);
    C.candidate_volume_mL = double(C.candidate_volume_mL);
    C.candidate_voxels = double(C.candidate_voxels);
    C.PET_path = string(C.PET_path);
    C.candidate_labelmap_path = string(C.candidate_labelmap_path);
    C.voxel_volume_mL = double(C.voxel_volume_mL);

    caseIDs = unique(C.case_id,'stable');

    outDir = fullfile(outRoot,cfg.cohort,cfg.split);
    checkpointDir = fullfile(outDir,"checkpoints");
    logDir = fullfile(outDir,"logs");

    if ~isfolder(outDir), mkdir(outDir); end
    if ~isfolder(checkpointDir), mkdir(checkpointDir); end
    if ~isfolder(logDir), mkdir(logDir); end

    candidateCheckpoint = fullfile(checkpointDir, ...
        "PSMA_objective9B_candidate_overlap_checkpoint.csv");

    lesionCheckpoint = fullfile(checkpointDir, ...
        "PSMA_objective9B_GT_lesion_coverage_checkpoint.csv");

    caseCheckpoint = fullfile(checkpointDir, ...
        "PSMA_objective9B_case_summary_checkpoint.csv");

    if isfile(candidateCheckpoint)
        CandidateOverlap = read_csv_strict(candidateCheckpoint);
    else
        CandidateOverlap = empty_candidate_overlap_table();
    end

    if isfile(lesionCheckpoint)
        LesionCoverage = read_csv_strict(lesionCheckpoint);
    else
        LesionCoverage = empty_lesion_coverage_table();
    end

    if isfile(caseCheckpoint)
        CaseSummary = read_csv_strict(caseCheckpoint);
    else
        CaseSummary = empty_case_summary_table();
    end

    logFile = fullfile(logDir, ...
        "PSMA_objective9B_" + cfg.cohort + "_" + cfg.split + ".log");
    diary(logFile);

    fprintf('Candidate source:\n%s\n',candidateIndexPath);
    fprintf('Unique cases: %d\n',numel(caseIDs));

    for i = 1:numel(caseIDs)

        caseID = caseIDs(i);

        fprintf('\n[%d/%d] %s\n',i,numel(caseIDs),caseID);

        completed = false;
        if ~overwriteSuccessfulCases && ~isempty(CaseSummary)
            completed = any(CaseSummary.case_id==caseID & ...
                CaseSummary.status=="success");
        end

        if completed
            fprintf('  Already completed successfully. Skipping.\n');
            continue;
        end

        started = datetime('now');

        try
            X = C(C.case_id==caseID,:);

            labelmapPath = unique(X.candidate_labelmap_path);
            petPath = unique(X.PET_path);

            if numel(labelmapPath)~=1 || numel(petPath)~=1
                error('Inconsistent candidate-labelmap or PET paths.');
            end

            gtPath = resolve_gt_path(cfg,caseID, ...
                internalSplitRoot,externalDatasetRoot);

            if ~isfile(labelmapPath)
                error('Candidate label map not found: %s',labelmapPath);
            end

            if ~isfile(gtPath)
                error('Ground-truth mask not found: %s',gtPath);
            end

            labelInfo = niftiinfo(char(labelmapPath));
            gtInfo = niftiinfo(char(gtPath));

            candidateMap = double(niftiread(char(labelmapPath)));
            GT = niftiread(char(gtPath))>0;

            if ~isequal(size(candidateMap),size(GT))
                error('Candidate-map and GT dimensions differ: %s vs %s.', ...
                    mat2str(size(candidateMap)),mat2str(size(GT)));
            end

            if ~isequal(double(labelInfo.ImageSize),double(gtInfo.ImageSize))
                error('Candidate-map and GT NIfTI image sizes differ.');
            end

            voxelVolume_mL = median(double(X.voxel_volume_mL),'omitnan');

            GTcc = bwconncomp(GT,CONNECTIVITY);
            GTlesionMap = zeros(size(GT),'uint32');

            gtLesionVolume_mL = zeros(GTcc.NumObjects,1);
            gtLesionVoxelCount = zeros(GTcc.NumObjects,1);

            for g = 1:GTcc.NumObjects
                idx = GTcc.PixelIdxList{g};
                GTlesionMap(idx) = uint32(g);
                gtLesionVoxelCount(g) = numel(idx);
                gtLesionVolume_mL(g) = numel(idx)*voxelVolume_mL;
            end

            candidateRows = cell(height(X),1);

            % Track lesion-level coverage by union of all PET candidates.
            allCandidateMask = candidateMap>0;
            gtCoveredVoxelsByAnyCandidate = zeros(GTcc.NumObjects,1);

            for g = 1:GTcc.NumObjects
                idx = GTcc.PixelIdxList{g};
                gtCoveredVoxelsByAnyCandidate(g) = nnz(allCandidateMask(idx));
            end

            for r = 1:height(X)

                candidateID = X.candidate_global_id(r);
                labelValue = X.candidate_label_value(r);

                candMask = candidateMap==labelValue;
                candidateVoxelsMeasured = nnz(candMask);

                if candidateVoxelsMeasured==0
                    error('Candidate label %d has zero voxels.',labelValue);
                end

                overlapMask = candMask & GT;
                overlapVoxels = nnz(overlapMask);
                overlapVolume_mL = overlapVoxels*voxelVolume_mL;
                candidateOverlapFraction = ...
                    overlapVoxels/candidateVoxelsMeasured;

                intersectedLesions = unique(double(GTlesionMap(overlapMask)));
                intersectedLesions(intersectedLesions==0) = [];
                intersectedCount = numel(intersectedLesions);

                maxGTLesionCoverage = 0;
                sumIntersectedGTVolume_mL = 0;
                largestIntersectedGTLesion_mL = 0;

                if intersectedCount>0
                    coverageValues = zeros(intersectedCount,1);

                    for q = 1:intersectedCount
                        lesionID = intersectedLesions(q);
                        lesionIdx = GTcc.PixelIdxList{lesionID};
                        coverageValues(q) = nnz(candMask(lesionIdx)) / ...
                            numel(lesionIdx);
                    end

                    maxGTLesionCoverage = max(coverageValues);
                    sumIntersectedGTVolume_mL = ...
                        sum(gtLesionVolume_mL(intersectedLesions));
                    largestIntersectedGTLesion_mL = ...
                        max(gtLesionVolume_mL(intersectedLesions));
                end

                labelAny = overlapVoxels>0;
                label10 = candidateOverlapFraction>=0.10;
                label25 = candidateOverlapFraction>=0.25;
                label50 = candidateOverlapFraction>=0.50;

                if overlapVoxels==0
                    overlapCategory = "no_overlap";
                elseif candidateOverlapFraction<0.10 && ...
                        maxGTLesionCoverage<0.10
                    overlapCategory = "limited_overlap";
                else
                    overlapCategory = "substantial_overlap";
                end

                candidateRows{r} = table( ...
                    candidateID,string(cfg.dataset_name), ...
                    string(cfg.cohort),string(cfg.split),caseID, ...
                    labelValue,X.candidate_voxels(r), ...
                    candidateVoxelsMeasured,X.candidate_volume_mL(r), ...
                    voxelVolume_mL,string(petPath),string(labelmapPath), ...
                    string(gtPath),overlapVoxels,overlapVolume_mL, ...
                    candidateOverlapFraction,intersectedCount, ...
                    maxGTLesionCoverage,sumIntersectedGTVolume_mL, ...
                    largestIntersectedGTLesion_mL, ...
                    double(labelAny),double(label10),double(label25), ...
                    double(label50),overlapCategory, ...
                    'VariableNames',{ ...
                    'candidate_global_id','dataset_name','cohort','split', ...
                    'case_id','candidate_label_value', ...
                    'candidate_voxels_from_9A', ...
                    'candidate_voxels_remeasured', ...
                    'candidate_volume_mL','voxel_volume_mL', ...
                    'PET_path','candidate_labelmap_path','GT_path', ...
                    'overlap_voxels','overlap_volume_mL', ...
                    'candidate_overlap_fraction', ...
                    'intersected_GT_lesion_count', ...
                    'maximum_GT_lesion_coverage', ...
                    'summed_intersected_GT_volume_mL', ...
                    'largest_intersected_GT_lesion_mL', ...
                    'candidate_label_any_overlap', ...
                    'candidate_label_overlap_10pct', ...
                    'candidate_label_overlap_25pct', ...
                    'candidate_label_overlap_50pct', ...
                    'overlap_category'});
            end

            caseCandidateTable = vertcat(candidateRows{:});

            lesionRows = cell(GTcc.NumObjects,1);

            for g = 1:GTcc.NumObjects
                idx = GTcc.PixelIdxList{g};
                covered = gtCoveredVoxelsByAnyCandidate(g);
                lesionCoverageFraction = covered/numel(idx);

                intersectingCandidateLabels = unique(candidateMap(idx));
                intersectingCandidateLabels(intersectingCandidateLabels==0) = [];

                lesionRows{g} = table( ...
                    string(cfg.dataset_name),string(cfg.cohort), ...
                    string(cfg.split),caseID,double(g), ...
                    gtLesionVoxelCount(g),gtLesionVolume_mL(g), ...
                    double(gtLesionVolume_mL(g)>=1), ...
                    covered,lesionCoverageFraction, ...
                    numel(intersectingCandidateLabels), ...
                    double(covered>0), ...
                    double(lesionCoverageFraction>=0.10), ...
                    double(lesionCoverageFraction>=0.25), ...
                    double(lesionCoverageFraction>=0.50), ...
                    'VariableNames',{ ...
                    'dataset_name','cohort','split','case_id', ...
                    'GT_lesion_id','GT_lesion_voxels', ...
                    'GT_lesion_volume_mL','robust_lesion_ge1mL', ...
                    'covered_voxels_by_any_PET_candidate', ...
                    'GT_lesion_coverage_fraction', ...
                    'intersecting_PET_candidate_count', ...
                    'detected_any_overlap','detected_10pct', ...
                    'detected_25pct','detected_50pct'});
            end

            if GTcc.NumObjects>0
                caseLesionTable = vertcat(lesionRows{:});
            else
                caseLesionTable = empty_lesion_coverage_table();
            end

            elapsedSeconds = seconds(datetime('now')-started);

            caseSummaryRow = table( ...
                string(cfg.dataset_name),string(cfg.cohort), ...
                string(cfg.split),caseID,"success","", ...
                height(X),sum(caseCandidateTable.candidate_label_any_overlap), ...
                sum(caseCandidateTable.candidate_label_any_overlap==0), ...
                sum(caseCandidateTable.candidate_label_overlap_10pct), ...
                sum(caseCandidateTable.candidate_label_overlap_25pct), ...
                sum(caseCandidateTable.candidate_label_overlap_50pct), ...
                GTcc.NumObjects, ...
                sum(gtLesionVolume_mL>=1), ...
                sum(caseLesionTable.detected_any_overlap), ...
                sum(caseLesionTable.detected_any_overlap & ...
                    caseLesionTable.robust_lesion_ge1mL), ...
                safe_div(sum(caseLesionTable.detected_any_overlap), ...
                    height(caseLesionTable)), ...
                safe_div(sum(caseLesionTable.detected_any_overlap & ...
                    caseLesionTable.robust_lesion_ge1mL), ...
                    sum(caseLesionTable.robust_lesion_ge1mL)), ...
                elapsedSeconds, ...
                string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
                'VariableNames',{ ...
                'dataset_name','cohort','split','case_id', ...
                'status','error_message','PET_candidates', ...
                'candidates_any_GT_overlap','candidates_no_GT_overlap', ...
                'candidates_overlap_ge10pct', ...
                'candidates_overlap_ge25pct', ...
                'candidates_overlap_ge50pct', ...
                'GT_lesions','robust_GT_lesions_ge1mL', ...
                'GT_lesions_detected_any_overlap', ...
                'robust_GT_lesions_detected_any_overlap', ...
                'all_lesion_candidate_recall', ...
                'robust_lesion_candidate_recall', ...
                'elapsed_seconds','completed_at'});

            % Replace any previous rows for this case.
            if ~isempty(CandidateOverlap)
                CandidateOverlap(CandidateOverlap.case_id==caseID,:) = [];
            end
            CandidateOverlap = [CandidateOverlap; caseCandidateTable]; %#ok<AGROW>

            if ~isempty(LesionCoverage)
                LesionCoverage(LesionCoverage.case_id==caseID,:) = [];
            end
            LesionCoverage = [LesionCoverage; caseLesionTable]; %#ok<AGROW>

            CaseSummary = replace_case_summary(CaseSummary,caseSummaryRow);

            writetable(CandidateOverlap,candidateCheckpoint);
            writetable(LesionCoverage,lesionCheckpoint);
            writetable(CaseSummary,caseCheckpoint);

            fprintf('  PET candidates: %d\n',height(X));
            fprintf('  Any-overlap positive candidates: %d\n', ...
                sum(caseCandidateTable.candidate_label_any_overlap));
            fprintf('  No-overlap candidates: %d\n', ...
                sum(caseCandidateTable.candidate_label_any_overlap==0));
            fprintf('  GT lesions: %d | Robust >=1 mL: %d\n', ...
                GTcc.NumObjects,sum(gtLesionVolume_mL>=1));
            fprintf('  Robust-lesion candidate recall: %.4f\n', ...
                caseSummaryRow.robust_lesion_candidate_recall);

        catch ME
            elapsedSeconds = seconds(datetime('now')-started);

            caseSummaryRow = table( ...
                string(cfg.dataset_name),string(cfg.cohort), ...
                string(cfg.split),caseID,"failed",string(ME.message), ...
                NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
                NaN,NaN,elapsedSeconds, ...
                string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
                'VariableNames',{ ...
                'dataset_name','cohort','split','case_id', ...
                'status','error_message','PET_candidates', ...
                'candidates_any_GT_overlap','candidates_no_GT_overlap', ...
                'candidates_overlap_ge10pct', ...
                'candidates_overlap_ge25pct', ...
                'candidates_overlap_ge50pct', ...
                'GT_lesions','robust_GT_lesions_ge1mL', ...
                'GT_lesions_detected_any_overlap', ...
                'robust_GT_lesions_detected_any_overlap', ...
                'all_lesion_candidate_recall', ...
                'robust_lesion_candidate_recall', ...
                'elapsed_seconds','completed_at'});

            CaseSummary = replace_case_summary(CaseSummary,caseSummaryRow);
            writetable(CaseSummary,caseCheckpoint);

            fprintf(2,'  FAILED: %s\n',ME.message);
        end
    end

    %% Final split outputs
    candidateCsv = fullfile(outDir, ...
        "PSMA_objective9B_candidate_overlap_table_" + runStamp + ".csv");

    lesionCsv = fullfile(outDir, ...
        "PSMA_objective9B_GT_lesion_coverage_" + runStamp + ".csv");

    caseCsv = fullfile(outDir, ...
        "PSMA_objective9B_case_summary_" + runStamp + ".csv");

    workbookPath = fullfile(outDir, ...
        "PSMA_objective9B_overlap_QC_" + runStamp + ".xlsx");

    matPath = fullfile(outDir, ...
        "PSMA_objective9B_workspace_" + runStamp + ".mat");

    reportPath = fullfile(outDir, ...
        "PSMA_objective9B_report_" + runStamp + ".txt");

    writetable(CandidateOverlap,candidateCsv);
    writetable(LesionCoverage,lesionCsv);
    writetable(CaseSummary,caseCsv);

    SplitSummary = make_split_summary(cfg,CandidateOverlap, ...
        LesionCoverage,CaseSummary);

    summaryCsv = fullfile(outDir, ...
        "PSMA_objective9B_split_summary_" + runStamp + ".csv");

    writetable(SplitSummary,summaryCsv);

    if saveExcelOutputs
        writetable(SplitSummary,workbookPath,'Sheet','Summary');
        writetable(CandidateOverlap,workbookPath,'Sheet','Candidates');
        writetable(LesionCoverage,workbookPath,'Sheet','GT_lesions');
        writetable(CaseSummary,workbookPath,'Sheet','Cases');
    end

    save(matPath,'CandidateOverlap','LesionCoverage','CaseSummary', ...
        'SplitSummary','cfg','candidateIndexPath','CONNECTIVITY','-v7.3');

    write_report(reportPath,SplitSummary);

    disp(SplitSummary);

    fprintf('\nCandidate overlap table:\n%s\n',candidateCsv);
    fprintf('GT lesion coverage table:\n%s\n',lesionCsv);
    fprintf('Case summary:\n%s\n',caseCsv);
    fprintf('Split summary:\n%s\n',summaryCsv);
    fprintf('Workbook:\n%s\n',workbookPath);
    fprintf('MAT workspace:\n%s\n',matPath);
    fprintf('Report:\n%s\n',reportPath);

    AllCandidateOverlap = [AllCandidateOverlap; CandidateOverlap]; %#ok<AGROW>
    AllLesionCoverage = [AllLesionCoverage; LesionCoverage]; %#ok<AGROW>
    AllCaseSummary = [AllCaseSummary; CaseSummary]; %#ok<AGROW>
    AllSplitSummary = [AllSplitSummary; SplitSummary]; %#ok<AGROW>

    diary off;
end

%% ========================================================================
% COMBINED OUTPUTS
% ========================================================================

combinedDir = fullfile(outRoot,"combined");
if ~isfolder(combinedDir), mkdir(combinedDir); end

allCandidateCsv = fullfile(combinedDir, ...
    "PSMA_objective9B_all_candidate_overlap_" + runStamp + ".csv");

allLesionCsv = fullfile(combinedDir, ...
    "PSMA_objective9B_all_GT_lesion_coverage_" + runStamp + ".csv");

allCaseCsv = fullfile(combinedDir, ...
    "PSMA_objective9B_all_case_summary_" + runStamp + ".csv");

allSummaryCsv = fullfile(combinedDir, ...
    "PSMA_objective9B_all_split_summary_" + runStamp + ".csv");

writetable(AllCandidateOverlap,allCandidateCsv);
writetable(AllLesionCoverage,allLesionCsv);
writetable(AllCaseSummary,allCaseCsv);
writetable(AllSplitSummary,allSummaryCsv);

save(fullfile(combinedDir, ...
    "PSMA_objective9B_all_workspace_" + runStamp + ".mat"), ...
    'AllCandidateOverlap','AllLesionCoverage', ...
    'AllCaseSummary','AllSplitSummary','-v7.3');

fprintf('\n============================================================\n');
fprintf('OBJECTIVE 9B COMPLETED\n');
fprintf('============================================================\n');
fprintf('Combined candidate overlap table:\n%s\n',allCandidateCsv);
fprintf('Combined GT lesion coverage table:\n%s\n',allLesionCsv);
fprintf('Combined case summary:\n%s\n',allCaseCsv);
fprintf('Combined split summary:\n%s\n\n',allSummaryCsv);

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function gtPath = resolve_gt_path(cfg,caseID,internalSplitRoot,externalRoot)

switch cfg.mode
    case "internal_split"
        labelsDir = fullfile(internalSplitRoot,cfg.split,"labelsTr");

        candidates = [
            fullfile(labelsDir,caseID+".nii.gz")
            fullfile(labelsDir,caseID+".nii")
        ];

        gtPath = first_existing(candidates);

    case "external_locked95"
        caseFolder = fullfile(externalRoot,caseID);
        gtFiles = dir(fullfile(caseFolder,"**","TTB.nii.gz"));

        if isempty(gtFiles)
            gtFiles = dir(fullfile(caseFolder,"**","TTB.nii"));
        end

        if isempty(gtFiles)
            error('No TTB.nii(.gz) found under %s.',caseFolder);
        end

        validPaths = strings(0,1);

        for i = 1:numel(gtFiles)
            folder = string(gtFiles(i).folder);

            hasPET = isfile(fullfile(folder,"PET.nii.gz")) || ...
                isfile(fullfile(folder,"PET.nii"));

            hasCT = isfile(fullfile(folder,"CT.nii.gz")) || ...
                isfile(fullfile(folder,"CT.nii"));

            if hasPET && hasCT
                validPaths(end+1,1) = string(fullfile( ...
                    gtFiles(i).folder,gtFiles(i).name)); %#ok<AGROW>
            end
        end

        if isempty(validPaths)
            error(['TTB copies were found, but none shared a folder with ' ...
                'the exact PET/CT trio under %s.'],caseFolder);
        end

        pathLengths = strlength(validPaths);
        shortest = sort(validPaths(pathLengths==min(pathLengths)));
        gtPath = shortest(1);

    otherwise
        error('Unknown GT-resolution mode: %s',cfg.mode);
end

end

function path = first_existing(candidates)

path = "";

for i = 1:numel(candidates)
    if isfile(candidates(i))
        path = candidates(i);
        return;
    end
end

error('None of the expected files exists:\n%s',strjoin(candidates,newline));

end

function T = replace_case_summary(T,newRow)

if ~isempty(T)
    same = T.cohort==newRow.cohort & ...
        T.split==newRow.split & ...
        T.case_id==newRow.case_id;
    T(same,:) = [];
end

T = [T; newRow];

end

function S = make_split_summary(cfg,C,L,CaseSummary)

success = CaseSummary.status=="success";
robust = L.robust_lesion_ge1mL==1;

SummaryMetric = [
    "Dataset"
    "Cohort"
    "Split"
    "Cases passed"
    "Cases failed"
    "PET-derived candidates"
    "Candidates with any GT overlap"
    "Candidates with no GT overlap"
    "Candidate positive fraction by any overlap"
    "Candidates with >=10% candidate overlap"
    "Candidates with >=25% candidate overlap"
    "Candidates with >=50% candidate overlap"
    "Reference GT lesions"
    "Robust GT lesions >=1 mL"
    "All-lesion candidate recall, any overlap"
    "Robust-lesion candidate recall, any overlap"
    "Robust-lesion candidate recall, >=10% coverage"
    "Robust-lesion candidate recall, >=25% coverage"
    "Robust-lesion candidate recall, >=50% coverage"
    "Ground truth used to generate candidates"
    "Primary exported label"
];

SummaryValue = [
    string(cfg.dataset_name)
    string(cfg.cohort)
    string(cfg.split)
    string(sum(success))
    string(sum(~success))
    string(height(C))
    string(sum(C.candidate_label_any_overlap))
    string(sum(C.candidate_label_any_overlap==0))
    string(safe_div(sum(C.candidate_label_any_overlap),height(C)))
    string(sum(C.candidate_label_overlap_10pct))
    string(sum(C.candidate_label_overlap_25pct))
    string(sum(C.candidate_label_overlap_50pct))
    string(height(L))
    string(sum(robust))
    string(safe_div(sum(L.detected_any_overlap),height(L)))
    string(safe_div(sum(L.detected_any_overlap & robust),sum(robust)))
    string(safe_div(sum(L.detected_10pct & robust),sum(robust)))
    string(safe_div(sum(L.detected_25pct & robust),sum(robust)))
    string(safe_div(sum(L.detected_50pct & robust),sum(robust)))
    "NO"
    "candidate_label_any_overlap (provisional)"
];

S = table(SummaryMetric,SummaryValue);

end

function T = empty_candidate_overlap_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{ ...
    'candidate_global_id','dataset_name','cohort','split', ...
    'case_id','candidate_label_value','candidate_voxels_from_9A', ...
    'candidate_voxels_remeasured','candidate_volume_mL', ...
    'voxel_volume_mL','PET_path','candidate_labelmap_path','GT_path', ...
    'overlap_voxels','overlap_volume_mL', ...
    'candidate_overlap_fraction','intersected_GT_lesion_count', ...
    'maximum_GT_lesion_coverage','summed_intersected_GT_volume_mL', ...
    'largest_intersected_GT_lesion_mL', ...
    'candidate_label_any_overlap','candidate_label_overlap_10pct', ...
    'candidate_label_overlap_25pct','candidate_label_overlap_50pct', ...
    'overlap_category'});

end

function T = empty_lesion_coverage_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{ ...
    'dataset_name','cohort','split','case_id','GT_lesion_id', ...
    'GT_lesion_voxels','GT_lesion_volume_mL', ...
    'robust_lesion_ge1mL','covered_voxels_by_any_PET_candidate', ...
    'GT_lesion_coverage_fraction','intersecting_PET_candidate_count', ...
    'detected_any_overlap','detected_10pct','detected_25pct', ...
    'detected_50pct'});

end

function T = empty_case_summary_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),strings(0,1), ...
    'VariableNames',{ ...
    'dataset_name','cohort','split','case_id','status', ...
    'error_message','PET_candidates','candidates_any_GT_overlap', ...
    'candidates_no_GT_overlap','candidates_overlap_ge10pct', ...
    'candidates_overlap_ge25pct','candidates_overlap_ge50pct', ...
    'GT_lesions','robust_GT_lesions_ge1mL', ...
    'GT_lesions_detected_any_overlap', ...
    'robust_GT_lesions_detected_any_overlap', ...
    'all_lesion_candidate_recall', ...
    'robust_lesion_candidate_recall','elapsed_seconds','completed_at'});

end

function filePath = locate_valid_9A_candidate_table(inputDir,requiredVars)
% Locate a valid Objective 9A candidate table.
% Search finalized tables first, then the per-split checkpoint.

patterns = [
    "PSMA_objective9A_PET_only_candidate_index_*.csv"
    fullfile("checkpoints","PSMA_objective9A_candidate_index_checkpoint.csv")
];

candidates = strings(0,1);

% Final timestamped files.
D = dir(fullfile(inputDir,patterns(1)));
if ~isempty(D)
    [~,order] = sort([D.datenum],'descend');
    for i = order
        candidates(end+1,1) = string(fullfile(D(i).folder,D(i).name)); %#ok<AGROW>
    end
end

% Checkpoint fallback.
checkpointPath = fullfile(inputDir,patterns(2));
if isfile(checkpointPath)
    candidates(end+1,1) = string(checkpointPath); %#ok<AGROW>
end

if isempty(candidates)
    error(['No Objective 9A candidate table or checkpoint found in:\n%s\n' ...
        'Run/complete Objective 9A first.'],inputDir);
end

diagnostics = strings(0,1);

for i = 1:numel(candidates)
    p = candidates(i);

    try
        T = read_csv_strict(p);

        vars = string(matlab.lang.makeValidName(T.Properties.VariableNames));
        missing = requiredVars(~ismember(requiredVars,vars));

        if isempty(missing) && height(T)>0
            filePath = p;
            fprintf('Using Objective 9A candidate table:\n%s\n',filePath);
            fprintf('Rows: %d | Columns: %d\n',height(T),width(T));
            return;
        end

        diagnostics(end+1,1) = ...
            p + " | rows=" + height(T) + ...
            " | columns=" + strjoin(vars,", ") + ...
            " | missing=" + strjoin(missing,", "); %#ok<AGROW>

    catch ME
        diagnostics(end+1,1) = p + " | read error: " + string(ME.message); %#ok<AGROW>
    end
end

error(['Objective 9A files were found, but none contained the required ' ...
    'candidate columns.\n\nInspected files:\n%s\n\n' ...
    'This usually means Objective 9A has not completed the split, or an ' ...
    'older malformed file is being selected.'],strjoin(diagnostics,newline));

end

function T = read_csv_strict(filePath)
% Read Objective 9 CSV files using an explicitly locked comma delimiter.
% MATLAB automatic delimiter detection may incorrectly interpret underscores
% in case IDs and Windows paths as delimiters.

opts = detectImportOptions(filePath, ...
    'FileType','text', ...
    'Delimiter',',', ...
    'VariableNamingRule','preserve');

opts.Delimiter = ',';
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];

% Preserve identifier/path/status columns as strings while retaining
% numeric inference for quantitative columns.
for i = 1:numel(opts.VariableNames)
    vn = lower(string(opts.VariableNames{i}));
    if contains(vn,"id") || contains(vn,"path") || ...
            contains(vn,"name") || contains(vn,"dataset") || ...
            contains(vn,"cohort") || contains(vn,"split") || ...
            contains(vn,"status") || contains(vn,"message") || ...
            contains(vn,"category") || contains(vn,"completed")
        try
            opts = setvartype(opts,opts.VariableNames{i},'string');
        catch
        end
    end
end

T = readtable(filePath,opts);

end

function filePath = latest_file(folderPath,pattern)

D = dir(fullfile(folderPath,pattern));

if isempty(D)
    error('No files matching %s in:\n%s',pattern,folderPath);
end

[~,idx] = max([D.datenum]);
filePath = string(fullfile(D(idx).folder,D(idx).name));

end

function must_have(T,names,source)

available = string(T.Properties.VariableNames);
missing = names(~ismember(names,available));

if ~isempty(missing)
    error('Missing required columns in %s: %s', ...
        source,strjoin(missing,', '));
end

end

function value = safe_div(a,b)

if b==0
    value = NaN;
else
    value = a/b;
end

end

function write_report(path,Summary)

fid = fopen(path,'w');

if fid<0
    warning('Could not create report: %s',path);
    return;
end

fprintf(fid,'PSMA Objective 9B: GT overlap measurement\n');
fprintf(fid,'Generated: %s\n\n',string(datetime('now')));

for i = 1:height(Summary)
    fprintf(fid,'%-62s %s\n', ...
        Summary.SummaryMetric(i),Summary.SummaryValue(i));
end

fclose(fid);

end
