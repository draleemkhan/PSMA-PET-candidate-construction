function PSMA_objective9B_external_PSMA_overlap_and_labels_ADDENDUM5()
% PSMA Objective 9B - CORRECTED EXTERNAL BRANCH ONLY
%
% Uses PET-derived candidates from corrected Objective 9A external PSMA branch
% and compares them with:
%   <case>\PSMA\TTB.nii.gz
%
% Candidate geometry remains entirely PET-derived.
% Ground truth is used only for overlap measurement and supervised labels.
%
% Outputs:
%   - candidate overlap table
%   - GT lesion coverage table
%   - case summary
%   - split summary
%   - workbook, MAT workspace, report
%
% Resumability:
%   - checkpoint CSVs after every completed case
%   - successful cases skipped on rerun

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 9B - corrected external PSMA overlap analysis\n');
fprintf('============================================================\n\n');

%% ========================================================================
% USER SETTINGS
% ========================================================================

newProjectRoot = ...
    string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));

externalDatasetRoot = psma_root("PSMA_EXTERNAL_DATA");

objective9AExternalDir = fullfile(newProjectRoot,"06_results", ...
    "Objective9_LockedExternal_Addendum_train0096_0100", ...
    "Objective9A_PET_Only_Candidate_Generation", ...
    "external","external");

outDir = fullfile(newProjectRoot,"06_results", ...
    "Objective9_LockedExternal_Addendum_train0096_0100", ...
    "Objective9B_GT_Overlap_and_Label_Assignment", ...
    "external","external");

checkpointDir = fullfile(outDir,"checkpoints");
logDir = fullfile(outDir,"logs");

if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(checkpointDir), mkdir(checkpointDir); end
if ~isfolder(logDir), mkdir(logDir); end

CONNECTIVITY = 26;
overwriteSuccessfulCases = false;
saveExcelOutputs = true;

runStamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

candidateCheckpoint = fullfile(checkpointDir, ...
    "PSMA_objective9B_candidate_overlap_checkpoint.csv");

lesionCheckpoint = fullfile(checkpointDir, ...
    "PSMA_objective9B_GT_lesion_coverage_checkpoint.csv");

caseCheckpoint = fullfile(checkpointDir, ...
    "PSMA_objective9B_case_summary_checkpoint.csv");

diary(fullfile(logDir,"PSMA_objective9B_external_PSMA_only.log"));

%% ========================================================================
% LOAD CORRECTED EXTERNAL OBJECTIVE 9A TABLE
% ========================================================================

candidateIndexPath = latest_valid_candidate_table(objective9AExternalDir);

C = read_csv_strict(candidateIndexPath);

requiredVars = [ ...
    "candidate_global_id","case_id","candidate_label_value", ...
    "candidate_volume_mL","candidate_voxels","PET_path", ...
    "candidate_labelmap_path","voxel_volume_mL"];

must_have(C,requiredVars,"Objective 9A corrected external candidate table");

C.candidate_global_id = string(C.candidate_global_id);
C.case_id = string(C.case_id);
C.candidate_label_value = double(C.candidate_label_value);
C.candidate_volume_mL = double(C.candidate_volume_mL);
C.candidate_voxels = double(C.candidate_voxels);
C.PET_path = string(C.PET_path);
C.candidate_labelmap_path = string(C.candidate_labelmap_path);
C.voxel_volume_mL = double(C.voxel_volume_mL);

% Safety stop: every external PET must come from PSMA subfolder.
badPET = ~contains(upper(C.PET_path),upper(filesep + "PSMA" + filesep));
if any(badPET)
    error('Safety stop: %d candidate rows do not point to PSMA PET paths.',sum(badPET));
end

caseIDs = unique(C.case_id,'stable');

fprintf('Candidate source:\n%s\n',candidateIndexPath);
fprintf('External PSMA cases: %d\n',numel(caseIDs));
fprintf('Candidate rows: %d\n',height(C));

%% ========================================================================
% LOAD OR INITIALIZE CHECKPOINTS
% ========================================================================

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

%% ========================================================================
% PROCESS CASES
% ========================================================================

for i = 1:numel(caseIDs)

    caseID = caseIDs(i);
    fprintf('\n[%d/%d] %s\n',i,numel(caseIDs),caseID);

    if ~overwriteSuccessfulCases && ~isempty(CaseSummary)
        if any(CaseSummary.case_id==caseID & CaseSummary.status=="success")
            fprintf('  Already completed successfully. Skipping.\n');
            continue;
        end
    end

    started = datetime('now');

    try
        X = C(C.case_id==caseID,:);

        labelmapPath = unique(X.candidate_labelmap_path);
        petPath = unique(X.PET_path);

        if numel(labelmapPath)~=1 || numel(petPath)~=1
            error('Inconsistent PET or label-map paths for %s.',caseID);
        end

        % Explicit PSMA-only GT path.
        gtGz = fullfile(externalDatasetRoot,caseID,"PSMA","TTB.nii.gz");
        gtNii = fullfile(externalDatasetRoot,caseID,"PSMA","TTB.nii");

        if isfile(gtGz)
            gtPath = string(gtGz);
        elseif isfile(gtNii)
            gtPath = string(gtNii);
        else
            error('PSMA TTB mask not found for %s.',caseID);
        end

        if ~contains(upper(gtPath),upper(filesep + "PSMA" + filesep))
            error('Safety stop: GT path is not in PSMA folder: %s',gtPath);
        end

        if ~isfile(labelmapPath)
            error('Candidate label map not found: %s',labelmapPath);
        end

        candidateMap = double(niftiread(char(labelmapPath)));
        GT = niftiread(char(gtPath))>0;

        if ~isequal(size(candidateMap),size(GT))
            error('Candidate map and PSMA GT dimensions differ: %s vs %s.', ...
                mat2str(size(candidateMap)),mat2str(size(GT)));
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

        allCandidateMask = candidateMap>0;
        gtCoveredVoxelsByAnyCandidate = zeros(GTcc.NumObjects,1);

        for g = 1:GTcc.NumObjects
            idx = GTcc.PixelIdxList{g};
            gtCoveredVoxelsByAnyCandidate(g) = nnz(allCandidateMask(idx));
        end

        candidateRows = cell(height(X),1);

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
            candidateOverlapFraction = overlapVoxels/candidateVoxelsMeasured;

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
                    coverageValues(q) = nnz(candMask(lesionIdx))/numel(lesionIdx);
                end

                maxGTLesionCoverage = max(coverageValues);
                sumIntersectedGTVolume_mL = sum(gtLesionVolume_mL(intersectedLesions));
                largestIntersectedGTLesion_mL = max(gtLesionVolume_mL(intersectedLesions));
            end

            labelAny = overlapVoxels>0;
            label10 = candidateOverlapFraction>=0.10;
            label25 = candidateOverlapFraction>=0.25;
            label50 = candidateOverlapFraction>=0.50;

            if overlapVoxels==0
                overlapCategory = "no_overlap";
            elseif candidateOverlapFraction<0.10 && maxGTLesionCoverage<0.10
                overlapCategory = "limited_overlap";
            else
                overlapCategory = "substantial_overlap";
            end

            candidateRows{r} = table( ...
                candidateID,"DEEP-PSMA","external","external",caseID, ...
                labelValue,X.candidate_voxels(r),candidateVoxelsMeasured, ...
                X.candidate_volume_mL(r),voxelVolume_mL,string(petPath), ...
                string(labelmapPath),gtPath,overlapVoxels,overlapVolume_mL, ...
                candidateOverlapFraction,intersectedCount,maxGTLesionCoverage, ...
                sumIntersectedGTVolume_mL,largestIntersectedGTLesion_mL, ...
                double(labelAny),double(label10),double(label25),double(label50), ...
                overlapCategory, ...
                'VariableNames',{ ...
                'candidate_global_id','dataset_name','cohort','split','case_id', ...
                'candidate_label_value','candidate_voxels_from_9A', ...
                'candidate_voxels_remeasured','candidate_volume_mL', ...
                'voxel_volume_mL','PET_path','candidate_labelmap_path','GT_path', ...
                'overlap_voxels','overlap_volume_mL','candidate_overlap_fraction', ...
                'intersected_GT_lesion_count','maximum_GT_lesion_coverage', ...
                'summed_intersected_GT_volume_mL','largest_intersected_GT_lesion_mL', ...
                'candidate_label_any_overlap','candidate_label_overlap_10pct', ...
                'candidate_label_overlap_25pct','candidate_label_overlap_50pct', ...
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
                "DEEP-PSMA","external","external",caseID,double(g), ...
                gtLesionVoxelCount(g),gtLesionVolume_mL(g), ...
                double(gtLesionVolume_mL(g)>=1),covered, ...
                lesionCoverageFraction,numel(intersectingCandidateLabels), ...
                double(covered>0),double(lesionCoverageFraction>=0.10), ...
                double(lesionCoverageFraction>=0.25), ...
                double(lesionCoverageFraction>=0.50), ...
                'VariableNames',{ ...
                'dataset_name','cohort','split','case_id','GT_lesion_id', ...
                'GT_lesion_voxels','GT_lesion_volume_mL','robust_lesion_ge1mL', ...
                'covered_voxels_by_any_PET_candidate','GT_lesion_coverage_fraction', ...
                'intersecting_PET_candidate_count','detected_any_overlap', ...
                'detected_10pct','detected_25pct','detected_50pct'});
        end

        if GTcc.NumObjects>0
            caseLesionTable = vertcat(lesionRows{:});
        else
            caseLesionTable = empty_lesion_coverage_table();
        end

        elapsedSeconds = seconds(datetime('now')-started);

        caseSummaryRow = table( ...
            "DEEP-PSMA","external","external",caseID,"success","", ...
            height(X),sum(caseCandidateTable.candidate_label_any_overlap), ...
            sum(caseCandidateTable.candidate_label_any_overlap==0), ...
            sum(caseCandidateTable.candidate_label_overlap_10pct), ...
            sum(caseCandidateTable.candidate_label_overlap_25pct), ...
            sum(caseCandidateTable.candidate_label_overlap_50pct), ...
            GTcc.NumObjects,sum(gtLesionVolume_mL>=1), ...
            sum(caseLesionTable.detected_any_overlap), ...
            sum(caseLesionTable.detected_any_overlap & ...
                caseLesionTable.robust_lesion_ge1mL), ...
            safe_div(sum(caseLesionTable.detected_any_overlap),height(caseLesionTable)), ...
            safe_div(sum(caseLesionTable.detected_any_overlap & ...
                caseLesionTable.robust_lesion_ge1mL), ...
                sum(caseLesionTable.robust_lesion_ge1mL)), ...
            elapsedSeconds, ...
            string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'VariableNames',{ ...
            'dataset_name','cohort','split','case_id','status','error_message', ...
            'PET_candidates','candidates_any_GT_overlap','candidates_no_GT_overlap', ...
            'candidates_overlap_ge10pct','candidates_overlap_ge25pct', ...
            'candidates_overlap_ge50pct','GT_lesions','robust_GT_lesions_ge1mL', ...
            'GT_lesions_detected_any_overlap', ...
            'robust_GT_lesions_detected_any_overlap', ...
            'all_lesion_candidate_recall','robust_lesion_candidate_recall', ...
            'elapsed_seconds','completed_at'});

        if ~isempty(CandidateOverlap)
            CandidateOverlap(CandidateOverlap.case_id==caseID,:) = [];
        end
        CandidateOverlap = [CandidateOverlap; caseCandidateTable];

        if ~isempty(LesionCoverage)
            LesionCoverage(LesionCoverage.case_id==caseID,:) = [];
        end
        LesionCoverage = [LesionCoverage; caseLesionTable];

        CaseSummary = replace_case_summary(CaseSummary,caseSummaryRow);

        writetable(CandidateOverlap,candidateCheckpoint);
        writetable(LesionCoverage,lesionCheckpoint);
        writetable(CaseSummary,caseCheckpoint);

        fprintf('  PET candidates: %d\n',height(X));
        fprintf('  Any-overlap positive candidates: %d\n', ...
            sum(caseCandidateTable.candidate_label_any_overlap));
        fprintf('  No-overlap candidates: %d\n', ...
            sum(caseCandidateTable.candidate_label_any_overlap==0));
        fprintf('  Robust-lesion candidate recall: %.4f\n', ...
            caseSummaryRow.robust_lesion_candidate_recall);

    catch ME
        elapsedSeconds = seconds(datetime('now')-started);

        caseSummaryRow = table( ...
            "DEEP-PSMA","external","external",caseID,"failed",string(ME.message), ...
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            elapsedSeconds,string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'VariableNames',{ ...
            'dataset_name','cohort','split','case_id','status','error_message', ...
            'PET_candidates','candidates_any_GT_overlap','candidates_no_GT_overlap', ...
            'candidates_overlap_ge10pct','candidates_overlap_ge25pct', ...
            'candidates_overlap_ge50pct','GT_lesions','robust_GT_lesions_ge1mL', ...
            'GT_lesions_detected_any_overlap', ...
            'robust_GT_lesions_detected_any_overlap', ...
            'all_lesion_candidate_recall','robust_lesion_candidate_recall', ...
            'elapsed_seconds','completed_at'});

        CaseSummary = replace_case_summary(CaseSummary,caseSummaryRow);
        writetable(CaseSummary,caseCheckpoint);

        fprintf(2,'  FAILED: %s\n',ME.message);
    end
end

%% ========================================================================
% FINAL OUTPUTS
% ========================================================================

candidateCsv = fullfile(outDir, ...
    "PSMA_objective9B_candidate_overlap_table_" + runStamp + ".csv");

lesionCsv = fullfile(outDir, ...
    "PSMA_objective9B_GT_lesion_coverage_" + runStamp + ".csv");

caseCsv = fullfile(outDir, ...
    "PSMA_objective9B_case_summary_" + runStamp + ".csv");

summaryCsv = fullfile(outDir, ...
    "PSMA_objective9B_split_summary_" + runStamp + ".csv");

workbookPath = fullfile(outDir, ...
    "PSMA_objective9B_overlap_QC_" + runStamp + ".xlsx");

matPath = fullfile(outDir, ...
    "PSMA_objective9B_workspace_" + runStamp + ".mat");

writetable(CandidateOverlap,candidateCsv);
writetable(LesionCoverage,lesionCsv);
writetable(CaseSummary,caseCsv);

success = CaseSummary.status=="success";
robust = LesionCoverage.robust_lesion_ge1mL==1;

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
    "External image folder"
    "Primary exported label"
];

SummaryValue = [
    "DEEP-PSMA"
    "external"
    "external"
    string(sum(success))
    string(sum(~success))
    string(height(CandidateOverlap))
    string(sum(CandidateOverlap.candidate_label_any_overlap))
    string(sum(CandidateOverlap.candidate_label_any_overlap==0))
    string(safe_div(sum(CandidateOverlap.candidate_label_any_overlap), ...
        height(CandidateOverlap)))
    string(sum(CandidateOverlap.candidate_label_overlap_10pct))
    string(sum(CandidateOverlap.candidate_label_overlap_25pct))
    string(sum(CandidateOverlap.candidate_label_overlap_50pct))
    string(height(LesionCoverage))
    string(sum(robust))
    string(safe_div(sum(LesionCoverage.detected_any_overlap), ...
        height(LesionCoverage)))
    string(safe_div(sum(LesionCoverage.detected_any_overlap & robust), ...
        sum(robust)))
    string(safe_div(sum(LesionCoverage.detected_10pct & robust),sum(robust)))
    string(safe_div(sum(LesionCoverage.detected_25pct & robust),sum(robust)))
    string(safe_div(sum(LesionCoverage.detected_50pct & robust),sum(robust)))
    "NO"
    "PSMA"
    "candidate_label_any_overlap (provisional)"
];

SplitSummary = table(SummaryMetric,SummaryValue);

writetable(SplitSummary,summaryCsv);

if saveExcelOutputs
    writetable(SplitSummary,workbookPath,'Sheet','Summary');
    writetable(CandidateOverlap,workbookPath,'Sheet','Candidates');
    writetable(LesionCoverage,workbookPath,'Sheet','GT_lesions');
    writetable(CaseSummary,workbookPath,'Sheet','Cases');
end

save(matPath,'CandidateOverlap','LesionCoverage','CaseSummary', ...
    'SplitSummary','candidateIndexPath','CONNECTIVITY','-v7.3');

disp(SplitSummary);

fprintf('\n============================================================\n');
fprintf('CORRECTED EXTERNAL OBJECTIVE 9B COMPLETED\n');
fprintf('============================================================\n');
fprintf('Candidate overlap table:\n%s\n',candidateCsv);
fprintf('GT lesion coverage table:\n%s\n',lesionCsv);
fprintf('Case summary:\n%s\n',caseCsv);
fprintf('Split summary:\n%s\n\n',summaryCsv);

diary off;

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function filePath = latest_valid_candidate_table(folderPath)

D = dir(fullfile(folderPath, ...
    "PSMA_objective9A_PET_only_candidate_index_*.csv"));

if isempty(D)
    error('No corrected Objective 9A candidate table found in:\n%s',folderPath);
end

[~,order] = sort([D.datenum],'descend');

for k = order
    p = string(fullfile(D(k).folder,D(k).name));
    T = read_csv_strict(p);
    if all(ismember(["candidate_global_id","case_id","PET_path"], ...
            string(T.Properties.VariableNames))) && height(T)>0
        filePath = p;
        return;
    end
end

error('No valid corrected Objective 9A candidate table found.');

end

function T = read_csv_strict(filePath)

opts = detectImportOptions(filePath, ...
    'FileType','text','Delimiter',',', ...
    'VariableNamingRule','preserve');

opts.Delimiter = ',';
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];

T = readtable(filePath,opts);

end

function must_have(T,names,source)

available = string(T.Properties.VariableNames);
missing = names(~ismember(names,available));

if ~isempty(missing)
    error('Missing required columns in %s: %s', ...
        source,strjoin(missing,', '));
end

end

function T = replace_case_summary(T,newRow)

if ~isempty(T)
    same = T.case_id==newRow.case_id;
    T(same,:) = [];
end

T = [T; newRow];

end

function value = safe_div(a,b)

if b==0
    value = NaN;
else
    value = a/b;
end

end

function T = empty_candidate_overlap_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    strings(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    strings(0,1), ...
    'VariableNames',{ ...
    'candidate_global_id','dataset_name','cohort','split','case_id', ...
    'candidate_label_value','candidate_voxels_from_9A', ...
    'candidate_voxels_remeasured','candidate_volume_mL','voxel_volume_mL', ...
    'PET_path','candidate_labelmap_path','GT_path','overlap_voxels', ...
    'overlap_volume_mL','candidate_overlap_fraction', ...
    'intersected_GT_lesion_count','maximum_GT_lesion_coverage', ...
    'summed_intersected_GT_volume_mL','largest_intersected_GT_lesion_mL', ...
    'candidate_label_any_overlap','candidate_label_overlap_10pct', ...
    'candidate_label_overlap_25pct','candidate_label_overlap_50pct', ...
    'overlap_category'});

end

function T = empty_lesion_coverage_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{ ...
    'dataset_name','cohort','split','case_id','GT_lesion_id', ...
    'GT_lesion_voxels','GT_lesion_volume_mL','robust_lesion_ge1mL', ...
    'covered_voxels_by_any_PET_candidate','GT_lesion_coverage_fraction', ...
    'intersecting_PET_candidate_count','detected_any_overlap', ...
    'detected_10pct','detected_25pct','detected_50pct'});

end

function T = empty_case_summary_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),strings(0,1), ...
    'VariableNames',{ ...
    'dataset_name','cohort','split','case_id','status','error_message', ...
    'PET_candidates','candidates_any_GT_overlap','candidates_no_GT_overlap', ...
    'candidates_overlap_ge10pct','candidates_overlap_ge25pct', ...
    'candidates_overlap_ge50pct','GT_lesions','robust_GT_lesions_ge1mL', ...
    'GT_lesions_detected_any_overlap', ...
    'robust_GT_lesions_detected_any_overlap', ...
    'all_lesion_candidate_recall','robust_lesion_candidate_recall', ...
    'elapsed_seconds','completed_at'});

end
