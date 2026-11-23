clear; clc; close all;

%% ============================================================
%  OBJECTIVE 4B
%  PATIENT-BASED TRAIN / VALIDATION / TEST SPLIT
%
%  Split unit:
%    SubjectID
%
%  This prevents leakage from the same patient appearing in both
%  training and testing.
%
%  Output:
%    candidate table with split labels
%    train table
%    validation table
%    test table
%    split summary
%  ============================================================

projectDir = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));

inputFile = fullfile(projectDir, ...
    "06_results", ...
    "Objective4_TRUE_vs_FALSE_candidate_dataset", ...
    "PSMA_objective4_candidate_features_TRUE_vs_FALSE.xlsx");

resultDir = fullfile(projectDir, ...
    "06_results", ...
    "Objective4_TRUE_vs_FALSE_candidate_dataset", ...
    "Patient_based_split");

if ~exist(resultDir, "dir")
    mkdir(resultDir);
end

outAllXlsx = fullfile(resultDir, ...
    "PSMA_objective4_candidates_with_patient_split.xlsx");

outAllCsv = fullfile(resultDir, ...
    "PSMA_objective4_candidates_with_patient_split.csv");

outTrainXlsx = fullfile(resultDir, ...
    "PSMA_objective4_TRAIN_candidates.xlsx");

outValXlsx = fullfile(resultDir, ...
    "PSMA_objective4_VALIDATION_candidates.xlsx");

outTestXlsx = fullfile(resultDir, ...
    "PSMA_objective4_TEST_candidates.xlsx");

outTrainCsv = fullfile(resultDir, ...
    "PSMA_objective4_TRAIN_candidates.csv");

outValCsv = fullfile(resultDir, ...
    "PSMA_objective4_VALIDATION_candidates.csv");

outTestCsv = fullfile(resultDir, ...
    "PSMA_objective4_TEST_candidates.csv");

outSummaryXlsx = fullfile(resultDir, ...
    "PSMA_objective4_patient_based_split_summary.xlsx");

%% ============================================================
%  SETTINGS
%  ============================================================

rng(2026);  % reproducible split

trainRatio = 0.70;
valRatio   = 0.15;
testRatio  = 0.15;

fprintf("Reading Objective 4 candidate table...\n");
T = readtable(inputFile);

%% Convert variables
T.SubjectID = string(T.SubjectID);
T.case_id = string(T.case_id);
T.candidate_label = double(T.candidate_label);

fprintf("Total candidates: %d\n", height(T));
fprintf("True candidates: %d\n", sum(T.candidate_label == 1));
fprintf("False candidates: %d\n", sum(T.candidate_label == 0));

%% ============================================================
%  CREATE SUBJECT-LEVEL SUMMARY
%  ============================================================

subjects = unique(T.SubjectID, 'stable');

subjectRows = {};

for i = 1:numel(subjects)

    sid = subjects(i);
    idx = T.SubjectID == sid;

    nCandidates = sum(idx);
    nTrue = sum(T.candidate_label(idx) == 1);
    nFalse = sum(T.candidate_label(idx) == 0);
    nCases = numel(unique(T.case_id(idx)));

    hasTrue = nTrue > 0;

    subjectRows(end+1,:) = {
        sid, ...
        nCases, ...
        nCandidates, ...
        nTrue, ...
        nFalse, ...
        hasTrue
        };

end

Tsub = cell2table(subjectRows, ...
    'VariableNames', {
    'SubjectID'
    'n_cases'
    'n_candidates'
    'n_true'
    'n_false'
    'has_true_candidate'
    });

fprintf("Total unique patients / SubjectIDs: %d\n", height(Tsub));

%% ============================================================
%  STRATIFIED PATIENT-LEVEL SPLIT
%
%  We stratify by has_true_candidate:
%    group 1: patients with at least one true tumor candidate
%    group 2: patients with only false hot candidates
%  ============================================================

Tsub.split_set = strings(height(Tsub),1);

groups = unique(Tsub.has_true_candidate);

for g = 1:numel(groups)

    groupValue = groups(g);
    idxGroup = find(Tsub.has_true_candidate == groupValue);

    idxGroup = idxGroup(randperm(numel(idxGroup)));

    nGroup = numel(idxGroup);

    nTrain = round(trainRatio * nGroup);
    nVal   = round(valRatio * nGroup);
    nTest  = nGroup - nTrain - nVal;

    trainIdx = idxGroup(1:nTrain);
    valIdx   = idxGroup(nTrain+1:nTrain+nVal);
    testIdx  = idxGroup(nTrain+nVal+1:end);

    Tsub.split_set(trainIdx) = "train";
    Tsub.split_set(valIdx)   = "validation";
    Tsub.split_set(testIdx)  = "test";

    fprintf("\nStratum has_true_candidate = %d\n", groupValue);
    fprintf("Patients: %d | train %d | validation %d | test %d\n", ...
        nGroup, nTrain, nVal, nTest);

end

%% ============================================================
%  APPLY SPLIT TO CANDIDATE TABLE
%  ============================================================

T.split_set = strings(height(T),1);

for i = 1:height(Tsub)

    sid = Tsub.SubjectID(i);
    splitName = Tsub.split_set(i);

    T.split_set(T.SubjectID == sid) = splitName;

end

%% Safety check
if any(T.split_set == "")
    error("Some candidate rows did not receive split labels.");
end

%% ============================================================
%  VERIFY NO PATIENT LEAKAGE
%  ============================================================

leakRows = {};

for i = 1:numel(subjects)

    sid = subjects(i);
    splitsForSubject = unique(T.split_set(T.SubjectID == sid));

    if numel(splitsForSubject) > 1
        leakRows(end+1,:) = {
            sid, ...
            strjoin(splitsForSubject, ",")
            };
    end

end

if isempty(leakRows)
    Tleak = table(strings(0,1), strings(0,1), ...
        'VariableNames', {'SubjectID','SplitSets'});
    fprintf("\nLeakage check: PASSED. No SubjectID appears in more than one split.\n");
else
    Tleak = cell2table(leakRows, ...
        'VariableNames', {'SubjectID','SplitSets'});
    warning("Leakage detected. Check Tleak.");
end

%% ============================================================
%  CREATE SPLIT TABLES
%  ============================================================

Ttrain = T(T.split_set == "train", :);
Tval   = T(T.split_set == "validation", :);
Ttest  = T(T.split_set == "test", :);

%% ============================================================
%  SUMMARY TABLE
%  ============================================================

TsplitSummary = makeSplitSummary(T);

fprintf("\nSplit summary:\n");
disp(TsplitSummary);

%% ============================================================
%  FEATURE LIST FOR NEXT ML STEP
%  ============================================================

FeatureName = [
    "candidate_volume_mL"
    "SUVmax"
    "SUVmean"
    "SUVmedian"
    "SUVstd"
    "CTmean"
    "CTmedian"
    "CTstd"
    "centroid_x_norm"
    "centroid_y_norm"
    "centroid_z_norm"
    "distance_from_center_mm"
    "lateral_zone"
    "AP_zone"
    "craniocaudal_zone"
    "scanner"
    "radionuclide"
    "contrast"
    ];

FeatureType = [
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "numeric"
    "categorical"
    "categorical"
    "categorical"
    "categorical"
    "categorical"
    "categorical"
    ];

Tfeatures = table(FeatureName, FeatureType);

%% ============================================================
%  SAVE OUTPUTS
%  ============================================================

fprintf("\nSaving split tables...\n");

writetable(T, outAllXlsx);
writetable(T, outAllCsv);

writetable(Ttrain, outTrainXlsx);
writetable(Tval, outValXlsx);
writetable(Ttest, outTestXlsx);

writetable(Ttrain, outTrainCsv);
writetable(Tval, outValCsv);
writetable(Ttest, outTestCsv);

writetable(TsplitSummary, outSummaryXlsx, "Sheet", "Split_summary");
writetable(Tsub, outSummaryXlsx, "Sheet", "Subject_split_table");
writetable(Tleak, outSummaryXlsx, "Sheet", "Leakage_check");
writetable(Tfeatures, outSummaryXlsx, "Sheet", "Feature_list");

fprintf("\nPatient-based split complete.\n");
fprintf("All candidates with split labels:\n%s\n", outAllXlsx);
fprintf("Train table:\n%s\n", outTrainXlsx);
fprintf("Validation table:\n%s\n", outValXlsx);
fprintf("Test table:\n%s\n", outTestXlsx);
fprintf("Summary:\n%s\n", outSummaryXlsx);

%% ============================================================
%  LOCAL FUNCTION
%  ============================================================

function Tsummary = makeSplitSummary(T)

    splitNames = ["train"; "validation"; "test"];

    rows = {};

    for i = 1:numel(splitNames)

        s = splitNames(i);
        idx = T.split_set == s;

        nCandidates = sum(idx);
        nTrue = sum(T.candidate_label(idx) == 1);
        nFalse = sum(T.candidate_label(idx) == 0);

        nSubjects = numel(unique(T.SubjectID(idx)));
        nCases = numel(unique(T.case_id(idx)));

        truePercent = 100 * nTrue / nCandidates;
        falsePercent = 100 * nFalse / nCandidates;

        rows(end+1,:) = {
            s, ...
            nSubjects, ...
            nCases, ...
            nCandidates, ...
            nTrue, ...
            nFalse, ...
            truePercent, ...
            falsePercent
            };

    end

    Tsummary = cell2table(rows, ...
        'VariableNames', {
        'Split'
        'N_Subjects'
        'N_Cases'
        'N_Candidates'
        'N_True'
        'N_False'
        'TruePercent'
        'FalsePercent'
        });

end