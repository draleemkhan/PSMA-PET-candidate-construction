%% PSMA_00B_create_physical_train_val_test_folders.m
% Create physical train / validation / test folders from the internal patient split.
%
% This script copies each case's CT, PET/SUV and SEG NIfTI files into:
%   D:\PSMA-PET-CT PROJECT\PSMA_AI_Project\00_data_split_physical\train
%   D:\PSMA-PET-CT PROJECT\PSMA_AI_Project\00_data_split_physical\validation
%   D:\PSMA-PET-CT PROJECT\PSMA_AI_Project\00_data_split_physical\test
%
% Folder structure created:
%   train\imagesTr\*_0000.nii.gz  CT
%   train\imagesTr\*_0001.nii.gz  PET/SUV
%   train\labelsTr\*.nii.gz       SEG
% same for validation and test.
%
% IMPORTANT:
% - This makes real physical copies, so it may use extra disk space.
% - It does not re-split the dataset. It uses your already locked patient-level split.

clear; clc;

fprintf('\n============================================================\n');
fprintf('PSMA physical train / validation / test folder creator\n');
fprintf('============================================================\n\n');

%% ---------------- USER PATHS ----------------
projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));
rawRoot = psma_root("PSMA_INTERNAL_DATA");

splitFile = fullfile(projectRoot, ...
    "06_results", ...
    "Objective4_TRUE_vs_FALSE_candidate_dataset", ...
    "Patient_based_split", ...
    "PSMA_objective4_candidates_with_patient_split.csv");

outRoot = fullfile(projectRoot, "00_data_split_physical");

imagesDir = fullfile(rawRoot, "imagesTr");
labelsDir = fullfile(rawRoot, "labelsTr");

%% ---------------- BASIC CHECKS ----------------
fprintf('Project root: %s\n', projectRoot);
fprintf('Raw root:     %s\n', rawRoot);
fprintf('Split file:   %s\n', splitFile);
fprintf('Output root:  %s\n\n', outRoot);

assert(isfolder(projectRoot), 'Project root not found.');
assert(isfolder(rawRoot), 'Raw dataset root not found.');
assert(isfolder(imagesDir), 'imagesTr folder not found.');
assert(isfolder(labelsDir), 'labelsTr folder not found.');
assert(isfile(splitFile), 'Split file not found.');

%% ---------------- READ SPLIT FILE ----------------
T = readtable(splitFile);

fprintf('Variables in split file:\n');
disp(T.Properties.VariableNames');

% Find case ID column
if ismember('case_id', T.Properties.VariableNames)
    caseCol = 'case_id';
elseif ismember('CaseID', T.Properties.VariableNames)
    caseCol = 'CaseID';
elseif ismember('subject_id', T.Properties.VariableNames)
    caseCol = 'subject_id';
else
    error('Could not find case ID column. Expected case_id or CaseID.');
end

% Find split column
if ismember('split_set', T.Properties.VariableNames)
    splitCol = 'split_set';
elseif ismember('split', T.Properties.VariableNames)
    splitCol = 'split';
elseif ismember('Split', T.Properties.VariableNames)
    splitCol = 'Split';
else
    error('Could not find split column. Expected split_set, split, or Split.');
end

caseIDs = string(T.(caseCol));
splitSet = lower(string(T.(splitCol)));

% Standardize validation names
splitSet(splitSet == "val") = "validation";
splitSet(splitSet == "valid") = "validation";
splitSet(splitSet == "validate") = "validation";

% One row per case
U = table(caseIDs, splitSet, 'VariableNames', {'case_id','split_set'});
U = unique(U, 'rows');

% Safety: each case must have one split only
[G, caseList] = findgroups(U.case_id);
nSplitsPerCase = splitapply(@(x) numel(unique(x)), U.split_set, G);
bad = caseList(nSplitsPerCase > 1);
if ~isempty(bad)
    disp(bad(1:min(10,numel(bad))));
    error('Some case IDs occur in more than one split. Stop to avoid leakage.');
end

% Final unique case table
[caseIDsUnique, ia] = unique(U.case_id, 'stable');
SplitsUnique = U.split_set(ia);
C = table(caseIDsUnique, SplitsUnique, 'VariableNames', {'case_id','split_set'});

fprintf('\nUnique cases found: %d\n', height(C));
disp(groupsummary(C, 'split_set'));

expectedSplits = ["train", "validation", "test"];
missingSplits = setdiff(expectedSplits, unique(C.split_set));
if ~isempty(missingSplits)
    warning('Missing expected split group(s): %s', strjoin(missingSplits, ', '));
end

%% ---------------- CREATE FOLDERS ----------------
if ~isfolder(outRoot)
    mkdir(outRoot);
end

for s = expectedSplits
    mkdir_if_needed(fullfile(outRoot, s));
    mkdir_if_needed(fullfile(outRoot, s, "imagesTr"));
    mkdir_if_needed(fullfile(outRoot, s, "labelsTr"));
end

%% ---------------- COPY FILES ----------------
copyLog = table('Size',[0 8], ...
    'VariableTypes', {'string','string','string','logical','logical','logical','string','string'}, ...
    'VariableNames', {'case_id','split_set','status','ct_exists','pet_exists','seg_exists','message','destination_folder'});

nCopiedFiles = 0;
nSkippedExisting = 0;

fprintf('\nStarting copy...\n');

for i = 1:height(C)
    caseID = C.case_id(i);
    splitName = C.split_set(i);

    if ~ismember(splitName, expectedSplits)
        msg = "Unknown split name; skipped.";
        copyLog = [copyLog; {caseID, splitName, "skipped", false, false, false, msg, ""}]; %#ok<AGROW>
        continue;
    end

    srcCT  = fullfile(imagesDir, caseID + "_0000.nii.gz");
    srcPET = fullfile(imagesDir, caseID + "_0001.nii.gz");
    srcSEG = fullfile(labelsDir, caseID + ".nii.gz");

    ctOK = isfile(srcCT);
    petOK = isfile(srcPET);
    segOK = isfile(srcSEG);

    destImages = fullfile(outRoot, splitName, "imagesTr");
    destLabels = fullfile(outRoot, splitName, "labelsTr");

    dstCT  = fullfile(destImages, caseID + "_0000.nii.gz");
    dstPET = fullfile(destImages, caseID + "_0001.nii.gz");
    dstSEG = fullfile(destLabels, caseID + ".nii.gz");

    if ~(ctOK && petOK && segOK)
        msg = sprintf('Missing files: CT=%d PET=%d SEG=%d', ctOK, petOK, segOK);
        fprintf('MISSING: %s | %s\n', caseID, msg);
        copyLog = [copyLog; {caseID, splitName, "missing", ctOK, petOK, segOK, string(msg), fullfile(outRoot, splitName)}]; %#ok<AGROW>
        continue;
    end

    try
        [c1, s1] = copy_if_needed(srcCT, dstCT);
        [c2, s2] = copy_if_needed(srcPET, dstPET);
        [c3, s3] = copy_if_needed(srcSEG, dstSEG);

        nCopiedFiles = nCopiedFiles + c1 + c2 + c3;
        nSkippedExisting = nSkippedExisting + s1 + s2 + s3;

        copyLog = [copyLog; {caseID, splitName, "ok", true, true, true, "", fullfile(outRoot, splitName)}]; %#ok<AGROW>

        if mod(i,25) == 0 || i == height(C)
            fprintf('Processed %d / %d cases...\n', i, height(C));
        end
    catch ME
        fprintf('ERROR copying %s: %s\n', caseID, ME.message);
        copyLog = [copyLog; {caseID, splitName, "copy_error", true, true, true, string(ME.message), fullfile(outRoot, splitName)}]; %#ok<AGROW>
    end
end

%% ---------------- VERIFY OUTPUT COUNTS ----------------
fprintf('\nVerifying output folders...\n');

summaryRows = table('Size',[0 5], ...
    'VariableTypes', {'string','double','double','double','double'}, ...
    'VariableNames', {'split_set','n_cases_expected','n_CT','n_PET','n_SEG'});

for s = expectedSplits
    nExpected = sum(C.split_set == s);
    dImg = fullfile(outRoot, s, "imagesTr");
    dLab = fullfile(outRoot, s, "labelsTr");

    nCT  = numel(dir(fullfile(dImg, "*_0000.nii.gz")));
    nPET = numel(dir(fullfile(dImg, "*_0001.nii.gz")));
    nSEG = numel(dir(fullfile(dLab, "*.nii.gz")));

    summaryRows = [summaryRows; {s, nExpected, nCT, nPET, nSEG}]; %#ok<AGROW>
end

disp(summaryRows);

%% ---------------- SAVE LOGS ----------------
logDir = fullfile(outRoot, "logs");
mkdir_if_needed(logDir);

summaryFileCSV = fullfile(logDir, "PSMA_physical_split_summary.csv");
summaryFileXLSX = fullfile(logDir, "PSMA_physical_split_summary.xlsx");
copyLogCSV = fullfile(logDir, "PSMA_physical_split_copy_log.csv");
copyLogXLSX = fullfile(logDir, "PSMA_physical_split_copy_log.xlsx");
caseListCSV = fullfile(logDir, "PSMA_physical_split_case_list.csv");
caseListXLSX = fullfile(logDir, "PSMA_physical_split_case_list.xlsx");

writetable(summaryRows, summaryFileCSV);
writetable(summaryRows, summaryFileXLSX);
writetable(copyLog, copyLogCSV);
writetable(copyLog, copyLogXLSX);
writetable(C, caseListCSV);
writetable(C, caseListXLSX);

fprintf('\n============================================================\n');
fprintf('DONE. Physical cohort folders created.\n');
fprintf('Output root: %s\n', outRoot);
fprintf('Copied files: %d\n', nCopiedFiles);
fprintf('Skipped existing files: %d\n', nSkippedExisting);
fprintf('Summary: %s\n', summaryFileXLSX);
fprintf('Copy log: %s\n', copyLogXLSX);
fprintf('Case list: %s\n', caseListXLSX);
fprintf('============================================================\n\n');

%% ---------------- LOCAL FUNCTIONS ----------------
function mkdir_if_needed(p)
    if ~isfolder(p)
        mkdir(p);
    end
end

function [nCopied, nSkipped] = copy_if_needed(src, dst)
    nCopied = 0;
    nSkipped = 0;

    if isfile(dst)
        srcInfo = dir(src);
        dstInfo = dir(dst);
        if ~isempty(srcInfo) && ~isempty(dstInfo) && srcInfo.bytes == dstInfo.bytes
            nSkipped = 1;
            return;
        end
    end

    [ok, msg] = copyfile(src, dst, 'f');
    if ~ok
        error('copyfile failed: %s', msg);
    end
    nCopied = 1;
end
