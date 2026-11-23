%% PSMA Objective 9F - Build Internal MASTER Candidate-Level Dataset
% Merges:
%   Objective 9B candidate labels
%   Objective 9C routine features
%   Objective 9D PET radiomics
%   Objective 9E CT radiomics
%
% Primary key:
%   candidate_global_id
%
% Expected internal total:
%   17,404 candidates
%
% Run this script only after PET and CT radiomics extraction has completed.
%
% Author: PSMA PET/CT Project
% Date: 2026-07-17

clear;
clc;

%% ------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------

projectRoot = ...
    string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));

resultsRoot = fullfile(projectRoot, "06_results", ...
    "Objective9_PET_Derived_Candidate_Pipeline");

outputDir = fullfile(resultsRoot, ...
    "Objective9F_MASTER_Candidate_Dataset", "internal");

expectedRows = 17404;
keyVar = "candidate_global_id";

% Set true only when you intentionally want to overwrite existing outputs.
overwriteExisting = true;

if ~isfolder(outputDir)
    mkdir(outputDir);
end

masterCsv = fullfile(outputDir, ...
    "PSMA_objective9F_internal_MASTER_candidate_dataset.csv");

masterMat = fullfile(outputDir, ...
    "PSMA_objective9F_internal_MASTER_candidate_dataset.mat");

qcCsv = fullfile(outputDir, ...
    "PSMA_objective9F_internal_MASTER_merge_QC.csv");

summaryTxt = fullfile(outputDir, ...
    "PSMA_objective9F_internal_MASTER_summary.txt");

if ~overwriteExisting
    existingOutputs = [isfile(masterCsv), isfile(masterMat), ...
        isfile(qcCsv), isfile(summaryTxt)];

    if any(existingOutputs)
        error("One or more Objective 9F outputs already exist. " + ...
            "Set overwriteExisting=true to replace them.");
    end
end

%% ------------------------------------------------------------------------
% LOCATE INPUT FILES
% -------------------------------------------------------------------------

fprintf("\n============================================================\n");
fprintf("OBJECTIVE 9F: BUILD INTERNAL MASTER CANDIDATE DATASET\n");
fprintf("============================================================\n");

splits = ["train", "validation", "test"];

labelsFiles = strings(numel(splits), 1);
routineFiles = strings(numel(splits), 1);

for s = 1:numel(splits)
    splitName = splits(s);

    labelsFiles(s) = findExactSplitCsv( ...
        resultsRoot, ...
        "Objective9B_GT_Overlap_and_Label_Assignment", ...
        splitName, ...
        "PSMA_objective9B_candidate_overlap_checkpoint.csv");

    routineFiles(s) = findSplitCsvByTokens( ...
        resultsRoot, ...
        "Objective9C", ...
        splitName, ...
        ["routine", "checkpoint"], ...
        ["case_summary", "failure", "failures"]);
end

petFile = findExactRadiomicsCheckpoint( ...
    resultsRoot, "PET", ...
    "PSMA_objective9PET_internal_PET_radiomics_checkpoint.csv");

ctFile = findExactRadiomicsCheckpoint( ...
    resultsRoot, "CT", ...
    "PSMA_objective9CT_internal_CT_radiomics_checkpoint.csv");

fprintf("\nInput files located:\n");
for s = 1:numel(splits)
    fprintf("Labels %-10s: %s\n", splits(s), labelsFiles(s));
end
for s = 1:numel(splits)
    fprintf("Routine %-9s: %s\n", splits(s), routineFiles(s));
end
fprintf("PET radiomics : %s\n", petFile);
fprintf("CT radiomics  : %s\n", ctFile);

%% ------------------------------------------------------------------------
% READ AND CONCATENATE TABLES
% -------------------------------------------------------------------------

fprintf("\nReading and concatenating split-level input tables...\n");

labels = table;
routine = table;

for s = 1:numel(splits)
    currentLabels = readCsvStrict(labelsFiles(s));

    currentRoutine = readCsvStrict(routineFiles(s));


    % Add split only when the source table does not already contain one.
    if ~hasSplitVariable(currentLabels)
        currentLabels.dataset_split = ...
            repmat(splits(s), height(currentLabels), 1);
    end

    if ~hasSplitVariable(currentRoutine)
        currentRoutine.dataset_split = ...
            repmat(splits(s), height(currentRoutine), 1);
    end

    if s == 1
        labels = currentLabels;
        routine = currentRoutine;
    else
        labels = verticallyConcatenateCompatible( ...
            labels, currentLabels, "Objective 9B labels");
        routine = verticallyConcatenateCompatible( ...
            routine, currentRoutine, "Objective 9C routine features");
    end
end

pet = readCsvStrict(petFile);

ct = readCsvStrict(ctFile);


fprintf("Labels rows : %d\n", height(labels));
fprintf("Routine rows: %d\n", height(routine));
fprintf("PET rows    : %d\n", height(pet));
fprintf("CT rows     : %d\n", height(ct));

%% ------------------------------------------------------------------------
% STANDARDIZE AND VALIDATE PRIMARY KEY
% -------------------------------------------------------------------------

labels = standardizeKey(labels, keyVar, "Labels");
routine = standardizeKey(routine, keyVar, "Routine");
pet = standardizeKey(pet, keyVar, "PET");
ct = standardizeKey(ct, keyVar, "CT");

validateUniqueKey(labels, keyVar, "Labels");
validateUniqueKey(routine, keyVar, "Routine");
validateUniqueKey(pet, keyVar, "PET");
validateUniqueKey(ct, keyVar, "CT");

%% ------------------------------------------------------------------------
% CHECK RADIOMICS COMPLETENESS BEFORE MERGE
% -------------------------------------------------------------------------

if height(pet) ~= expectedRows
    warning("PET checkpoint contains %d rows; expected %d.", ...
        height(pet), expectedRows);
end

if height(ct) ~= expectedRows
    warning("CT checkpoint contains %d rows; expected %d.", ...
        height(ct), expectedRows);
end

%% ------------------------------------------------------------------------
% PREPARE TABLES FOR SAFE JOINING
% -------------------------------------------------------------------------

% Labels table is treated as the anchor table.
% Remove duplicate non-key columns from the other tables before joining.
routine = removeOverlappingNonKeyVariables(routine, labels, keyVar);

% Add the routine variables to the reference list before PET preparation.
referenceAfterRoutine = outerjoin(labels, routine, ...
    "Keys", keyVar, ...
    "MergeKeys", true, ...
    "Type", "left");

pet = removeOverlappingNonKeyVariables(pet, referenceAfterRoutine, keyVar);
pet = prefixFeatureVariables(pet, keyVar, "PET_rad__");

referenceAfterPET = outerjoin(referenceAfterRoutine, pet, ...
    "Keys", keyVar, ...
    "MergeKeys", true, ...
    "Type", "left");

ct = removeOverlappingNonKeyVariables(ct, referenceAfterPET, keyVar);
ct = prefixFeatureVariables(ct, keyVar, "CT_rad__");

%% ------------------------------------------------------------------------
% MERGE
% -------------------------------------------------------------------------

fprintf("\nMerging tables using %s...\n", keyVar);

master = outerjoin(labels, routine, ...
    "Keys", keyVar, ...
    "MergeKeys", true, ...
    "Type", "left");

master = outerjoin(master, pet, ...
    "Keys", keyVar, ...
    "MergeKeys", true, ...
    "Type", "left");

master = outerjoin(master, ct, ...
    "Keys", keyVar, ...
    "MergeKeys", true, ...
    "Type", "left");

% Preserve the anchor-table ordering where possible.
[~, anchorOrder] = ismember(labels.(keyVar), master.(keyVar));
if all(anchorOrder > 0)
    master = master(anchorOrder, :);
end

%% ------------------------------------------------------------------------
% POST-MERGE QUALITY CONTROL
% -------------------------------------------------------------------------

validateUniqueKey(master, keyVar, "MASTER");

missingFromRoutine = setdiff(labels.(keyVar), routine.(keyVar));
missingFromPET = setdiff(labels.(keyVar), pet.(keyVar));
missingFromCT = setdiff(labels.(keyVar), ct.(keyVar));

extraInRoutine = setdiff(routine.(keyVar), labels.(keyVar));
extraInPET = setdiff(pet.(keyVar), labels.(keyVar));
extraInCT = setdiff(ct.(keyVar), labels.(keyVar));

qcMetric = [
    "Expected MASTER rows"
    "Labels rows"
    "Routine rows"
    "PET rows"
    "CT rows"
    "MASTER rows"
    "Missing labels->routine"
    "Missing labels->PET"
    "Missing labels->CT"
    "Extra routine IDs"
    "Extra PET IDs"
    "Extra CT IDs"
    "MASTER variables"
];

qcValue = [
    expectedRows
    height(labels)
    height(routine)
    height(pet)
    height(ct)
    height(master)
    numel(missingFromRoutine)
    numel(missingFromPET)
    numel(missingFromCT)
    numel(extraInRoutine)
    numel(extraInPET)
    numel(extraInCT)
    width(master)
];

qc = table(qcMetric, qcValue, ...
    'VariableNames', {'QC_metric', 'Value'});

disp(qc);

%% ------------------------------------------------------------------------
% SPLIT COUNTS
% -------------------------------------------------------------------------

splitVar = locateSplitVariable(master);

if strlength(splitVar) > 0
    splitValues = lower(strtrim(string(master.(splitVar))));
    [splitGroups, ~, splitIndex] = unique(splitValues);
    splitCounts = accumarray(splitIndex, 1);

    splitTable = table(splitGroups, splitCounts, ...
        'VariableNames', {'Split', 'Candidate_count'});

    fprintf("\nCandidate counts by split:\n");
    disp(splitTable);
else
    splitTable = table;
    warning("No train/validation/test split variable was identified.");
end

%% ------------------------------------------------------------------------
% FINAL PASS/FAIL CONDITIONS
% -------------------------------------------------------------------------

problems = strings(0, 1);

if height(master) ~= expectedRows
    problems(end + 1) = ...
        "MASTER row count is " + height(master) + ...
        "; expected " + expectedRows + ".";
end

if ~isempty(missingFromRoutine)
    problems(end + 1) = ...
        numel(missingFromRoutine) + ...
        " label IDs are missing from routine features.";
end

if ~isempty(missingFromPET)
    problems(end + 1) = ...
        numel(missingFromPET) + ...
        " label IDs are missing from PET radiomics.";
end

if ~isempty(missingFromCT)
    problems(end + 1) = ...
        numel(missingFromCT) + ...
        " label IDs are missing from CT radiomics.";
end

if ~isempty(extraInRoutine)
    problems(end + 1) = ...
        numel(extraInRoutine) + ...
        " routine-feature IDs are absent from labels.";
end

if ~isempty(extraInPET)
    problems(end + 1) = ...
        numel(extraInPET) + ...
        " PET-radiomics IDs are absent from labels.";
end

if ~isempty(extraInCT)
    problems(end + 1) = ...
        numel(extraInCT) + ...
        " CT-radiomics IDs are absent from labels.";
end

%% ------------------------------------------------------------------------
% SAVE OUTPUTS
% -------------------------------------------------------------------------

fprintf("\nSaving Objective 9F outputs...\n");

writetable(master, masterCsv);
save(masterMat, "master", "qc", "splitTable", "-v7.3");
writetable(qc, qcCsv);

fid = fopen(summaryTxt, "w");
if fid == -1
    error("Unable to create summary file: %s", summaryTxt);
end
cleanupObj = onCleanup(@() fclose(fid));

fprintf(fid, "PSMA OBJECTIVE 9F INTERNAL MASTER DATASET\n");
fprintf(fid, "Generated: %s\n\n", ...
    char(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")));

fprintf(fid, "LABEL FILES\n");
for s = 1:numel(splits)
    fprintf(fid, "%s: %s\n", splits(s), labelsFiles(s));
end
fprintf(fid, "\nROUTINE FEATURE FILES\n");
for s = 1:numel(splits)
    fprintf(fid, "%s: %s\n", splits(s), routineFiles(s));
end
fprintf(fid, "\nPET file:\n%s\n\n", petFile);
fprintf(fid, "CT file:\n%s\n\n", ctFile);

fprintf(fid, "Expected rows: %d\n", expectedRows);
fprintf(fid, "MASTER rows: %d\n", height(master));
fprintf(fid, "MASTER variables: %d\n", width(master));
fprintf(fid, "Missing from routine: %d\n", numel(missingFromRoutine));
fprintf(fid, "Missing from PET: %d\n", numel(missingFromPET));
fprintf(fid, "Missing from CT: %d\n", numel(missingFromCT));
fprintf(fid, "Extra routine IDs: %d\n", numel(extraInRoutine));
fprintf(fid, "Extra PET IDs: %d\n", numel(extraInPET));
fprintf(fid, "Extra CT IDs: %d\n\n", numel(extraInCT));

if ~isempty(splitTable)
    fprintf(fid, "SPLIT COUNTS\n");
    for i = 1:height(splitTable)
        fprintf(fid, "%s: %d\n", ...
            splitTable.Split(i), splitTable.Candidate_count(i));
    end
    fprintf(fid, "\n");
end

if isempty(problems)
    fprintf(fid, "FINAL STATUS: PASS\n");
else
    fprintf(fid, "FINAL STATUS: REVIEW REQUIRED\n");
    for i = 1:numel(problems)
        fprintf(fid, "- %s\n", problems(i));
    end
end

clear cleanupObj;

fprintf("\n============================================================\n");

if isempty(problems)
    fprintf("OBJECTIVE 9F COMPLETED SUCCESSFULLY\n");
    fprintf("MASTER rows: %d/%d\n", height(master), expectedRows);
else
    fprintf("OBJECTIVE 9F COMPLETED WITH QC WARNINGS\n");
    for i = 1:numel(problems)
        fprintf("- %s\n", problems(i));
    end
end

fprintf("MASTER CSV: %s\n", masterCsv);
fprintf("MASTER MAT: %s\n", masterMat);
fprintf("QC CSV    : %s\n", qcCsv);
fprintf("Summary   : %s\n", summaryTxt);
fprintf("============================================================\n");

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================



function tbl = readCsvStrict(filePath)
    % Read a CSV while explicitly forcing the first line to be the header.
    % This prevents MATLAB from replacing genuine names with Var1, Var2, ...
    % when some later rows contain extra trailing fields.

    filePath = string(filePath);

    fid = fopen(filePath, "r");
    if fid == -1
        error("Unable to open CSV file: %s", filePath);
    end
    cleanupObj = onCleanup(@() fclose(fid));

    headerLine = fgetl(fid);
    if ~ischar(headerLine) && ~isstring(headerLine)
        error("CSV file is empty: %s", filePath);
    end

    headerNames = split(string(headerLine), ",");
    headerNames = erase(headerNames, '"');
    headerNames = strtrim(headerNames);

    if any(strlength(headerNames) == 0)
        error("CSV header contains an empty variable name: %s", filePath);
    end

    clear cleanupObj;

    opts = detectImportOptions(filePath, ...
        "FileType", "text", ...
        "Delimiter", ",", ...
        "VariableNamingRule", "preserve");

    opts.VariableNamesLine = 1;
    opts.DataLines = [2, Inf];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Force the exact names from the first line.
    if numel(opts.VariableNames) ~= numel(headerNames)
        opts = delimitedTextImportOptions( ...
            "NumVariables", numel(headerNames), ...
            "Delimiter", ",");

        opts.VariableNames = cellstr(headerNames);
        opts.VariableNamesLine = 1;
        opts.DataLines = [2, Inf];
        opts.VariableTypes = repmat("string", 1, numel(headerNames));
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
    else
        opts.VariableNames = cellstr(headerNames);
    end

    tbl = readtable(filePath, opts);

    % Final verification: the primary header must not be imported as data.
    if any(startsWith(string(tbl.Properties.VariableNames), "Var"))
        warning("Some autogenerated variable names remain in %s.", filePath);
    end
end


function filePath = findExactSplitCsv( ...
        resultsRoot, objectiveFolderToken, splitName, exactFileName)

    listing = dir(fullfile(resultsRoot, "**", exactFileName));

    if isempty(listing)
        error("Required file not found: %s for split %s", ...
            exactFileName, splitName);
    end

    matches = strings(0, 1);

    for i = 1:numel(listing)
        candidate = string(fullfile(listing(i).folder, listing(i).name));
        candidateLower = lower(candidate);

        if contains(candidateLower, lower(objectiveFolderToken)) && ...
                contains(candidateLower, ...
                lower(filesep + "internal" + filesep + splitName + filesep))
            matches(end + 1, 1) = candidate; %#ok<AGROW>
        end
    end

    if numel(matches) ~= 1
        fprintf("\nMatches for %s / %s:\n", objectiveFolderToken, splitName);
        for i = 1:numel(matches)
            fprintf("  [%d] %s\n", i, matches(i));
        end
        error("Expected exactly one file but found %d.", numel(matches));
    end

    filePath = matches(1);
end


function filePath = findSplitCsvByTokens( ...
        resultsRoot, objectiveToken, splitName, ...
        requiredTokens, excludedTokens)

    listing = dir(fullfile(resultsRoot, "**", "*.csv"));
    matches = strings(0, 1);

    for i = 1:numel(listing)
        candidate = string(fullfile(listing(i).folder, listing(i).name));
        candidateLower = lower(candidate);

        if ~contains(candidateLower, lower(objectiveToken))
            continue;
        end

        splitPattern = lower(filesep + "internal" + filesep + ...
            splitName + filesep);

        if ~contains(candidateLower, splitPattern)
            continue;
        end

        include = true;

        for t = 1:numel(requiredTokens)
            if ~contains(candidateLower, lower(requiredTokens(t)))
                include = false;
                break;
            end
        end

        if ~include
            continue;
        end

        for t = 1:numel(excludedTokens)
            if contains(candidateLower, lower(excludedTokens(t)))
                include = false;
                break;
            end
        end

        if include
            matches(end + 1, 1) = candidate; %#ok<AGROW>
        end
    end

    if numel(matches) ~= 1
        fprintf("\nRoutine-feature matches for split %s:\n", splitName);
        for i = 1:numel(matches)
            fprintf("  [%d] %s\n", i, matches(i));
        end

        error("Expected exactly one Objective 9C routine-feature " + ...
            "checkpoint for split %s, but found %d.", ...
            splitName, numel(matches));
    end

    filePath = matches(1);
end


function filePath = findExactRadiomicsCheckpoint( ...
        resultsRoot, modality, exactFileName)

    expectedFolder = fullfile(resultsRoot, ...
        "Objective9D9E_Direct_PET_CT_Radiomics", ...
        "internal", modality);

    filePath = string(fullfile(expectedFolder, exactFileName));

    if ~isfile(filePath)
        error("Radiomics checkpoint was not found: %s", filePath);
    end
end


function tf = hasSplitVariable(tbl)

    vars = string(tbl.Properties.VariableNames);
    normalizedVars = lower(regexprep(vars, "[^a-zA-Z0-9]", ""));

    candidates = [
        "split"
        "datasetsplit"
        "cohortsplit"
        "internalsplit"
        "partition"
        "set"
    ];

    tf = false;

    for i = 1:numel(candidates)
        if any(normalizedVars == candidates(i))
            tf = true;
            return;
        end
    end
end



function tbl = removeEmptyAutoVariables(tbl, tableDescription)

    vars = string(tbl.Properties.VariableNames);
    autoMask = ~cellfun(@isempty, regexp(cellstr(vars), '^Var\d+$', 'once'));

    removeMask = false(size(vars));

    for i = find(autoMask)
        column = tbl.(vars(i));

        if isstring(column) || ischar(column) || iscellstr(column)
            values = strtrim(string(column));
            isEmptyColumn = all(ismissing(values) | strlength(values) == 0);
        elseif isnumeric(column)
            isEmptyColumn = all(isnan(column));
        elseif islogical(column)
            isEmptyColumn = false;
        elseif isdatetime(column)
            isEmptyColumn = all(isnat(column));
        elseif iscategorical(column)
            isEmptyColumn = all(isundefined(column));
        elseif iscell(column)
            isEmptyColumn = all(cellfun(@isempty, column));
        else
            try
                isEmptyColumn = all(ismissing(column));
            catch
                isEmptyColumn = false;
            end
        end

        if isEmptyColumn
            removeMask(i) = true;
        end
    end

    removedVars = vars(removeMask);

    if ~isempty(removedVars)
        fprintf("Removed %d empty autogenerated columns from %s: %s\\n", ...
            numel(removedVars), tableDescription, ...
            strjoin(removedVars, ", "));

        tbl(:, cellstr(removedVars)) = [];
    end
end


function combined = verticallyConcatenateCompatible( ...
        firstTable, secondTable, tableDescription)

    vars1 = string(firstTable.Properties.VariableNames);
    vars2 = string(secondTable.Properties.VariableNames);

    missingFromSecond = setdiff(vars1, vars2, "stable");
    extraInSecond = setdiff(vars2, vars1, "stable");

    if ~isempty(missingFromSecond) || ~isempty(extraInSecond)
        fprintf("\nVariable mismatch while combining %s.\n", ...
            tableDescription);

        if ~isempty(missingFromSecond)
            fprintf("Missing from later split: %s\n", ...
                strjoin(missingFromSecond, ", "));
        end

        if ~isempty(extraInSecond)
            fprintf("Extra in later split: %s\n", ...
                strjoin(extraInSecond, ", "));
        end

        error("Split-level tables are not structurally compatible.");
    end

    secondTable = secondTable(:, cellstr(vars1));
    combined = [firstTable; secondTable];
end


function filePath = findOneCsv(rootDir, requiredTokens)

    listing = dir(fullfile(rootDir, "**", "*.csv"));

    if isempty(listing)
        error("No CSV files were found below: %s", rootDir);
    end

    fullPaths = strings(numel(listing), 1);
    keep = true(numel(listing), 1);

    for i = 1:numel(listing)
        fullPaths(i) = string(fullfile( ...
            listing(i).folder, listing(i).name));

        candidateLower = lower(fullPaths(i));

        for t = 1:numel(requiredTokens)
            tokenLower = lower(string(requiredTokens{t}));

            if ~contains(candidateLower, tokenLower)
                keep(i) = false;
                break;
            end
        end
    end

    matches = fullPaths(keep);

    % Exclude failure, QC, summary, or previously generated MASTER files.
    exclusionTokens = [
        "failure"
        "failures"
        "merge_qc"
        "master_candidate_dataset"
    ];

    for e = 1:numel(exclusionTokens)
        matches(contains(lower(matches), exclusionTokens(e))) = [];
    end

    if isempty(matches)
        error("No CSV matched required tokens: %s", ...
            strjoin(string(requiredTokens), ", "));
    end

    if numel(matches) > 1
        fprintf("\nMultiple CSV files matched these tokens: %s\n", ...
            strjoin(string(requiredTokens), ", "));

        for i = 1:numel(matches)
            fprintf("  [%d] %s\n", i, matches(i));
        end

        error("Input-file match is ambiguous. " + ...
            "Edit the required tokens or set the path explicitly.");
    end

    filePath = char(matches(1));
end


function tbl = standardizeKey(tbl, requiredKey, tableName)

    vars = string(tbl.Properties.VariableNames);

    exactIndex = find(strcmpi(vars, requiredKey), 1);

    if isempty(exactIndex)
        normalizedVars = lower(regexprep(vars, "[^a-zA-Z0-9]", ""));
        normalizedRequired = lower(regexprep(requiredKey, ...
            "[^a-zA-Z0-9]", ""));

        exactIndex = find(normalizedVars == normalizedRequired, 1);
    end

    if isempty(exactIndex)
        error("%s table does not contain required key: %s", ...
            tableName, requiredKey);
    end

    currentName = vars(exactIndex);

    if currentName ~= requiredKey
        tbl.Properties.VariableNames{exactIndex} = char(requiredKey);
    end

    tbl.(requiredKey) = strtrim(string(tbl.(requiredKey)));

    invalidKey = ismissing(tbl.(requiredKey)) | ...
        strlength(tbl.(requiredKey)) == 0;

    if any(invalidKey)
        error("%s table contains %d missing/empty key values.", ...
            tableName, nnz(invalidKey));
    end
end


function validateUniqueKey(tbl, keyVar, tableName)

    keys = string(tbl.(keyVar));

    [uniqueKeys, ~, groupIndex] = unique(keys);
    counts = accumarray(groupIndex, 1);
    duplicateKeys = uniqueKeys(counts > 1);

    if ~isempty(duplicateKeys)
        previewCount = min(10, numel(duplicateKeys));
        preview = strjoin(duplicateKeys(1:previewCount), ", ");

        error("%s table contains %d duplicated %s values. " + ...
            "Examples: %s", ...
            tableName, numel(duplicateKeys), keyVar, preview);
    end
end


function tbl = removeOverlappingNonKeyVariables(tbl, referenceTbl, keyVar)

    vars = string(tbl.Properties.VariableNames);
    referenceVars = string(referenceTbl.Properties.VariableNames);

    overlap = intersect(vars, referenceVars, "stable");
    overlap(overlap == keyVar) = [];

    if ~isempty(overlap)
        tbl(:, cellstr(overlap)) = [];
    end
end


function tbl = prefixFeatureVariables(tbl, keyVar, prefix)

    vars = string(tbl.Properties.VariableNames);
    featureVars = vars(vars ~= keyVar);

    newNames = strings(size(featureVars));

    for i = 1:numel(featureVars)
        safeName = matlab.lang.makeValidName(prefix + featureVars(i));
        newNames(i) = safeName;
    end

    if ~isempty(featureVars)
        tbl = renamevars(tbl, featureVars, newNames);
    end
end


function splitVar = locateSplitVariable(tbl)

    vars = string(tbl.Properties.VariableNames);
    normalizedVars = lower(regexprep(vars, "[^a-zA-Z0-9]", ""));

    candidates = [
        "split"
        "datasetsplit"
        "cohortsplit"
        "internal_split"
        "partition"
        "set"
    ];

    splitVar = "";

    for i = 1:numel(candidates)
        normalizedCandidate = lower(regexprep(candidates(i), ...
            "[^a-zA-Z0-9]", ""));

        index = find(normalizedVars == normalizedCandidate, 1);

        if ~isempty(index)
            splitVar = vars(index);
            return;
        end
    end
end
