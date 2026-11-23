function PSMA_objective9A_generate_PET_only_candidates()
% PSMA Objective 9A
% PET-ONLY candidate generation for:
%   1) Existing internal train / validation / test splits
%   2) Locked 95-case DEEP-PSMA external cohort
%
% SCIENTIFIC RULES (IDENTICAL IN EVERY COHORT)
%   PET threshold      : SUV >= 3
%   Connectivity       : 26
%   Minimum volume     : >= 1 mL
%   Maximum volume     : none
%   Ground-truth use   : NONE in Objective 9A
%
% IMPORTANT
%   This script never reads or uses a ground-truth tumor mask.
%   Ground-truth overlap and supervised labels belong to Objective 9B.
%
% OUTPUTS
%   - one PET-derived multi-label candidate map per case
%   - candidate-level checkpoint CSV
%   - case-status checkpoint CSV
%   - final candidate index CSV/XLSX
%   - final case summary CSV/XLSX
%   - text report and MAT workspace
%
% RESUMABILITY
%   Results are checkpointed after every case. Successfully completed cases
%   are skipped on rerun unless overwriteSuccessfulCases is set to true.

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 9A: PET-only candidate generation\n');
fprintf('No ground-truth masks are read or used\n');
fprintf('============================================================\n\n');

%% ========================================================================
% USER SETTINGS
% ========================================================================

% Existing input project containing the already locked physical split.
sourceProjectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));

% Brand-new Version-2 project root.
newProjectRoot = ...
    string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));

internalSplitRoot = fullfile(sourceProjectRoot,"00_data_split_physical");
externalDatasetRoot = psma_root("PSMA_EXTERNAL_DATA");

% Exact locked external case IDs are train_0001 ... train_0095.
externalCaseIDs = "train_" + compose("%04d",(1:95)');

% Candidate-generation rules.
SUV_THRESHOLD = 3.0;
MIN_VOLUME_ML = 1.0;
CONNECTIVITY = 26;

% Processing choices.
processInternal = true;
processExternal = true;
overwriteSuccessfulCases = false;

% Checkpointing is performed after each case.
saveExcelOutputs = true;

%% ========================================================================
% CREATE NEW PROJECT STRUCTURE
% ========================================================================

codeRoot = fullfile(newProjectRoot,"07_code");
resultsRoot = fullfile(newProjectRoot,"06_results", ...
    "Objective9_PET_Derived_Candidate_Pipeline", ...
    "Objective9A_PET_Only_Candidate_Generation");

if ~isfolder(newProjectRoot), mkdir(newProjectRoot); end
if ~isfolder(codeRoot), mkdir(codeRoot); end
if ~isfolder(resultsRoot), mkdir(resultsRoot); end

runStamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

fprintf('Source project root:\n%s\n\n',sourceProjectRoot);
fprintf('New Version-2 project root:\n%s\n\n',newProjectRoot);
fprintf('Objective 9A output root:\n%s\n\n',resultsRoot);

fprintf('Locked candidate rules:\n');
fprintf('  SUV threshold: %.2f\n',SUV_THRESHOLD);
fprintf('  Minimum volume: %.2f mL\n',MIN_VOLUME_ML);
fprintf('  Connectivity: %d\n',CONNECTIVITY);
fprintf('  Maximum-volume exclusion: NONE\n');
fprintf('  Ground-truth use: NONE\n\n');

%% ========================================================================
% BUILD DATASET CONFIGURATIONS
% ========================================================================

configs = {};

if processInternal
    internalSplits = ["train","validation","test"];

    for s = 1:numel(internalSplits)
        splitName = internalSplits(s);

        cfg = struct();
        cfg.dataset_name = "PSMA-PET-CT-Lesions";
        cfg.cohort = "internal";
        cfg.split = splitName;
        cfg.mode = "internal_split";
        cfg.root = fullfile(internalSplitRoot,splitName);
        cfg.case_ids = strings(0,1);

        configs{end+1} = cfg; %#ok<AGROW>
    end
end

if processExternal
    cfg = struct();
    cfg.dataset_name = "DEEP-PSMA";
    cfg.cohort = "external";
    cfg.split = "external";
    cfg.mode = "external_locked95";
    cfg.root = externalDatasetRoot;
    cfg.case_ids = externalCaseIDs;

    configs{end+1} = cfg; %#ok<AGROW>
end

%% ========================================================================
% PROCESS EACH CONFIGURATION
% ========================================================================

AllCandidates = table();
AllCaseStatus = table();
AllSummaries = table();

for d = 1:numel(configs)
    cfg = configs{d};

    fprintf('\n============================================================\n');
    fprintf('DATASET: %s | COHORT: %s | SPLIT: %s\n', ...
        cfg.dataset_name,cfg.cohort,cfg.split);
    fprintf('============================================================\n');

    outDir = fullfile(resultsRoot,cfg.cohort,cfg.split);
    labelmapDir = fullfile(outDir,"candidate_label_maps");
    checkpointDir = fullfile(outDir,"checkpoints");
    logDir = fullfile(outDir,"logs");

    if ~isfolder(outDir), mkdir(outDir); end
    if ~isfolder(labelmapDir), mkdir(labelmapDir); end
    if ~isfolder(checkpointDir), mkdir(checkpointDir); end
    if ~isfolder(logDir), mkdir(logDir); end

    candidateCheckpoint = fullfile(checkpointDir, ...
        "PSMA_objective9A_candidate_index_checkpoint.csv");

    statusCheckpoint = fullfile(checkpointDir, ...
        "PSMA_objective9A_case_status_checkpoint.csv");

    logFile = fullfile(logDir, ...
        "PSMA_objective9A_" + cfg.cohort + "_" + cfg.split + ".log");

    diary(logFile);

    caseManifest = build_case_manifest(cfg);

    fprintf('Cases discovered/locked: %d\n',height(caseManifest));

    validate_manifest(caseManifest,cfg);

    if isfile(candidateCheckpoint)
        CandidateIndex = readtable(candidateCheckpoint, ...
            'TextType','string', ...
            'VariableNamingRule','preserve');
    else
        CandidateIndex = empty_candidate_table();
    end

    if isfile(statusCheckpoint)
        CaseStatus = readtable(statusCheckpoint, ...
            'TextType','string', ...
            'VariableNamingRule','preserve');
    else
        CaseStatus = empty_status_table();
    end

    for i = 1:height(caseManifest)
        caseID = string(caseManifest.case_id(i));
        petPath = string(caseManifest.PET_path(i));

        fprintf('\n[%d/%d] %s\n',i,height(caseManifest),caseID);

        outputLabelmap = fullfile(labelmapDir, ...
            safe_filename(caseID) + "_Objective9A_PET_candidate_labelmap.nii.gz");

        alreadySuccessful = false;

        if ~overwriteSuccessfulCases && ~isempty(CaseStatus)
            row = find(CaseStatus.case_id==caseID & ...
                CaseStatus.status=="success",1,'last');

            hasCandidateRows = any(CandidateIndex.case_id==caseID);

            if ~isempty(row) && isfile(outputLabelmap) && ...
                    (hasCandidateRows || CaseStatus.retained_candidates(row)==0)
                alreadySuccessful = true;
            end
        end

        if alreadySuccessful
            fprintf('  Already completed successfully. Skipping.\n');
            continue;
        end

        started = datetime('now');

        try
            if ~isfile(petPath)
                error('PET file not found: %s',petPath);
            end

            petInfo = niftiinfo(char(petPath));
            PET = double(niftiread(char(petPath)));

            if ndims(PET)~=3
                error('Expected a 3-D PET image; found size %s.', ...
                    mat2str(size(PET)));
            end

            voxelSize_mm = double(petInfo.PixelDimensions(1:3));
            voxelVolume_mL = prod(voxelSize_mm)/1000;
            minimumVoxels = max(1,ceil(MIN_VOLUME_ML/voxelVolume_mL));

            validPET = isfinite(PET);
            hotMask = validPET & PET>=SUV_THRESHOLD;

            rawCC = bwconncomp(hotMask,CONNECTIVITY);
            rawComponents = rawCC.NumObjects;

            % Remove only components below the fixed minimum volume.
            hotMaskFiltered = bwareaopen( ...
                hotMask,minimumVoxels,CONNECTIVITY);

            retainedCC = bwconncomp( ...
                hotMaskFiltered,CONNECTIVITY);

            retainedComponents = retainedCC.NumObjects;
            excludedBelowMinimum = rawComponents-retainedComponents;

            if retainedComponents>double(intmax('uint16'))
                error(['Case contains %d candidates, exceeding uint16 ' ...
                    'label capacity.'],retainedComponents);
            end

            candidateLabelMap = zeros(size(PET),'uint16');
            caseRows = cell(retainedComponents,1);

            for k = 1:retainedComponents
                idx = retainedCC.PixelIdxList{k};
                candidateLabelMap(idx) = uint16(k);

                nVox = numel(idx);
                vol_mL = nVox*voxelVolume_mL;

                petVals = PET(idx);
                [x,y,z] = ind2sub(size(PET),idx);

                candidateGlobalID = ...
                    cfg.cohort + "|" + cfg.split + "|" + ...
                    caseID + "|" + compose("%05d",k);

                caseRows{k} = table( ...
                    candidateGlobalID, ...
                    string(cfg.dataset_name), ...
                    string(cfg.cohort), ...
                    string(cfg.split), ...
                    caseID, ...
                    double(k), ...
                    double(nVox), ...
                    double(vol_mL), ...
                    double(max(petVals,[],'omitnan')), ...
                    double(mean(petVals,'omitnan')), ...
                    double(median(petVals,'omitnan')), ...
                    double(std(petVals,'omitnan')), ...
                    double(mean(x)), ...
                    double(mean(y)), ...
                    double(mean(z)), ...
                    double(mean(x)/size(PET,1)), ...
                    double(mean(y)/size(PET,2)), ...
                    double(mean(z)/size(PET,3)), ...
                    double((max(x)-min(x)+1)*voxelSize_mm(1)), ...
                    double((max(y)-min(y)+1)*voxelSize_mm(2)), ...
                    double((max(z)-min(z)+1)*voxelSize_mm(3)), ...
                    string(petPath), ...
                    string(outputLabelmap), ...
                    SUV_THRESHOLD, ...
                    MIN_VOLUME_ML, ...
                    CONNECTIVITY, ...
                    voxelVolume_mL, ...
                    'VariableNames',{ ...
                    'candidate_global_id','dataset_name','cohort','split', ...
                    'case_id','candidate_label_value','candidate_voxels', ...
                    'candidate_volume_mL','SUVmax_QC','SUVmean_QC', ...
                    'SUVmedian_QC','SUVstd_QC', ...
                    'centroid_x_voxel','centroid_y_voxel', ...
                    'centroid_z_voxel','centroid_x_norm', ...
                    'centroid_y_norm','centroid_z_norm', ...
                    'bbox_x_mm','bbox_y_mm','bbox_z_mm', ...
                    'PET_path','candidate_labelmap_path', ...
                    'SUV_threshold','minimum_volume_mL', ...
                    'connectivity','voxel_volume_mL'});
            end

            if retainedComponents>0
                caseTable = vertcat(caseRows{:});
            else
                caseTable = empty_candidate_table();
            end

            write_labelmap(candidateLabelMap,outputLabelmap,petInfo);

            % Replace any prior rows for this case.
            if ~isempty(CandidateIndex)
                CandidateIndex(CandidateIndex.case_id==caseID,:) = [];
            end
            CandidateIndex = [CandidateIndex; caseTable]; %#ok<AGROW>

            elapsedSeconds = seconds(datetime('now')-started);

            statusRow = table( ...
                string(cfg.dataset_name),string(cfg.cohort), ...
                string(cfg.split),caseID,string(petPath), ...
                string(outputLabelmap),"success","", ...
                rawComponents,excludedBelowMinimum,retainedComponents, ...
                max_or_nan(caseTable.candidate_volume_mL), ...
                median_or_nan(caseTable.candidate_volume_mL), ...
                voxelVolume_mL,elapsedSeconds, ...
                string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
                'VariableNames',{ ...
                'dataset_name','cohort','split','case_id','PET_path', ...
                'candidate_labelmap_path','status','error_message', ...
                'raw_hot_components','excluded_below_1mL', ...
                'retained_candidates','maximum_candidate_volume_mL', ...
                'median_candidate_volume_mL','voxel_volume_mL', ...
                'elapsed_seconds','completed_at'});

            CaseStatus = replace_status_row(CaseStatus,statusRow);

            writetable(CandidateIndex,candidateCheckpoint);
            writetable(CaseStatus,statusCheckpoint);

            fprintf('  Raw hot components: %d\n',rawComponents);
            fprintf('  Excluded below 1 mL: %d\n',excludedBelowMinimum);
            fprintf('  Retained PET candidates: %d\n',retainedComponents);
            fprintf('  Label map: %s\n',outputLabelmap);

        catch ME
            elapsedSeconds = seconds(datetime('now')-started);

            statusRow = table( ...
                string(cfg.dataset_name),string(cfg.cohort), ...
                string(cfg.split),caseID,string(petPath), ...
                string(outputLabelmap),"failed",string(ME.message), ...
                NaN,NaN,NaN,NaN,NaN,NaN,elapsedSeconds, ...
                string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
                'VariableNames',{ ...
                'dataset_name','cohort','split','case_id','PET_path', ...
                'candidate_labelmap_path','status','error_message', ...
                'raw_hot_components','excluded_below_1mL', ...
                'retained_candidates','maximum_candidate_volume_mL', ...
                'median_candidate_volume_mL','voxel_volume_mL', ...
                'elapsed_seconds','completed_at'});

            CaseStatus = replace_status_row(CaseStatus,statusRow);
            writetable(CaseStatus,statusCheckpoint);

            fprintf(2,'  FAILED: %s\n',ME.message);
        end
    end

    % Final outputs for this split.
    finalCandidateCsv = fullfile(outDir, ...
        "PSMA_objective9A_PET_only_candidate_index_" + runStamp + ".csv");

    finalCandidateXlsx = fullfile(outDir, ...
        "PSMA_objective9A_PET_only_candidate_index_" + runStamp + ".xlsx");

    finalStatusCsv = fullfile(outDir, ...
        "PSMA_objective9A_case_status_" + runStamp + ".csv");

    finalStatusXlsx = fullfile(outDir, ...
        "PSMA_objective9A_case_status_" + runStamp + ".xlsx");

    writetable(CandidateIndex,finalCandidateCsv);
    writetable(CaseStatus,finalStatusCsv);

    if saveExcelOutputs
        writetable(CandidateIndex,finalCandidateXlsx);
        writetable(CaseStatus,finalStatusXlsx);
    end

    Summary = make_summary(cfg,caseManifest,CandidateIndex,CaseStatus, ...
        SUV_THRESHOLD,MIN_VOLUME_ML,CONNECTIVITY);

    summaryCsv = fullfile(outDir, ...
        "PSMA_objective9A_summary_" + runStamp + ".csv");

    summaryXlsx = fullfile(outDir, ...
        "PSMA_objective9A_summary_" + runStamp + ".xlsx");

    writetable(Summary,summaryCsv);

    if saveExcelOutputs
        writetable(Summary,summaryXlsx);
    end

    matPath = fullfile(outDir, ...
        "PSMA_objective9A_workspace_" + runStamp + ".mat");

    save(matPath,'CandidateIndex','CaseStatus','Summary', ...
        'caseManifest','cfg','SUV_THRESHOLD','MIN_VOLUME_ML', ...
        'CONNECTIVITY','-v7.3');

    reportPath = fullfile(outDir, ...
        "PSMA_objective9A_report_" + runStamp + ".txt");

    write_report(reportPath,Summary);

    disp(Summary);

    fprintf('\nFinal candidate CSV:\n%s\n',finalCandidateCsv);
    fprintf('Final case-status CSV:\n%s\n',finalStatusCsv);
    fprintf('Summary:\n%s\n',summaryCsv);
    fprintf('Workspace:\n%s\n',matPath);
    fprintf('Report:\n%s\n',reportPath);

    AllCandidates = [AllCandidates; CandidateIndex]; %#ok<AGROW>
    AllCaseStatus = [AllCaseStatus; CaseStatus]; %#ok<AGROW>
    AllSummaries = [AllSummaries; Summary]; %#ok<AGROW>

    diary off;
end

%% ========================================================================
% OVERALL COMBINED OUTPUT
% ========================================================================

combinedDir = fullfile(resultsRoot,"combined");
if ~isfolder(combinedDir), mkdir(combinedDir); end

combinedCandidateCsv = fullfile(combinedDir, ...
    "PSMA_objective9A_all_cohorts_candidate_index_" + runStamp + ".csv");

combinedStatusCsv = fullfile(combinedDir, ...
    "PSMA_objective9A_all_cohorts_case_status_" + runStamp + ".csv");

combinedSummaryCsv = fullfile(combinedDir, ...
    "PSMA_objective9A_all_cohorts_summary_" + runStamp + ".csv");

writetable(AllCandidates,combinedCandidateCsv);
writetable(AllCaseStatus,combinedStatusCsv);
writetable(AllSummaries,combinedSummaryCsv);

save(fullfile(combinedDir, ...
    "PSMA_objective9A_all_cohorts_workspace_" + runStamp + ".mat"), ...
    'AllCandidates','AllCaseStatus','AllSummaries', ...
    'SUV_THRESHOLD','MIN_VOLUME_ML','CONNECTIVITY','-v7.3');

fprintf('\n============================================================\n');
fprintf('OBJECTIVE 9A COMPLETED\n');
fprintf('============================================================\n');
fprintf('Combined candidate table:\n%s\n',combinedCandidateCsv);
fprintf('Combined status table:\n%s\n',combinedStatusCsv);
fprintf('Combined summary:\n%s\n\n',combinedSummaryCsv);

end

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function M = build_case_manifest(cfg)

switch cfg.mode
    case "internal_split"
        imagesDir = fullfile(cfg.root,"imagesTr");

        if ~isfolder(imagesDir)
            error('Internal imagesTr folder not found: %s',imagesDir);
        end

        petFiles = dir(fullfile(imagesDir,"*_0001.nii.gz"));

        if isempty(petFiles)
            petFiles = dir(fullfile(imagesDir,"*_0001.nii"));
        end

        n = numel(petFiles);
        case_id = strings(n,1);
        PET_path = strings(n,1);

        for i = 1:n
            PET_path(i) = string(fullfile( ...
                petFiles(i).folder,petFiles(i).name));

            case_id(i) = regexprep( ...
                string(petFiles(i).name), ...
                "_0001\.nii(\.gz)?$","");
        end

        M = table(case_id,PET_path);
        M = sortrows(M,'case_id');

    case "external_locked95"
        n = numel(cfg.case_ids);
        case_id = string(cfg.case_ids(:));
        PET_path = strings(n,1);

        for i = 1:n
            caseFolder = fullfile(cfg.root,case_id(i));

            if ~isfolder(caseFolder)
                error('External case folder not found: %s',caseFolder);
            end

            PET_path(i) = resolve_external_pet(caseFolder);
        end

        M = table(case_id,PET_path);
        M = sortrows(M,'case_id');

    otherwise
        error('Unknown manifest mode: %s',cfg.mode);
end

end

function petPath = resolve_external_pet(caseFolder)
% Find directories containing the exact trio PET.nii.gz, CT.nii.gz, TTB.nii.gz.
% Only PET is used in Objective 9A. Requiring the trio helps avoid selecting
% unrelated derivative PET copies.

petFiles = dir(fullfile(caseFolder,"**","PET.nii.gz"));

if isempty(petFiles)
    petFiles = dir(fullfile(caseFolder,"**","PET.nii"));
end

if isempty(petFiles)
    error('No exact PET.nii(.gz) file found under %s.',caseFolder);
end

validPaths = strings(0,1);

for i = 1:numel(petFiles)
    folder = string(petFiles(i).folder);

    hasCT = isfile(fullfile(folder,"CT.nii.gz")) || ...
        isfile(fullfile(folder,"CT.nii"));

    hasTTB = isfile(fullfile(folder,"TTB.nii.gz")) || ...
        isfile(fullfile(folder,"TTB.nii"));

    if hasCT && hasTTB
        validPaths(end+1,1) = string(fullfile( ...
            petFiles(i).folder,petFiles(i).name)); %#ok<AGROW>
    end
end

if isempty(validPaths)
    error(['PET copies were found, but no folder contained the exact ' ...
        'PET/CT/TTB trio under %s.'],caseFolder);
end

if numel(validPaths)>1
    % Deterministic choice: shortest full path, then alphabetical order.
    pathLengths = strlength(validPaths);
    minLength = min(pathLengths);
    shortest = sort(validPaths(pathLengths==minLength));
    petPath = shortest(1);

    warning(['Multiple valid PET/CT/TTB trios found for %s. ' ...
        'Using deterministic shortest path:\n%s'],caseFolder,petPath);
else
    petPath = validPaths(1);
end

end

function validate_manifest(M,cfg)

if isempty(M)
    error('No cases found for %s / %s.',cfg.cohort,cfg.split);
end

if numel(unique(M.case_id))~=height(M)
    error('Duplicate case IDs detected in %s / %s.', ...
        cfg.cohort,cfg.split);
end

missingPET = ~isfile(M.PET_path);

if any(missingPET)
    error('Missing PET files in %s / %s: %d.', ...
        cfg.cohort,cfg.split,sum(missingPET));
end

end

function T = empty_candidate_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1), ...
    strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{ ...
    'candidate_global_id','dataset_name','cohort','split', ...
    'case_id','candidate_label_value','candidate_voxels', ...
    'candidate_volume_mL','SUVmax_QC','SUVmean_QC', ...
    'SUVmedian_QC','SUVstd_QC', ...
    'centroid_x_voxel','centroid_y_voxel', ...
    'centroid_z_voxel','centroid_x_norm', ...
    'centroid_y_norm','centroid_z_norm', ...
    'bbox_x_mm','bbox_y_mm','bbox_z_mm', ...
    'PET_path','candidate_labelmap_path', ...
    'SUV_threshold','minimum_volume_mL', ...
    'connectivity','voxel_volume_mL'});

end

function T = empty_status_table()

T = table( ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{ ...
    'dataset_name','cohort','split','case_id','PET_path', ...
    'candidate_labelmap_path','status','error_message', ...
    'raw_hot_components','excluded_below_1mL', ...
    'retained_candidates','maximum_candidate_volume_mL', ...
    'median_candidate_volume_mL','voxel_volume_mL', ...
    'elapsed_seconds','completed_at'});

end

function T = replace_status_row(T,newRow)

if ~isempty(T)
    same = T.cohort==newRow.cohort & ...
        T.split==newRow.split & ...
        T.case_id==newRow.case_id;

    T(same,:) = [];
end

T = [T; newRow];

end

function write_labelmap(labelMap,outPath,referenceInfo)

outInfo = referenceInfo;
outInfo.Datatype = 'uint16';
outInfo.BitsPerPixel = 16;

if isfield(outInfo,'Filename'), outInfo.Filename = ''; end
if isfield(outInfo,'Filemoddate'), outInfo.Filemoddate = ''; end
if isfield(outInfo,'Filesize'), outInfo.Filesize = 0; end

outPath = string(outPath);

% niftiwrite expects the base path and adds the NIfTI extension.
if endsWith(outPath,".nii.gz")
    basePath = extractBefore(outPath,strlength(outPath)-6);
    niftiwrite(labelMap,char(basePath),outInfo,'Compressed',true);
elseif endsWith(outPath,".nii")
    basePath = extractBefore(outPath,strlength(outPath)-3);
    niftiwrite(labelMap,char(basePath),outInfo,'Compressed',false);
else
    niftiwrite(labelMap,char(outPath),outInfo,'Compressed',true);
end

end

function S = make_summary(cfg,M,CandidateIndex,CaseStatus, ...
    SUVthr,minVolume,connectivity)

successRows = CaseStatus.status=="success";
failedRows = CaseStatus.status=="failed";

SummaryMetric = [
    "Dataset"
    "Cohort"
    "Split"
    "Cases expected"
    "Cases passed"
    "Cases failed"
    "Raw SUV-hot components"
    "Components excluded below 1 mL"
    "Retained PET-derived candidates"
    "Maximum retained candidate volume mL"
    "Median retained candidate volume mL"
    "SUV threshold"
    "Minimum retained volume mL"
    "Maximum-volume exclusion"
    "Connectivity"
    "Ground-truth masks used"
];

SummaryValue = [
    string(cfg.dataset_name)
    string(cfg.cohort)
    string(cfg.split)
    string(height(M))
    string(sum(successRows))
    string(sum(failedRows))
    string(sum(CaseStatus.raw_hot_components(successRows),'omitnan'))
    string(sum(CaseStatus.excluded_below_1mL(successRows),'omitnan'))
    string(height(CandidateIndex))
    string(max_or_nan(CandidateIndex.candidate_volume_mL))
    string(median_or_nan(CandidateIndex.candidate_volume_mL))
    string(SUVthr)
    string(minVolume)
    "NONE"
    string(connectivity)
    "NO"
];

S = table(SummaryMetric,SummaryValue);

end

function write_report(path,Summary)

fid = fopen(path,'w');

if fid<0
    warning('Could not create report: %s',path);
    return;
end

fprintf(fid,'PSMA Objective 9A: PET-only candidate generation\n');
fprintf(fid,'Generated: %s\n\n',string(datetime('now')));

for i = 1:height(Summary)
    fprintf(fid,'%-50s %s\n', ...
        Summary.SummaryMetric(i),Summary.SummaryValue(i));
end

fclose(fid);

end

function value = max_or_nan(x)

if isempty(x)
    value = NaN;
else
    value = max(double(x),[],'omitnan');
end

end

function value = median_or_nan(x)

if isempty(x)
    value = NaN;
else
    value = median(double(x),'omitnan');
end

end

function name = safe_filename(value)

name = string(regexprep(char(value),'[^a-zA-Z0-9_\-]','_'));

end
