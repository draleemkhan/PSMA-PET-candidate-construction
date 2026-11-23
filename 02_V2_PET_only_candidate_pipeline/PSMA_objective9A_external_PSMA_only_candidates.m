function PSMA_objective9A_external_PSMA_only_candidates()
% PSMA Objective 9A - corrected external branch only
clc; close all;

newProjectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
externalDatasetRoot = psma_root("PSMA_EXTERNAL_DATA");

outDir = fullfile(newProjectRoot,"06_results", ...
    "Objective9_PET_Derived_Candidate_Pipeline", ...
    "Objective9A_PET_Only_Candidate_Generation","external","external");
labelmapDir = fullfile(outDir,"candidate_label_maps");
checkpointDir = fullfile(outDir,"checkpoints");
logDir = fullfile(outDir,"logs");
if ~isfolder(outDir), mkdir(outDir); end
if ~isfolder(labelmapDir), mkdir(labelmapDir); end
if ~isfolder(checkpointDir), mkdir(checkpointDir); end
if ~isfolder(logDir), mkdir(logDir); end

SUV_THRESHOLD = 3.0;
MIN_VOLUME_ML = 1.0;
CONNECTIVITY = 26;
overwriteSuccessfulCases = false;
saveExcelOutputs = true;

caseIDs = "train_" + compose("%04d",(1:95)');
runStamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

candidateCheckpoint = fullfile(checkpointDir,"PSMA_objective9A_candidate_index_checkpoint.csv");
statusCheckpoint = fullfile(checkpointDir,"PSMA_objective9A_case_status_checkpoint.csv");

diary(fullfile(logDir,"PSMA_objective9A_external_PSMA_only.log"));

PET_path = strings(numel(caseIDs),1);
for i = 1:numel(caseIDs)
    caseFolder = fullfile(externalDatasetRoot,caseIDs(i),"PSMA");
    petGz = fullfile(caseFolder,"PET.nii.gz");
    petNii = fullfile(caseFolder,"PET.nii");
    if isfile(petGz)
        PET_path(i) = string(petGz);
    elseif isfile(petNii)
        PET_path(i) = string(petNii);
    else
        error('PSMA PET not found for %s',caseIDs(i));
    end
    if ~contains(upper(PET_path(i)), upper(filesep + "PSMA" + filesep))
        error('Safety stop: external PET is not in PSMA folder: %s',PET_path(i));
    end
end
Manifest = table(caseIDs,PET_path,'VariableNames',{'case_id','PET_path'});

if isfile(candidateCheckpoint)
    CandidateIndex = readtable(candidateCheckpoint,'Delimiter',',','VariableNamingRule','preserve');
else
    CandidateIndex = empty_candidate_table();
end
if isfile(statusCheckpoint)
    CaseStatus = readtable(statusCheckpoint,'Delimiter',',','VariableNamingRule','preserve');
else
    CaseStatus = empty_status_table();
end

for i = 1:height(Manifest)
    caseID = string(Manifest.case_id(i));
    petPath = string(Manifest.PET_path(i));
    fprintf('\n[%d/%d] %s\n',i,height(Manifest),caseID);

    outputLabelmap = fullfile(labelmapDir,caseID + "_Objective9A_PET_candidate_labelmap.nii.gz");

    alreadySuccessful = false;
    if ~overwriteSuccessfulCases && ~isempty(CaseStatus)
        row = find(CaseStatus.case_id==caseID & CaseStatus.status=="success",1,'last');
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
        info = niftiinfo(char(petPath));
        PET = double(niftiread(char(petPath)));
        voxelSize_mm = double(info.PixelDimensions(1:3));
        voxelVolume_mL = prod(voxelSize_mm)/1000;
        minimumVoxels = max(1,ceil(MIN_VOLUME_ML/voxelVolume_mL));

        hotMask = isfinite(PET) & PET>=SUV_THRESHOLD;
        rawCC = bwconncomp(hotMask,CONNECTIVITY);
        filteredMask = bwareaopen(hotMask,minimumVoxels,CONNECTIVITY);
        retainedCC = bwconncomp(filteredMask,CONNECTIVITY);

        rawComponents = rawCC.NumObjects;
        retainedComponents = retainedCC.NumObjects;
        excludedBelowMinimum = rawComponents-retainedComponents;

        candidateLabelMap = zeros(size(PET),'uint16');
        caseRows = cell(retainedComponents,1);

        for k = 1:retainedComponents
            idx = retainedCC.PixelIdxList{k};
            candidateLabelMap(idx) = uint16(k);
            nVox = numel(idx);
            vol_mL = nVox*voxelVolume_mL;
            vals = PET(idx);
            [x,y,z] = ind2sub(size(PET),idx);
            globalID = "external|external|" + caseID + "|" + compose("%05d",k);

            caseRows{k} = table(globalID,"DEEP-PSMA","external","external",caseID, ...
                double(k),double(nVox),double(vol_mL), ...
                double(max(vals,[],'omitnan')),double(mean(vals,'omitnan')), ...
                double(median(vals,'omitnan')),double(std(vals,'omitnan')), ...
                double(mean(x)),double(mean(y)),double(mean(z)), ...
                double(mean(x)/size(PET,1)),double(mean(y)/size(PET,2)), ...
                double(mean(z)/size(PET,3)), ...
                double((max(x)-min(x)+1)*voxelSize_mm(1)), ...
                double((max(y)-min(y)+1)*voxelSize_mm(2)), ...
                double((max(z)-min(z)+1)*voxelSize_mm(3)), ...
                petPath,string(outputLabelmap),SUV_THRESHOLD,MIN_VOLUME_ML, ...
                CONNECTIVITY,voxelVolume_mL, ...
                'VariableNames',{'candidate_global_id','dataset_name','cohort','split', ...
                'case_id','candidate_label_value','candidate_voxels','candidate_volume_mL', ...
                'SUVmax_QC','SUVmean_QC','SUVmedian_QC','SUVstd_QC', ...
                'centroid_x_voxel','centroid_y_voxel','centroid_z_voxel', ...
                'centroid_x_norm','centroid_y_norm','centroid_z_norm', ...
                'bbox_x_mm','bbox_y_mm','bbox_z_mm','PET_path', ...
                'candidate_labelmap_path','SUV_threshold','minimum_volume_mL', ...
                'connectivity','voxel_volume_mL'});
        end

        if retainedComponents>0
            caseTable = vertcat(caseRows{:});
        else
            caseTable = empty_candidate_table();
        end

        write_labelmap(candidateLabelMap,outputLabelmap,info);

        if ~isempty(CandidateIndex)
            CandidateIndex(CandidateIndex.case_id==caseID,:) = [];
        end
        CandidateIndex = [CandidateIndex; caseTable];

        elapsedSeconds = seconds(datetime('now')-started);
        statusRow = table("DEEP-PSMA","external","external",caseID,petPath, ...
            string(outputLabelmap),"success","",rawComponents, ...
            excludedBelowMinimum,retainedComponents, ...
            max_or_nan(caseTable.candidate_volume_mL), ...
            median_or_nan(caseTable.candidate_volume_mL), ...
            voxelVolume_mL,elapsedSeconds, ...
            string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'VariableNames',{'dataset_name','cohort','split','case_id','PET_path', ...
            'candidate_labelmap_path','status','error_message','raw_hot_components', ...
            'excluded_below_1mL','retained_candidates','maximum_candidate_volume_mL', ...
            'median_candidate_volume_mL','voxel_volume_mL','elapsed_seconds','completed_at'});

        CaseStatus = replace_status_row(CaseStatus,statusRow);
        writetable(CandidateIndex,candidateCheckpoint);
        writetable(CaseStatus,statusCheckpoint);

        fprintf('  PET source: %s\n',petPath);
        fprintf('  Retained candidates: %d\n',retainedComponents);

    catch ME
        elapsedSeconds = seconds(datetime('now')-started);
        statusRow = table("DEEP-PSMA","external","external",caseID,petPath, ...
            string(outputLabelmap),"failed",string(ME.message), ...
            NaN,NaN,NaN,NaN,NaN,NaN,elapsedSeconds, ...
            string(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
            'VariableNames',{'dataset_name','cohort','split','case_id','PET_path', ...
            'candidate_labelmap_path','status','error_message','raw_hot_components', ...
            'excluded_below_1mL','retained_candidates','maximum_candidate_volume_mL', ...
            'median_candidate_volume_mL','voxel_volume_mL','elapsed_seconds','completed_at'});
        CaseStatus = replace_status_row(CaseStatus,statusRow);
        writetable(CaseStatus,statusCheckpoint);
        fprintf(2,'  FAILED: %s\n',ME.message);
    end
end

candidateCsv = fullfile(outDir,"PSMA_objective9A_PET_only_candidate_index_" + runStamp + ".csv");
statusCsv = fullfile(outDir,"PSMA_objective9A_case_status_" + runStamp + ".csv");
writetable(CandidateIndex,candidateCsv);
writetable(CaseStatus,statusCsv);

if saveExcelOutputs
    writetable(CandidateIndex,replace(candidateCsv,".csv",".xlsx"));
    writetable(CaseStatus,replace(statusCsv,".csv",".xlsx"));
end

success = CaseStatus.status=="success";
SummaryMetric = ["Dataset";"Cases passed";"Cases failed";"Retained PET-derived candidates"; ...
    "SUV threshold";"Minimum retained volume mL";"Maximum-volume exclusion"; ...
    "Connectivity";"Ground-truth masks used";"External image folder"];
SummaryValue = ["DEEP-PSMA";string(sum(success));string(sum(~success)); ...
    string(height(CandidateIndex));string(SUV_THRESHOLD);string(MIN_VOLUME_ML); ...
    "NONE";string(CONNECTIVITY);"NO";"PSMA"];
Summary = table(SummaryMetric,SummaryValue);
summaryCsv = fullfile(outDir,"PSMA_objective9A_summary_" + runStamp + ".csv");
writetable(Summary,summaryCsv);
disp(Summary);

diary off;
end

function T = empty_candidate_table()
T = table(strings(0,1),strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1),strings(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    'VariableNames',{'candidate_global_id','dataset_name','cohort','split','case_id', ...
    'candidate_label_value','candidate_voxels','candidate_volume_mL','SUVmax_QC', ...
    'SUVmean_QC','SUVmedian_QC','SUVstd_QC','centroid_x_voxel','centroid_y_voxel', ...
    'centroid_z_voxel','centroid_x_norm','centroid_y_norm','centroid_z_norm', ...
    'bbox_x_mm','bbox_y_mm','bbox_z_mm','PET_path','candidate_labelmap_path', ...
    'SUV_threshold','minimum_volume_mL','connectivity','voxel_volume_mL'});
end

function T = empty_status_table()
T = table(strings(0,1),strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
    strings(0,1),strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
    zeros(0,1),zeros(0,1),zeros(0,1),zeros(0,1),strings(0,1), ...
    'VariableNames',{'dataset_name','cohort','split','case_id','PET_path', ...
    'candidate_labelmap_path','status','error_message','raw_hot_components', ...
    'excluded_below_1mL','retained_candidates','maximum_candidate_volume_mL', ...
    'median_candidate_volume_mL','voxel_volume_mL','elapsed_seconds','completed_at'});
end

function T = replace_status_row(T,newRow)
if ~isempty(T)
    same = T.cohort==newRow.cohort & T.split==newRow.split & T.case_id==newRow.case_id;
    T(same,:) = [];
end
T = [T; newRow];
end

function write_labelmap(labelMap,outPath,referenceInfo)
outInfo = referenceInfo;
outInfo.Datatype = 'uint16';
outInfo.BitsPerPixel = 16;
outPath = string(outPath);
basePath = extractBefore(outPath,strlength(outPath)-6);
niftiwrite(labelMap,char(basePath),outInfo,'Compressed',true);
end

function value = max_or_nan(x)
if isempty(x), value = NaN; else, value = max(double(x),[],'omitnan'); end
end

function value = median_or_nan(x)
if isempty(x), value = NaN; else, value = median(double(x),'omitnan'); end
end
