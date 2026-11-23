function PSMA_objective8B_build_exact_internal_style_external_candidate_table()
% PSMA Objective 8B
% Build the corrected external candidate dataset by reproducing the original
% internal Objective 4 candidate-definition logic exactly:
%
% TRUE candidates:
%   connected components of the ground-truth TTB mask
%   26-connectivity
%   volume >= 1 mL
%
% FALSE candidates:
%   PET SUV >= 3 hot regions
%   outside a 3x3x3 dilation of the TTB mask
%   26-connectivity
%   volume >= 1 mL
%
% NO upper-volume exclusion is applied before feature extraction/inference.
%
% IMPORTANT:
% Objective 8A generated the independent SUV-hot component pool and remains
% useful as an audit. However, the original internal Objective 4 training
% table did NOT simply label all SUV-hot components by overlap. It explicitly
% created TRUE candidates from GT-mask components and FALSE candidates from
% hot regions outside a dilated GT mask. This Objective 8B reproduces that
% exact internal candidate-definition strategy.
%
% Outputs:
%   - Exact internal-style external candidate table
%   - Per-case integer candidate label maps
%   - Case-level QC and class-balance summaries
%   - Large-candidate audit
%   - MAT workspace and text report

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 8B: Exact internal-style external candidates\n');
fprintf('TRUE = GT components | FALSE = SUV-hot outside dilated GT\n');
fprintf('SUV >=3 | volume >=1 mL | no upper limit | 26-connectivity\n');
fprintf('============================================================\n\n');

%% ---------------- SETTINGS ----------------
projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));

SUV_THRESHOLD = 3.0;
MIN_VOLUME_ML = 1.0;
CONNECTIVITY = 26;
DILATE_TUMOUR_MASK = true;

lockedMat = fullfile(projectRoot, ...
    "06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7B1_GroundTruth_Mask_Analysis", ...
    "PSMA_objective7B1_ground_truth_mask_analysis_20260713_081112.mat");

if ~isfile(lockedMat)
    [f,p] = uigetfile('*.mat', ...
        'Select Objective 7B1 MAT containing CaseSummary');
    if isequal(f,0)
        error('Locked 95-case MAT was not selected.');
    end
    lockedMat = string(fullfile(p,f));
end

S = load(lockedMat);

if ~isfield(S,'CaseSummary') || ~istable(S.CaseSummary)
    error('Selected MAT does not contain CaseSummary.');
end

Cases = S.CaseSummary;

if height(Cases)~=95
    error('Safety stop: expected 95 locked cases, found %d.',height(Cases));
end

required = ["case_id","PET_path","TTB_path"];
if ~all(ismember(required,string(Cases.Properties.VariableNames)))
    error('CaseSummary must contain case_id, PET_path and TTB_path.');
end

Cases.case_id = string(Cases.case_id);
Cases.PET_path = string(Cases.PET_path);
Cases.TTB_path = string(Cases.TTB_path);

if numel(unique(Cases.case_id))~=95
    error('Safety stop: locked case IDs are not unique.');
end

outRoot = fullfile(projectRoot, ...
    "06_results", ...
    "Objective8_Corrected_External_Validation", ...
    "Objective8B_Exact_InternalStyle_Candidate_Dataset");

labelMapDir = fullfile(outRoot,"candidate_label_maps");
if ~isfolder(outRoot), mkdir(outRoot); end
if ~isfolder(labelMapDir), mkdir(labelMapDir); end

stamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

candidateCsv = fullfile(outRoot, ...
    "PSMA_objective8B_exact_internal_style_candidates_" + stamp + ".csv");
candidateXlsx = fullfile(outRoot, ...
    "PSMA_objective8B_exact_internal_style_candidates_" + stamp + ".xlsx");
caseCsv = fullfile(outRoot, ...
    "PSMA_objective8B_case_QC_" + stamp + ".csv");
largeCsv = fullfile(outRoot, ...
    "PSMA_objective8B_candidates_above_1000mL_" + stamp + ".csv");
matPath = fullfile(outRoot, ...
    "PSMA_objective8B_workspace_" + stamp + ".mat");
txtPath = fullfile(outRoot, ...
    "PSMA_objective8B_report_" + stamp + ".txt");

fprintf('Locked case source:\n%s\n\n',lockedMat);
fprintf('Output root:\n%s\n\n',outRoot);

%% ---------------- STORAGE ----------------
rows = cell(0,38);
caseRows = cell(0,18);
problemRows = cell(0,4);

globalCandidateID = 0;

%% ---------------- MAIN LOOP ----------------
for i = 1:height(Cases)
    caseID = strtrim(Cases.case_id(i));
    petPath = strtrim(Cases.PET_path(i));
    ttbPath = strtrim(Cases.TTB_path(i));

    fprintf('Case %d/%d: %s\n',i,height(Cases),caseID);

    try
        if ~isfile(petPath)
            error('PET file not found: %s',petPath);
        end
        if ~isfile(ttbPath)
            error('TTB file not found: %s',ttbPath);
        end

        PET = double(niftiread(char(petPath)));
        GT = niftiread(char(ttbPath))>0;

        petInfo = niftiinfo(char(petPath));
        gtInfo = niftiinfo(char(ttbPath));

        if ~isequal(size(PET),size(GT))
            error('PET and TTB dimensions differ.');
        end

        voxelSize_mm = double(petInfo.PixelDimensions(1:3));
        gtVoxelSize_mm = double(gtInfo.PixelDimensions(1:3));

        if any(abs(voxelSize_mm-gtVoxelSize_mm)>1e-4)
            error('PET and TTB voxel spacing differ.');
        end

        voxelVolume_mL = prod(voxelSize_mm)/1000;
        imageSize = size(PET);

        labelMap = zeros(imageSize,'uint32');
        localLabel = 0;

        %% TRUE CANDIDATES FROM GT COMPONENTS
        CCtrue = bwconncomp(GT,CONNECTIVITY);
        nTrueAll = CCtrue.NumObjects;
        nTrueKept = 0;

        for c = 1:CCtrue.NumObjects
            idx = CCtrue.PixelIdxList{c};
            volume_mL = numel(idx)*voxelVolume_mL;

            if volume_mL<MIN_VOLUME_ML
                continue;
            end

            localLabel = localLabel+1;
            globalCandidateID = globalCandidateID+1;
            nTrueKept = nTrueKept+1;
            labelMap(idx) = uint32(localLabel);

            row = make_candidate_row( ...
                globalCandidateID,caseID,localLabel,c, ...
                "true_tumor_mask_component",1,"true_tumor", ...
                nTrueKept,idx,PET,imageSize,voxelSize_mm, ...
                voxelVolume_mL,petPath,ttbPath,"");

            rows(end+1,:) = row; %#ok<AGROW>
        end

        %% FALSE CANDIDATES FROM HOT REGIONS OUTSIDE DILATED GT
        if DILATE_TUMOUR_MASK
            exclusionMask = imdilate(GT,true(3,3,3));
        else
            exclusionMask = GT;
        end

        hotMask = isfinite(PET) & PET>=SUV_THRESHOLD;
        falseHotMask = hotMask & ~exclusionMask;

        CCfalse = bwconncomp(falseHotMask,CONNECTIVITY);
        nFalseAll = CCfalse.NumObjects;
        nFalseKept = 0;

        for c = 1:CCfalse.NumObjects
            idx = CCfalse.PixelIdxList{c};
            volume_mL = numel(idx)*voxelVolume_mL;

            if volume_mL<MIN_VOLUME_ML
                continue;
            end

            localLabel = localLabel+1;
            globalCandidateID = globalCandidateID+1;
            nFalseKept = nFalseKept+1;
            labelMap(idx) = uint32(localLabel);

            row = make_candidate_row( ...
                globalCandidateID,caseID,localLabel,c, ...
                "false_non_tumor_hot_region",0,"false_hot_region", ...
                nFalseKept,idx,PET,imageSize,voxelSize_mm, ...
                voxelVolume_mL,petPath,ttbPath,"");

            rows(end+1,:) = row; %#ok<AGROW>
        end

        outBase = fullfile(labelMapDir, ...
            caseID+"_Objective8B_candidate_labelmap");
        write_label_map(labelMap,outBase,petInfo);
        labelMapPath = outBase+".nii.gz";

        % Populate output label-map path for this case.
        nCaseCandidates = nTrueKept+nFalseKept;
        if nCaseCandidates>0
            firstRow = size(rows,1)-nCaseCandidates+1;
            for rr = firstRow:size(rows,1)
                rows{rr,38} = labelMapPath;
            end
        end

        caseRows(end+1,:) = { ...
            caseID,petPath,ttbPath,labelMapPath, ...
            nTrueAll,nTrueKept,nFalseAll,nFalseKept, ...
            nCaseCandidates, ...
            sum(GT(:)),sum(hotMask(:)),sum(falseHotMask(:)), ...
            voxelVolume_mL,SUV_THRESHOLD,MIN_VOLUME_ML, ...
            double(DILATE_TUMOUR_MASK),"PASS",""}; %#ok<AGROW>

    catch ME
        warning('Case failed: %s | %s',caseID,ME.message);

        caseRows(end+1,:) = { ...
            caseID,petPath,ttbPath,"", ...
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            SUV_THRESHOLD,MIN_VOLUME_ML, ...
            double(DILATE_TUMOUR_MASK),"FAIL",string(ME.message)}; %#ok<AGROW>

        problemRows(end+1,:) = { ...
            caseID,fileparts(petPath),string(ME.identifier), ...
            string(ME.message)}; %#ok<AGROW>
    end
end

%% ---------------- TABLES ----------------
CandidateTable = cell2table(rows, ...
    'VariableNames',{ ...
    'candidate_global_id','case_id','candidate_label_value', ...
    'raw_component_id','candidate_source','candidate_label_binary', ...
    'candidate_label_name','candidate_number_within_source', ...
    'candidate_voxels','candidate_volume_mL', ...
    'SUVmax','SUVmean','SUVmedian','SUVstd', ...
    'centroid_x_voxel','centroid_y_voxel','centroid_z_voxel', ...
    'centroid_x_norm','centroid_y_norm','centroid_z_norm', ...
    'centroid_x_mm','centroid_y_mm','centroid_z_mm', ...
    'bbox_x_mm','bbox_y_mm','bbox_z_mm', ...
    'distance_from_center_mm', ...
    'lateral_zone','AP_zone','craniocaudal_zone', ...
    'anatomical_spatial_zone', ...
    'voxel_volume_mL','PET_path','TTB_path', ...
    'SUV_threshold','minimum_volume_mL', ...
    'maximum_volume_exclusion','candidate_label_map_path'});

CaseQC = cell2table(caseRows, ...
    'VariableNames',{ ...
    'case_id','PET_path','TTB_path','candidate_label_map_path', ...
    'n_true_components_all','n_true_candidates_kept', ...
    'n_false_hot_components_all','n_false_candidates_kept', ...
    'n_total_candidates_kept','gt_voxel_count','hot_voxel_count', ...
    'false_hot_voxel_count','voxel_volume_mL','SUV_threshold', ...
    'minimum_volume_mL','tumour_mask_dilation_used', ...
    'status','error_message'});

if isempty(problemRows)
    Problems = table(strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
        'VariableNames',{'case_id','case_folder','error_identifier','error_message'});
else
    Problems = cell2table(problemRows, ...
        'VariableNames',{'case_id','case_folder','error_identifier','error_message'});
end

LargeCandidates = CandidateTable(CandidateTable.candidate_volume_mL>1000,:);

%% ---------------- SAFETY CHECKS ----------------
if any(CandidateTable.candidate_volume_mL<MIN_VOLUME_ML)
    error('Safety stop: retained candidate below 1 mL.');
end

if any(CandidateTable.candidate_label_binary~=0 & ...
       CandidateTable.candidate_label_binary~=1)
    error('Safety stop: invalid binary class labels.');
end

if height(CandidateTable)~=numel(unique(CandidateTable.candidate_global_id))
    error('Safety stop: duplicate global candidate IDs.');
end

if sum(CaseQC.status=="PASS")~=95
    error('Safety stop: not all 95 cases passed.');
end

%% ---------------- SUMMARY ----------------
SummaryMetric = [
    "Locked cases"
    "Cases passed"
    "Cases failed"
    "True GT components before 1 mL filter"
    "True candidates retained"
    "False hot components before 1 mL filter"
    "False candidates retained"
    "Total candidates retained"
    "Candidates above 1000 mL"
    "True candidates above 1000 mL"
    "False candidates above 1000 mL"
    "Maximum retained candidate volume mL"
    "Median retained candidate volume mL"
    "SUV threshold"
    "Minimum retained volume mL"
    "Maximum-volume exclusion"
    "Connectivity"
    "Tumour-mask dilation"
];

SummaryValue = [
    string(height(Cases))
    string(sum(CaseQC.status=="PASS"))
    string(sum(CaseQC.status=="FAIL"))
    string(sum(CaseQC.n_true_components_all,'omitnan'))
    string(sum(CandidateTable.candidate_label_binary==1))
    string(sum(CaseQC.n_false_hot_components_all,'omitnan'))
    string(sum(CandidateTable.candidate_label_binary==0))
    string(height(CandidateTable))
    string(height(LargeCandidates))
    string(sum(LargeCandidates.candidate_label_binary==1))
    string(sum(LargeCandidates.candidate_label_binary==0))
    string(max(CandidateTable.candidate_volume_mL,[],'omitnan'))
    string(median(CandidateTable.candidate_volume_mL,'omitnan'))
    string(SUV_THRESHOLD)
    string(MIN_VOLUME_ML)
    "NONE"
    string(CONNECTIVITY)
    "3x3x3 logical dilation"
];

Summary = table(SummaryMetric,SummaryValue);

%% ---------------- SAVE ----------------
writetable(CandidateTable,candidateCsv);
writetable(CandidateTable,candidateXlsx);
writetable(CaseQC,caseCsv);
writetable(LargeCandidates,largeCsv);

writetable(Summary,candidateXlsx,'Sheet','Summary');
writetable(CaseQC,candidateXlsx,'Sheet','Case_QC');
writetable(LargeCandidates,candidateXlsx,'Sheet','Candidates_gt1000mL');
writetable(Problems,candidateXlsx,'Sheet','Problems');

save(matPath,'CandidateTable','CaseQC','LargeCandidates','Problems', ...
    'Summary','Cases','SUV_THRESHOLD','MIN_VOLUME_ML', ...
    'CONNECTIVITY','DILATE_TUMOUR_MASK','lockedMat','-v7.3');

fid = fopen(txtPath,'w');
if fid>=0
    fprintf(fid,'PSMA Objective 8B exact internal-style candidate dataset\n');
    fprintf(fid,'Generated: %s\n\n',string(datetime('now')));
    for i = 1:height(Summary)
        fprintf(fid,'%-55s %s\n', ...
            Summary.SummaryMetric(i),Summary.SummaryValue(i));
    end
    fclose(fid);
end

fprintf('\n============================================================\n');
fprintf('OBJECTIVE 8B COMPLETED\n');
fprintf('============================================================\n');
disp(Summary);

fprintf('\nCandidate table:\n%s\n',candidateCsv);
fprintf('\nWorkbook:\n%s\n',candidateXlsx);
fprintf('\nCase QC:\n%s\n',caseCsv);
fprintf('\nLarge-candidate audit:\n%s\n',largeCsv);
fprintf('\nMAT workspace:\n%s\n',matPath);
fprintf('\nText report:\n%s\n',txtPath);
fprintf('\nCandidate label maps:\n%s\n\n',labelMapDir);

fprintf('Next: Objective 8C routine and radiomics feature extraction.\n');

end

%% ========================================================================
function row = make_candidate_row( ...
    globalID,caseID,localLabel,rawComponentID, ...
    source,labelBinary,labelName,numberWithinSource, ...
    idx,PET,imageSize,voxelSize,voxelVolume_mL, ...
    petPath,ttbPath,labelMapPath)

vals = PET(idx);
vals = vals(isfinite(vals));

candidate_voxels = numel(idx);
candidate_volume_mL = candidate_voxels*voxelVolume_mL;

SUVmax = max(vals);
SUVmean = mean(vals,'omitnan');
SUVmedian = median(vals,'omitnan');
SUVstd = std(vals,'omitnan');

[x,y,z] = ind2sub(imageSize,idx);

cx = mean(x);
cy = mean(y);
cz = mean(z);

cxn = cx/imageSize(1);
cyn = cy/imageSize(2);
czn = cz/imageSize(3);

cxmm = cx*voxelSize(1);
cymm = cy*voxelSize(2);
czmm = cz*voxelSize(3);

bboxX = (max(x)-min(x)+1)*voxelSize(1);
bboxY = (max(y)-min(y)+1)*voxelSize(2);
bboxZ = (max(z)-min(z)+1)*voxelSize(3);

centerX = imageSize(1)/2;
centerY = imageSize(2)/2;
centerZ = imageSize(3)/2;

dx = (cx-centerX)*voxelSize(1);
dy = (cy-centerY)*voxelSize(2);
dz = (cz-centerZ)*voxelSize(3);
distance = sqrt(dx^2+dy^2+dz^2);

if cxn<0.45
    lateral = "right";
elseif cxn<=0.55
    lateral = "midline";
else
    lateral = "left";
end

if cyn<0.45
    ap = "posterior";
elseif cyn<=0.55
    ap = "central";
else
    ap = "anterior";
end

if czn<1/3
    cranio = "superior";
elseif czn<2/3
    cranio = "middle";
else
    cranio = "inferior";
end

zone = lateral+"_"+ap+"_"+cranio;

row = { ...
    globalID,caseID,localLabel,rawComponentID, ...
    source,labelBinary,labelName,numberWithinSource, ...
    candidate_voxels,candidate_volume_mL, ...
    SUVmax,SUVmean,SUVmedian,SUVstd, ...
    cx,cy,cz,cxn,cyn,czn,cxmm,cymm,czmm, ...
    bboxX,bboxY,bboxZ,distance, ...
    lateral,ap,cranio,zone, ...
    voxelVolume_mL,petPath,ttbPath, ...
    3.0,1.0,"NONE",labelMapPath};
end

function write_label_map(labelMap,outBase,referenceInfo)
infoOut = referenceInfo;
infoOut.Datatype = 'uint32';
infoOut.BitsPerPixel = 32;

if isfield(infoOut,'Filename'), infoOut.Filename=''; end
if isfield(infoOut,'Filemoddate'), infoOut.Filemoddate=''; end
if isfield(infoOut,'Filesize'), infoOut.Filesize=0; end

try
    niftiwrite(labelMap,char(outBase),infoOut,'Compressed',true);
catch
    niftiwrite(labelMap,char(outBase), ...
        'Compressed',true, ...
        'PixelDimensions',referenceInfo.PixelDimensions(1:3));
end
end
