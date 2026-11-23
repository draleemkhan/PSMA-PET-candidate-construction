function PSMA_objective7B1_DEEP_PSMA_ground_truth_mask_analysis_CORRECTED()
% PSMA Objective 7B1 — DEEP-PSMA ground-truth mask analysis
%
% Purpose:
%   Characterize all PSMA tumour-burden masks (TTB.nii.gz) in the locked
%   95-case external-validation cohort before model application.
%
% Inputs:
%   1) Locked eligibility CSV:
%      ...\Objective7A4_Locked_Eligibility\
%      PSMA_objective7_external_validation_eligibility_LOCKED_95.csv
%
%   2) DEEP-PSMA dataset root:
%      D:\DEEP PSMA\DEEP PSMA Extracted
%
% Per-case outputs:
%   - binary-mask validity
%   - tumour voxel count
%   - total tumour volume (mL)
%   - connected-component count
%   - whole-body SUVmean, SUVmedian, SUVmax
%   - total lesion PSMA (SUV*mL)
%   - edge-touching component count
%   - lesion burden category
%
% Per-component outputs:
%   - lesion/component volume
%   - lesion SUVmean, SUVmedian, SUVmax
%   - TL-PSMA
%   - centroid
%   - bounding box
%   - edge-touching status
%   - volume bin
%
% Important:
%   - Only PSMA files are used.
%   - Excluded cases are not analyzed.
%   - No model predictions are used.
%   - No source files are modified.

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 7B1: Ground-truth mask analysis\n');
fprintf('============================================================\n\n');

%% ---------------- PATHS ----------------
datasetRoot = psma_root("PSMA_EXTERNAL_DATA");

projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));

eligibilityCsv = fullfile(projectRoot,"06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7A4_Locked_Eligibility", ...
    "PSMA_objective7_external_validation_eligibility_LOCKED_95.csv");

outRoot = fullfile(projectRoot,"06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7B1_GroundTruth_Mask_Analysis");

if ~isfolder(datasetRoot)
    error('Dataset root not found:\n%s',datasetRoot);
end

if ~isfile(eligibilityCsv)
    error('Locked eligibility CSV not found:\n%s',eligibilityCsv);
end

if ~isfolder(outRoot), mkdir(outRoot); end

stamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

xlsxPath = fullfile(outRoot, ...
    "PSMA_objective7B1_ground_truth_mask_analysis_" + stamp + ".xlsx");

csvCasePath = fullfile(outRoot, ...
    "PSMA_objective7B1_case_summary_" + stamp + ".csv");

csvLesionPath = fullfile(outRoot, ...
    "PSMA_objective7B1_lesion_summary_" + stamp + ".csv");

matPath = fullfile(outRoot, ...
    "PSMA_objective7B1_ground_truth_mask_analysis_" + stamp + ".mat");

txtPath = fullfile(outRoot, ...
    "PSMA_objective7B1_ground_truth_mask_analysis_" + stamp + ".txt");

%% ---------------- SETTINGS ----------------
connectivity = 26;

% Lesion-volume bins.
volumeBinEdges = [0 0.1 0.5 1 2 5 10 50 100 inf];
volumeBinLabels = [ ...
    "<0.1 mL", ...
    "0.1-<0.5 mL", ...
    "0.5-<1 mL", ...
    "1-<2 mL", ...
    "2-<5 mL", ...
    "5-<10 mL", ...
    "10-<50 mL", ...
    "50-<100 mL", ...
    ">=100 mL"];

% Descriptive burden categories based on reference total tumour volume.
lowBurdenUpper_mL = 10;
moderateBurdenUpper_mL = 100;

%% ---------------- LOAD LOCKED COHORT ----------------
E = readtable(eligibilityCsv,'VariableNamingRule','preserve');
E.Properties.VariableNames = matlab.lang.makeValidName(E.Properties.VariableNames);

must_have(E,["case_id","include_external_validation"],"Locked eligibility CSV");

E.case_id = string(E.case_id);

if islogical(E.include_external_validation)
    includeMask = E.include_external_validation;
elseif isnumeric(E.include_external_validation)
    includeMask = E.include_external_validation~=0;
else
    includeMask = lower(strtrim(string(E.include_external_validation)))=="true" | ...
                  strtrim(string(E.include_external_validation))=="1";
end

caseList = sort(E.case_id(includeMask));

if numel(caseList)~=95
    error('Safety stop: expected 95 included cases, found %d.',numel(caseList));
end

fprintf('Locked PSMA cases to analyze: %d\n\n',numel(caseList));

%% ---------------- ACCUMULATORS ----------------
nCases = numel(caseList);
CaseRows = cell(nCases,1);
LesionRows = cell(nCases,1);
ProblemRows = cell(nCases,1);

%% ---------------- PROCESS CASES ----------------
for c = 1:nCases
    caseID = caseList(c);
    psmaFolder = fullfile(datasetRoot,caseID,"PSMA");

    petPath = fullfile(psmaFolder,"PET.nii.gz");
    ttbPath = fullfile(psmaFolder,"TTB.nii.gz");

    fprintf('Processing %d/%d: %s\n',c,nCases,caseID);

    try
        if ~isfile(petPath), error('Missing PSMA PET.'); end
        if ~isfile(ttbPath), error('Missing PSMA TTB.'); end

        petInfo = niftiinfo(char(petPath));
        ttbInfo = niftiinfo(char(ttbPath));

        pet = double(niftiread(char(petPath)));
        ttbRaw = niftiread(char(ttbPath));
        ttb = ttbRaw~=0;

        if ~isequal(size(pet),size(ttb))
            error('PET and TTB dimensions differ.');
        end

        if ~same_affine(petInfo,ttbInfo,1e-4)
            error('PET and TTB affine transforms differ.');
        end

        spacing = double(ttbInfo.PixelDimensions(1:3));
        voxelVolume_mL = prod(spacing)/1000;

        uniqueValues = unique(double(ttbRaw(:)));
        isBinary = all(ismember(uniqueValues,[0 1]));

        tumourVoxelCount = nnz(ttb);
        tumourVolume_mL = tumourVoxelCount*voxelVolume_mL;

        if tumourVoxelCount>0
            tumourSUV = pet(ttb & isfinite(pet));
            globalSUVmean = mean(tumourSUV,'omitnan');
            globalSUVmedian = median(tumourSUV,'omitnan');
            globalSUVmax = max(tumourSUV,[],'omitnan');
            totalLesionPSMA_SUVmL = sum(tumourSUV,'omitnan')*voxelVolume_mL;
        else
            globalSUVmean = NaN;
            globalSUVmedian = NaN;
            globalSUVmax = NaN;
            totalLesionPSMA_SUVmL = 0;
        end

        CC = bwconncomp(ttb,connectivity);
        nComponents = CC.NumObjects;

        componentVolumes = zeros(nComponents,1);
        edgeTouchCount = 0;
        sub1mLCount = 0;
        ge1mLCount = 0;

        lesionTableCase = table();

        for k = 1:nComponents
            idx = CC.PixelIdxList{k};
            [r,col,s] = ind2sub(size(ttb),idx);

            lesionVolume_mL = numel(idx)*voxelVolume_mL;
            lesionSUV = pet(idx);
            lesionSUV = lesionSUV(isfinite(lesionSUV));

            if isempty(lesionSUV)
                lesionSUVmean = NaN;
                lesionSUVmedian = NaN;
                lesionSUVmax = NaN;
                lesionTLPSMA = NaN;
            else
                lesionSUVmean = mean(lesionSUV,'omitnan');
                lesionSUVmedian = median(lesionSUV,'omitnan');
                lesionSUVmax = max(lesionSUV,[],'omitnan');
                lesionTLPSMA = sum(lesionSUV,'omitnan')*voxelVolume_mL;
            end

            touchesEdge = any(r==1 | r==size(ttb,1) | ...
                              col==1 | col==size(ttb,2) | ...
                              s==1 | s==size(ttb,3));

            if touchesEdge, edgeTouchCount = edgeTouchCount+1; end
            if lesionVolume_mL<1
                sub1mLCount = sub1mLCount+1;
            else
                ge1mLCount = ge1mLCount+1;
            end

            componentVolumes(k) = lesionVolume_mL;

            volumeBin = assign_volume_bin( ...
                lesionVolume_mL,volumeBinEdges,volumeBinLabels);

            centroid_row = mean(r);
            centroid_col = mean(col);
            centroid_slice = mean(s);

            bbox_r_min = min(r); bbox_r_max = max(r);
            bbox_c_min = min(col); bbox_c_max = max(col);
            bbox_s_min = min(s); bbox_s_max = max(s);

            lesionRow = table( ...
                caseID,k,numel(idx),lesionVolume_mL, ...
                lesionSUVmean,lesionSUVmedian,lesionSUVmax,lesionTLPSMA, ...
                centroid_row,centroid_col,centroid_slice, ...
                bbox_r_min,bbox_r_max,bbox_c_min,bbox_c_max, ...
                bbox_s_min,bbox_s_max,touchesEdge,volumeBin, ...
                'VariableNames',{ ...
                'case_id','component_id','voxel_count','volume_mL', ...
                'SUVmean','SUVmedian','SUVmax','TL_PSMA_SUV_mL', ...
                'centroid_row','centroid_col','centroid_slice', ...
                'bbox_row_min','bbox_row_max','bbox_col_min','bbox_col_max', ...
                'bbox_slice_min','bbox_slice_max','touches_image_edge', ...
                'volume_bin'});

            lesionTableCase = [lesionTableCase;lesionRow]; %#ok<AGROW>
        end

        if nComponents>0
            smallestLesion_mL = min(componentVolumes);
            largestLesion_mL = max(componentVolumes);
            medianLesion_mL = median(componentVolumes);
            meanLesion_mL = mean(componentVolumes);
        else
            smallestLesion_mL = NaN;
            largestLesion_mL = NaN;
            medianLesion_mL = NaN;
            meanLesion_mL = NaN;
        end

        if tumourVolume_mL<lowBurdenUpper_mL
            burdenCategory = "LOW_<10mL";
        elseif tumourVolume_mL<moderateBurdenUpper_mL
            burdenCategory = "MODERATE_10-<100mL";
        else
            burdenCategory = "HIGH_>=100mL";
        end

        maskStatus = "PASS";
        problemText = "";

        if ~isBinary
            maskStatus = "REVIEW";
            problemText = "TTB contains values other than 0 and 1.";
        elseif tumourVoxelCount==0
            maskStatus = "REVIEW";
            problemText = "Empty tumour mask.";
        end

        CaseRows{c} = table( ...
            caseID,petPath,ttbPath, ...
            isBinary,string(uniqueValues_text(uniqueValues)), ...
            tumourVoxelCount,voxelVolume_mL,tumourVolume_mL,nComponents, ...
            sub1mLCount,ge1mLCount,edgeTouchCount, ...
            smallestLesion_mL,medianLesion_mL,meanLesion_mL,largestLesion_mL, ...
            globalSUVmean,globalSUVmedian,globalSUVmax,totalLesionPSMA_SUVmL, ...
            burdenCategory,maskStatus,problemText, ...
            'VariableNames',{ ...
            'case_id','PET_path','TTB_path', ...
            'is_binary','unique_mask_values', ...
            'tumour_voxel_count','voxel_volume_mL','total_tumour_volume_mL', ...
            'connected_component_count','components_below_1mL', ...
            'components_1mL_or_greater','edge_touching_components', ...
            'smallest_component_mL','median_component_mL', ...
            'mean_component_mL','largest_component_mL', ...
            'global_SUVmean','global_SUVmedian','global_SUVmax', ...
            'total_lesion_PSMA_SUV_mL','burden_category', ...
            'mask_status','problem_summary'});

        LesionRows{c} = lesionTableCase;

        if maskStatus~="PASS"
            ProblemRows{c} = CaseRows{c};
        end

    catch ME
        CaseRows{c} = table( ...
            caseID,petPath,ttbPath, ...
            false,"",NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            "NOT_ASSESSABLE","FAIL",string(ME.message), ...
            'VariableNames',{ ...
            'case_id','PET_path','TTB_path', ...
            'is_binary','unique_mask_values', ...
            'tumour_voxel_count','voxel_volume_mL','total_tumour_volume_mL', ...
            'connected_component_count','components_below_1mL', ...
            'components_1mL_or_greater','edge_touching_components', ...
            'smallest_component_mL','median_component_mL', ...
            'mean_component_mL','largest_component_mL', ...
            'global_SUVmean','global_SUVmedian','global_SUVmax', ...
            'total_lesion_PSMA_SUV_mL','burden_category', ...
            'mask_status','problem_summary'});

        ProblemRows{c} = CaseRows{c};
        LesionRows{c} = table();
    end
end

%% ---------------- COMBINE TABLES ----------------
CaseSummary = vertcat(CaseRows{:});

nonEmptyLesions = LesionRows(~cellfun(@isempty,LesionRows));
if isempty(nonEmptyLesions)
    LesionSummary = table();
else
    LesionSummary = vertcat(nonEmptyLesions{:});
end

nonEmptyProblems = ProblemRows(~cellfun(@isempty,ProblemRows));
if isempty(nonEmptyProblems)
    Problems = CaseSummary([],:);
else
    Problems = vertcat(nonEmptyProblems{:});
end

%% ---------------- VOLUME-BIN SUMMARY ----------------
if isempty(LesionSummary)
    VolumeBinSummary = table();
else
    [binGroups,~,binIdx] = unique(LesionSummary.volume_bin,'stable');
    lesionCount = accumarray(binIdx,1,[],@sum);
    totalVolume = accumarray(binIdx,LesionSummary.volume_mL,[],@sum);
    medianVolume = accumarray(binIdx,LesionSummary.volume_mL,[],@median);
    meanSUVmax = accumarray(binIdx,LesionSummary.SUVmax,[],@(x) mean(x,'omitnan'));
    meanSUVmean = accumarray(binIdx,LesionSummary.SUVmean,[],@(x) mean(x,'omitnan'));

    VolumeBinSummary = table(binGroups,lesionCount,totalVolume, ...
        medianVolume,meanSUVmean,meanSUVmax, ...
        'VariableNames',{'volume_bin','lesion_count','total_volume_mL', ...
        'median_lesion_volume_mL','mean_lesion_SUVmean','mean_lesion_SUVmax'});
end

%% ---------------- OVERALL SUMMARY ----------------
SummaryMetric = [ ...
    "Locked cases analyzed"; ...
    "Cases passing mask QC"; ...
    "Cases requiring review or failing"; ...
    "Binary masks"; ...
    "Empty masks"; ...
    "Total connected components"; ...
    "Components below 1 mL"; ...
    "Components 1 mL or greater"; ...
    "Edge-touching components"; ...
    "Median case tumour volume mL"; ...
    "Mean case tumour volume mL"; ...
    "Maximum case tumour volume mL"; ...
    "Median components per case"; ...
    "Mean components per case"; ...
    "Median global SUVmean"; ...
    "Median global SUVmax"; ...
    "Median TL-PSMA SUV*mL"];

SummaryValue = [ ...
    height(CaseSummary); ...
    sum(CaseSummary.mask_status=="PASS"); ...
    sum(CaseSummary.mask_status~="PASS"); ...
    sum(CaseSummary.is_binary); ...
    sum(CaseSummary.tumour_voxel_count==0); ...
    sum(CaseSummary.connected_component_count,'omitnan'); ...
    sum(CaseSummary.components_below_1mL,'omitnan'); ...
    sum(CaseSummary.components_1mL_or_greater,'omitnan'); ...
    sum(CaseSummary.edge_touching_components,'omitnan'); ...
    median(CaseSummary.total_tumour_volume_mL,'omitnan'); ...
    mean(CaseSummary.total_tumour_volume_mL,'omitnan'); ...
    max(CaseSummary.total_tumour_volume_mL,[],'omitnan'); ...
    median(CaseSummary.connected_component_count,'omitnan'); ...
    mean(CaseSummary.connected_component_count,'omitnan'); ...
    median(CaseSummary.global_SUVmean,'omitnan'); ...
    median(CaseSummary.global_SUVmax,'omitnan'); ...
    median(CaseSummary.total_lesion_PSMA_SUV_mL,'omitnan')];

Summary = table(SummaryMetric,SummaryValue);

%% ---------------- SAVE OUTPUTS ----------------
fprintf('\nWriting outputs...\n');

writetable(Summary,xlsxPath,'Sheet','Summary');
writetable(CaseSummary,xlsxPath,'Sheet','Case_summary');

if ~isempty(LesionSummary)
    writetable(LesionSummary,xlsxPath,'Sheet','Lesion_summary');
end

if ~isempty(VolumeBinSummary)
    writetable(VolumeBinSummary,xlsxPath,'Sheet','Volume_bins');
end

writetable(Problems,xlsxPath,'Sheet','Problems');

writetable(CaseSummary,csvCasePath);

if ~isempty(LesionSummary)
    writetable(LesionSummary,csvLesionPath);
end

save(matPath, ...
    'CaseSummary','LesionSummary','VolumeBinSummary','Problems','Summary', ...
    'volumeBinEdges','volumeBinLabels','connectivity','-v7.3');

fid = fopen(txtPath,'w');

if fid>=0
    fprintf(fid,'PSMA Objective 7B1 ground-truth mask analysis\n');
    fprintf(fid,'Generated: %s\n',string(datetime('now')));
    fprintf(fid,'Dataset root: %s\n',datasetRoot);
    fprintf(fid,'Locked eligibility CSV: %s\n\n',eligibilityCsv);

    for i = 1:height(Summary)
        fprintf(fid,'%-45s %g\n', ...
            Summary.SummaryMetric(i),Summary.SummaryValue(i));
    end

    fprintf(fid,'\nCases requiring review or failing:\n');

    if isempty(Problems)
        fprintf(fid,'None\n');
    else
        for i = 1:height(Problems)
            fprintf(fid,'%s | %s | %s\n', ...
                Problems.case_id(i), ...
                Problems.mask_status(i), ...
                Problems.problem_summary(i));
        end
    end

    fclose(fid);
end

fprintf('\n============================================================\n');
fprintf('OBJECTIVE 7B1 COMPLETED\n');
fprintf('============================================================\n');

disp(Summary);

fprintf('\nExcel workbook:\n%s\n',xlsxPath);
fprintf('\nCase CSV:\n%s\n',csvCasePath);
fprintf('\nLesion CSV:\n%s\n',csvLesionPath);
fprintf('\nMAT workspace:\n%s\n',matPath);
fprintf('\nText report:\n%s\n',txtPath);
fprintf('\nOutput folder:\n%s\n\n',outRoot);

end

%% ========================================================================
function tf = same_affine(infoA,infoB,tol)
tf = false;

try
    A = double(infoA.Transform.T);
    B = double(infoB.Transform.T);

    tf = isequal(size(A),size(B)) && all(abs(A(:)-B(:))<=tol);
catch
    tf = false;
end
end

function label = assign_volume_bin(v,edges,labels)
idx = find(v>=edges(1:end-1) & v<edges(2:end),1,'first');

if isempty(idx)
    label = "UNCLASSIFIED";
else
    label = labels(idx);
end
end

function txt = uniqueValues_text(v)
v = double(v(:))';
txt = "[" + strjoin(compose("%.9g",v)," ") + "]";
end

function must_have(T,names,source)
available = string(T.Properties.VariableNames);
missing = names(~ismember(names,available));

if ~isempty(missing)
    fprintf('\nColumns detected in %s:\n',source);
    disp(available');

    error('Missing required columns in %s: %s', ...
        source,strjoin(missing,', '));
end
end
