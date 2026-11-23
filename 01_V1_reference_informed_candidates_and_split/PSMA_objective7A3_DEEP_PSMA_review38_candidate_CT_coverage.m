function PSMA_objective7A3_DEEP_PSMA_review38_candidate_CT_coverage()
% PSMA Objective 7A3 — Candidate-level CT coverage audit for review cases
%
% Dataset root:
%   D:\DEEP PSMA\DEEP PSMA Extracted
%
% Purpose:
%   Reassess the DEEP-PSMA cases that had:
%       - PET tumour mask coverage by CT >=99.9%
%       - but total PET world-volume coverage by CT <90%
%
% For those cases, this script applies the locked candidate-generation rules:
%       SUV threshold:             >=3
%       Minimum candidate volume:  >=1 mL
%       Maximum candidate volume:  <=1000 mL
%
% For every retained SUV-hot candidate, it calculates:
%       - candidate volume
%       - overlap with reference TTB
%       - centroid inside CT grid
%       - bounding-box corners inside CT grid
%       - proportion of candidate voxels inside CT grid
%       - volume outside CT coverage
%
% Candidate coverage categories:
%       FULLY_COVERED       >=99.9%
%       NEARLY_COVERED      >=95% and <99.9%
%       PARTLY_UNCOVERED    >0% and <95%
%       NOT_COVERED         0%
%
% Case recommendation:
%       PASS_ALL_CANDIDATES_FULLY_COVERED
%       PASS_ONLY_EDGE_NON_TUMOUR_CANDIDATES_AFFECTED
%       REVIEW_TUMOUR_CANDIDATE_NOT_FULLY_COVERED
%       REVIEW_MULTIPLE_CANDIDATES_NOT_FULLY_COVERED
%
% The script does not modify source files.

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 7A3: Candidate-level CT coverage audit\n');
fprintf('============================================================\n\n');

%% ---------------- SETTINGS ----------------
rootFolder = psma_root("PSMA_EXTERNAL_DATA");

if ~isfolder(rootFolder)
    error('Dataset root folder not found:\n%s',rootFolder);
end

projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));
outRoot = fullfile(projectRoot,"06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7A3_Review38_Candidate_CT_Coverage");

if ~isfolder(outRoot), mkdir(outRoot); end

stamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

xlsxPath = fullfile(outRoot, ...
    "PSMA_objective7A3_review38_candidate_CT_coverage_" + stamp + ".xlsx");

txtPath = fullfile(outRoot, ...
    "PSMA_objective7A3_review38_candidate_CT_coverage_" + stamp + ".txt");

matPath = fullfile(outRoot, ...
    "PSMA_objective7A3_review38_candidate_CT_coverage_" + stamp + ".mat");

% Locked candidate-generation settings.
suvThreshold = 3;
minCandidateVolume_mL = 1;
maxCandidateVolume_mL = 1000;
connectivity = 26;

% World-space review selection rules, matching Objective 7A2.
minPETCoverageByCT = 0.90;
minTumourCoverageByCT = 0.999;

% Candidate coverage categories.
fullCoverageThreshold = 0.999;
nearCoverageThreshold = 0.95;

fprintf('Dataset root:\n%s\n\n',rootFolder);
fprintf('Output folder:\n%s\n\n',outRoot);

%% ---------------- CASE FOLDERS ----------------
caseDirs = dir(fullfile(rootFolder,"train_*"));
caseDirs = caseDirs([caseDirs.isdir]);

if isempty(caseDirs)
    error('No train_* case folders found.');
end

caseNames = sort(string({caseDirs.name})');
nAllCases = numel(caseNames);

fprintf('Total cases found: %d\n',nAllCases);
fprintf('Identifying the Objective 7A2 review cases...\n\n');

%% ---------------- OUTPUT ACCUMULATORS ----------------
CandidateResults = table();
CaseRows = cell(nAllCases,1);
reviewCaseCounter = 0;

%% ---------------- PROCESS ALL CASES ----------------
for c = 1:nAllCases
    caseID = caseNames(c);
    psmaFolder = fullfile(rootFolder,caseID,"PSMA");

    petPath = fullfile(psmaFolder,"PET.nii.gz");
    ctPath  = fullfile(psmaFolder,"CT.nii.gz");
    ttbPath = fullfile(psmaFolder,"TTB.nii.gz");

    if ~isfile(petPath) || ~isfile(ctPath) || ~isfile(ttbPath)
        continue;
    end

    try
        petInfo = niftiinfo(char(petPath));
        ctInfo = niftiinfo(char(ctPath));
        ttbInfo = niftiinfo(char(ttbPath));

        A_pet = affine_matrix(petInfo);
        A_ct = affine_matrix(ctInfo);
        A_ttb = affine_matrix(ttbInfo);

        % PET-TTB must be exactly aligned.
        petTtbOK = isequal(double(petInfo.ImageSize),double(ttbInfo.ImageSize)) && ...
            all(abs(A_pet(:)-A_ttb(:))<1e-4);

        if ~petTtbOK
            continue;
        end

        % Reproduce Objective 7A2 case selection.
        petBox = world_bbox(petInfo,A_pet);
        ctBox = world_bbox(ctInfo,A_ct);

        [ox,~] = interval_overlap(petBox(1,:),ctBox(1,:));
        [oy,~] = interval_overlap(petBox(2,:),ctBox(2,:));
        [oz,~] = interval_overlap(petBox(3,:),ctBox(3,:));

        petExtent = petBox(:,2)-petBox(:,1);
        petVolumeCoverage = prod(max([ox;oy;oz],0))/prod(petExtent);

        ttb = niftiread(char(ttbPath));
        tumourIdx = find(ttb~=0);

        if isempty(tumourIdx)
            tumourCoverage = 1;
        else
            [rr,cc,ss] = ind2sub(size(ttb),tumourIdx);
            petVoxelCoords = [double(cc),double(rr),double(ss), ...
                ones(numel(tumourIdx),1)];
            worldCoords = petVoxelCoords*A_pet;
            ctVoxelCoords = worldCoords/A_ct;

            tumourInside = inside_ct_grid(ctVoxelCoords,ctInfo.ImageSize);
            tumourCoverage = mean(tumourInside);
        end

        isReviewCase = tumourCoverage>=minTumourCoverageByCT && ...
                       petVolumeCoverage<minPETCoverageByCT;

        if ~isReviewCase
            continue;
        end

        reviewCaseCounter = reviewCaseCounter+1;

        fprintf('Review case %d: %s\n',reviewCaseCounter,caseID);

        %% Candidate generation
        pet = double(niftiread(char(petPath)));

        if ~isequal(size(pet),size(ttb))
            error('PET and TTB dimensions differ.');
        end

        hotMask = isfinite(pet) & pet>=suvThreshold;
        CC = bwconncomp(hotMask,connectivity);

        voxelVolume_mL = prod(double(petInfo.PixelDimensions(1:3)))/1000;

        nRawComponents = CC.NumObjects;
        retainedCount = 0;

        nFully = 0;
        nNearly = 0;
        nPartial = 0;
        nNone = 0;
        nTumourCandidates = 0;
        nTumourNotFullyCovered = 0;
        nNonTumourNotFullyCovered = 0;

        totalCandidateVolume = 0;
        totalOutsideVolume = 0;
        maximumOutsideFraction = 0;

        for k = 1:CC.NumObjects
            idx = CC.PixelIdxList{k};
            volume_mL = numel(idx)*voxelVolume_mL;

            if volume_mL<minCandidateVolume_mL || ...
                    volume_mL>maxCandidateVolume_mL
                continue;
            end

            retainedCount = retainedCount+1;

            [r,col,s] = ind2sub(size(pet),idx);

            % Candidate voxel centres in PET voxel coordinates.
            petVox = [double(col),double(r),double(s),ones(numel(idx),1)];
            world = petVox*A_pet;
            ctVox = world/A_ct;

            inside = inside_ct_grid(ctVox,ctInfo.ImageSize);
            coverageFraction = mean(inside);
            outsideFraction = 1-coverageFraction;
            outsideVolume_mL = outsideFraction*volume_mL;

            % Candidate centroid.
            centroidPET = [mean(col),mean(r),mean(s),1];
            centroidWorld = centroidPET*A_pet;
            centroidCT = centroidWorld/A_ct;
            centroidInside = inside_ct_grid(centroidCT,ctInfo.ImageSize);

            % Candidate bounding-box corners.
            x = [min(col)-0.5,max(col)+0.5];
            y = [min(r)-0.5,max(r)+0.5];
            z = [min(s)-0.5,max(s)+0.5];

            bboxPET = zeros(8,4);
            q = 0;
            for ix = 1:2
                for iy = 1:2
                    for iz = 1:2
                        q = q+1;
                        bboxPET(q,:) = [x(ix),y(iy),z(iz),1];
                    end
                end
            end

            bboxWorld = bboxPET*A_pet;
            bboxCT = bboxWorld/A_ct;
            bboxAllInside = all(inside_ct_grid(bboxCT,ctInfo.ImageSize));

            % Ground-truth overlap.
            tumourOverlapVoxels = nnz(ttb(idx)~=0);
            tumourOverlapFractionCandidate = tumourOverlapVoxels/numel(idx);
            isTumourCandidate = tumourOverlapVoxels>0;

            if isTumourCandidate
                nTumourCandidates = nTumourCandidates+1;
            end

            if coverageFraction>=fullCoverageThreshold
                coverageCategory = "FULLY_COVERED";
                nFully = nFully+1;
            elseif coverageFraction>=nearCoverageThreshold
                coverageCategory = "NEARLY_COVERED";
                nNearly = nNearly+1;
            elseif coverageFraction>0
                coverageCategory = "PARTLY_UNCOVERED";
                nPartial = nPartial+1;
            else
                coverageCategory = "NOT_COVERED";
                nNone = nNone+1;
            end

            if coverageFraction<fullCoverageThreshold
                if isTumourCandidate
                    nTumourNotFullyCovered = nTumourNotFullyCovered+1;
                else
                    nNonTumourNotFullyCovered = nNonTumourNotFullyCovered+1;
                end
            end

            totalCandidateVolume = totalCandidateVolume+volume_mL;
            totalOutsideVolume = totalOutsideVolume+outsideVolume_mL;
            maximumOutsideFraction = max(maximumOutsideFraction,outsideFraction);

            candidateRow = table( ...
                caseID,retainedCount,k, ...
                numel(idx),volume_mL, ...
                min(pet(idx)),mean(pet(idx)),max(pet(idx)), ...
                tumourOverlapVoxels,tumourOverlapFractionCandidate, ...
                isTumourCandidate, ...
                coverageFraction,outsideFraction,outsideVolume_mL, ...
                logical(centroidInside),bboxAllInside,coverageCategory, ...
                'VariableNames',{ ...
                'CaseID','CandidateLabel','RawComponentIndex', ...
                'VoxelCount','CandidateVolume_mL', ...
                'SUVmin','SUVmean','SUVmax', ...
                'TumourOverlapVoxels','TumourOverlapFractionOfCandidate', ...
                'IsTumourCandidate', ...
                'CTCoverageFraction','OutsideCTFraction','OutsideCTVolume_mL', ...
                'CentroidInsideCT','BoundingBoxFullyInsideCT','CoverageCategory'});

            CandidateResults = [CandidateResults;candidateRow]; %#ok<AGROW>
        end

        nNotFully = nNearly+nPartial+nNone;

        if retainedCount==0
            recommendation = "REVIEW_NO_RETAINED_CANDIDATES";
        elseif nNotFully==0
            recommendation = "PASS_ALL_CANDIDATES_FULLY_COVERED";
        elseif nTumourNotFullyCovered>0
            recommendation = "REVIEW_TUMOUR_CANDIDATE_NOT_FULLY_COVERED";
        elseif nNotFully<=2
            recommendation = "PASS_ONLY_EDGE_NON_TUMOUR_CANDIDATES_AFFECTED";
        else
            recommendation = "REVIEW_MULTIPLE_NON_TUMOUR_CANDIDATES_NOT_FULLY_COVERED";
        end

        CaseRows{reviewCaseCounter} = table( ...
            caseID,petPath,ctPath,ttbPath, ...
            petVolumeCoverage,tumourCoverage, ...
            nRawComponents,retainedCount,nTumourCandidates, ...
            nFully,nNearly,nPartial,nNone,nNotFully, ...
            nTumourNotFullyCovered,nNonTumourNotFullyCovered, ...
            totalCandidateVolume,totalOutsideVolume,maximumOutsideFraction, ...
            recommendation, ...
            'VariableNames',{ ...
            'CaseID','PET_Path','CT_Path','TTB_Path', ...
            'PETWorldVolumeCoverageByCT','TumourVoxelCoverageByCT', ...
            'RawHotComponents','RetainedCandidates','TumourCandidates', ...
            'FullyCoveredCandidates','NearlyCoveredCandidates', ...
            'PartlyUncoveredCandidates','NotCoveredCandidates', ...
            'NotFullyCoveredCandidates', ...
            'TumourCandidatesNotFullyCovered', ...
            'NonTumourCandidatesNotFullyCovered', ...
            'TotalCandidateVolume_mL','TotalOutsideCTVolume_mL', ...
            'MaximumOutsideFraction','Recommendation'});

    catch ME
        reviewCaseCounter = reviewCaseCounter+1;
        fprintf('  ERROR in %s: %s\n',caseID,ME.message);

        CaseRows{reviewCaseCounter} = table( ...
            caseID,petPath,ctPath,ttbPath, ...
            NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN, ...
            NaN,NaN,NaN,"ERROR", ...
            'VariableNames',{ ...
            'CaseID','PET_Path','CT_Path','TTB_Path', ...
            'PETWorldVolumeCoverageByCT','TumourVoxelCoverageByCT', ...
            'RawHotComponents','RetainedCandidates','TumourCandidates', ...
            'FullyCoveredCandidates','NearlyCoveredCandidates', ...
            'PartlyUncoveredCandidates','NotCoveredCandidates', ...
            'NotFullyCoveredCandidates', ...
            'TumourCandidatesNotFullyCovered', ...
            'NonTumourCandidatesNotFullyCovered', ...
            'TotalCandidateVolume_mL','TotalOutsideCTVolume_mL', ...
            'MaximumOutsideFraction','Recommendation'});
    end
end

%% ---------------- COMBINE CASE TABLE ----------------
CaseRows = CaseRows(1:reviewCaseCounter);
CaseRows = CaseRows(~cellfun(@isempty,CaseRows));

if isempty(CaseRows)
    error('No Objective 7A2 review cases were identified.');
end

CaseSummary = vertcat(CaseRows{:});

%% ---------------- SUMMARY ----------------
SummaryMetric = [ ...
    "Review cases identified"; ...
    "Cases all candidates fully covered"; ...
    "Cases only edge non-tumour candidates affected"; ...
    "Cases with tumour candidate not fully covered"; ...
    "Cases with multiple non-tumour candidates affected"; ...
    "Cases with processing error"; ...
    "Total retained candidates"; ...
    "Fully covered candidates"; ...
    "Nearly covered candidates"; ...
    "Partly uncovered candidates"; ...
    "Not covered candidates"; ...
    "Tumour candidates not fully covered"; ...
    "Non-tumour candidates not fully covered"];

SummaryValue = [ ...
    height(CaseSummary); ...
    sum(CaseSummary.Recommendation=="PASS_ALL_CANDIDATES_FULLY_COVERED"); ...
    sum(CaseSummary.Recommendation=="PASS_ONLY_EDGE_NON_TUMOUR_CANDIDATES_AFFECTED"); ...
    sum(CaseSummary.Recommendation=="REVIEW_TUMOUR_CANDIDATE_NOT_FULLY_COVERED"); ...
    sum(CaseSummary.Recommendation=="REVIEW_MULTIPLE_NON_TUMOUR_CANDIDATES_NOT_FULLY_COVERED"); ...
    sum(CaseSummary.Recommendation=="ERROR"); ...
    sum(CaseSummary.RetainedCandidates,'omitnan'); ...
    sum(CaseSummary.FullyCoveredCandidates,'omitnan'); ...
    sum(CaseSummary.NearlyCoveredCandidates,'omitnan'); ...
    sum(CaseSummary.PartlyUncoveredCandidates,'omitnan'); ...
    sum(CaseSummary.NotCoveredCandidates,'omitnan'); ...
    sum(CaseSummary.TumourCandidatesNotFullyCovered,'omitnan'); ...
    sum(CaseSummary.NonTumourCandidatesNotFullyCovered,'omitnan')];

Summary = table(SummaryMetric,SummaryValue);

CasesNeedingReview = CaseSummary(~startsWith(CaseSummary.Recommendation,"PASS_"),:);
CandidatesNotFullyCovered = CandidateResults( ...
    CandidateResults.CTCoverageFraction<fullCoverageThreshold,:);

%% ---------------- SAVE OUTPUTS ----------------
fprintf('\nWriting outputs...\n');

writetable(Summary,xlsxPath,'Sheet','Summary');
writetable(CaseSummary,xlsxPath,'Sheet','Case_summary');
writetable(CandidateResults,xlsxPath,'Sheet','All_candidates');
writetable(CandidatesNotFullyCovered,xlsxPath,'Sheet','Candidates_not_fully_covered');
writetable(CasesNeedingReview,xlsxPath,'Sheet','Cases_needing_review');

fid = fopen(txtPath,'w');

if fid>=0
    fprintf(fid,'PSMA Objective 7A3 candidate-level CT coverage audit\n');
    fprintf(fid,'Generated: %s\n',string(datetime('now')));
    fprintf(fid,'Dataset root: %s\n\n',rootFolder);

    fprintf(fid,'Locked candidate rules:\n');
    fprintf(fid,'  SUV threshold: >= %.2f\n',suvThreshold);
    fprintf(fid,'  Minimum volume: %.2f mL\n',minCandidateVolume_mL);
    fprintf(fid,'  Maximum volume: %.2f mL\n\n',maxCandidateVolume_mL);

    for i = 1:height(Summary)
        fprintf(fid,'%-55s %g\n', ...
            Summary.SummaryMetric(i),Summary.SummaryValue(i));
    end

    fprintf(fid,'\nCases needing review:\n');

    if isempty(CasesNeedingReview)
        fprintf(fid,'None\n');
    else
        for i = 1:height(CasesNeedingReview)
            fprintf(fid,'%s | %s | not fully covered candidates=%g | tumour candidates affected=%g\n', ...
                CasesNeedingReview.CaseID(i), ...
                CasesNeedingReview.Recommendation(i), ...
                CasesNeedingReview.NotFullyCoveredCandidates(i), ...
                CasesNeedingReview.TumourCandidatesNotFullyCovered(i));
        end
    end

    fclose(fid);
end

save(matPath, ...
    'rootFolder','Summary','CaseSummary','CandidateResults', ...
    'CandidatesNotFullyCovered','CasesNeedingReview', ...
    'suvThreshold','minCandidateVolume_mL','maxCandidateVolume_mL', ...
    'fullCoverageThreshold','nearCoverageThreshold','-v7.3');

fprintf('\n============================================================\n');
fprintf('OBJECTIVE 7A3 COMPLETED\n');
fprintf('============================================================\n');

disp(Summary);

fprintf('\nExcel workbook:\n%s\n',xlsxPath);
fprintf('\nText report:\n%s\n',txtPath);
fprintf('\nMAT workspace:\n%s\n',matPath);
fprintf('\nOutput folder:\n%s\n\n',outRoot);

end

%% ========================================================================
function A = affine_matrix(info)
if ~isfield(info,'Transform') || isempty(info.Transform)
    error('NIfTI header has no usable affine transform.');
end

A = double(info.Transform.T);

if ~isequal(size(A),[4 4])
    error('Unexpected affine matrix size.');
end
end

function box = world_bbox(info,A)
sz = double(info.ImageSize(1:3));

x = [0.5,sz(2)+0.5];
y = [0.5,sz(1)+0.5];
z = [0.5,sz(3)+0.5];

corners = zeros(8,4);
k = 0;

for ix = 1:2
    for iy = 1:2
        for iz = 1:2
            k = k+1;
            corners(k,:) = [x(ix),y(iy),z(iz),1];
        end
    end
end

world = corners*A;

box = [min(world(:,1)),max(world(:,1)); ...
       min(world(:,2)),max(world(:,2)); ...
       min(world(:,3)),max(world(:,3))];
end

function [overlap,coverageOfA] = interval_overlap(a,b)
a = sort(double(a(:)));
b = sort(double(b(:)));

overlap = max(0,min(a(2),b(2))-max(a(1),b(1)));
extentA = a(2)-a(1);

if extentA>0
    coverageOfA = overlap/extentA;
else
    coverageOfA = NaN;
end
end

function inside = inside_ct_grid(ctVoxelCoordinates,ctImageSize)
coords = double(ctVoxelCoordinates);
sz = double(ctImageSize);

if isvector(coords)
    coords = reshape(coords,1,[]);
end

inside = coords(:,1)>=0.5 & coords(:,1)<=sz(2)+0.5 & ...
         coords(:,2)>=0.5 & coords(:,2)<=sz(1)+0.5 & ...
         coords(:,3)>=0.5 & coords(:,3)<=sz(3)+0.5;
end
