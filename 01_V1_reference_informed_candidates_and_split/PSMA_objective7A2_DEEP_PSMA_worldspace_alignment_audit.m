function PSMA_objective7A2_DEEP_PSMA_worldspace_alignment_audit()
% PSMA Objective 7A2 — DEEP-PSMA PET/CT world-space alignment audit
%
% Dataset root:
%   D:\DEEP PSMA\DEEP PSMA Extracted
%
% Purpose:
%   Evaluate PSMA PET/CT alignment using full NIfTI affine transforms rather
%   than simple matrix-size x voxel-spacing comparisons.
%
% Per case, the script:
%   1. Reads PSMA PET, CT, and TTB headers.
%   2. Computes PET and CT world-coordinate bounding boxes from all 8 corners.
%   3. Calculates world-space overlap and PET field-of-view coverage by CT.
%   4. Calculates the proportion of TTB tumour voxels that fall inside the
%      CT world-space grid.
%   5. Flags cases suitable for CT-to-PET resampling.
%
% Important:
%   - No image is modified.
%   - PET and CT are allowed to have different matrix sizes and spacing.
%   - TTB must remain aligned to PET.
%
% Outputs:
%   Excel workbook, text report, and MAT workspace.

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 7A2: PET/CT world-space alignment audit\n');
fprintf('============================================================\n\n');

%% ---------------- SETTINGS ----------------
rootFolder = psma_root("PSMA_EXTERNAL_DATA");

if ~isfolder(rootFolder)
    error('Dataset root folder not found:\n%s',rootFolder);
end

projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));
outRoot = fullfile(projectRoot,"06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7A2_WorldSpace_Alignment_Audit");

if ~isfolder(outRoot), mkdir(outRoot); end

stamp = string(datetime('now','Format','yyyyMMdd_HHmmss'));

xlsxPath = fullfile(outRoot, ...
    "PSMA_objective7A2_worldspace_alignment_audit_" + stamp + ".xlsx");

txtPath = fullfile(outRoot, ...
    "PSMA_objective7A2_worldspace_alignment_audit_" + stamp + ".txt");

matPath = fullfile(outRoot, ...
    "PSMA_objective7A2_worldspace_alignment_audit_" + stamp + ".mat");

% Decision rules.
minPETCoverageByCT = 0.90;      % CT covers at least 90% of PET world volume
minTumourCoverageByCT = 0.999;  % at least 99.9% of tumour voxels lie in CT grid
minAxisOverlap = 0.90;          % at least 90% overlap along each axis

fprintf('Dataset root:\n%s\n\n',rootFolder);
fprintf('Output folder:\n%s\n\n',outRoot);

%% ---------------- CASE FOLDERS ----------------
caseDirs = dir(fullfile(rootFolder,"train_*"));
caseDirs = caseDirs([caseDirs.isdir]);

if isempty(caseDirs)
    error('No train_* case folders found.');
end

caseNames = string({caseDirs.name})';
caseNames = sort(caseNames);
nCases = numel(caseNames);

fprintf('Cases found: %d\n\n',nCases);

%% ---------------- PREALLOCATE ----------------
CaseID = caseNames;

PET_Path = strings(nCases,1);
CT_Path = strings(nCases,1);
TTB_Path = strings(nCases,1);

PET_Size = strings(nCases,1);
CT_Size = strings(nCases,1);
PET_Spacing = strings(nCases,1);
CT_Spacing = strings(nCases,1);

PET_Xmin = nan(nCases,1); PET_Xmax = nan(nCases,1);
PET_Ymin = nan(nCases,1); PET_Ymax = nan(nCases,1);
PET_Zmin = nan(nCases,1); PET_Zmax = nan(nCases,1);

CT_Xmin = nan(nCases,1); CT_Xmax = nan(nCases,1);
CT_Ymin = nan(nCases,1); CT_Ymax = nan(nCases,1);
CT_Zmin = nan(nCases,1); CT_Zmax = nan(nCases,1);

Overlap_X_mm = nan(nCases,1);
Overlap_Y_mm = nan(nCases,1);
Overlap_Z_mm = nan(nCases,1);

PET_X_Coverage = nan(nCases,1);
PET_Y_Coverage = nan(nCases,1);
PET_Z_Coverage = nan(nCases,1);
PET_Volume_Coverage_By_CT = nan(nCases,1);

TumourVoxels = nan(nCases,1);
TumourVoxelsInsideCT = nan(nCases,1);
TumourCoverageByCT = nan(nCases,1);

PET_TTB_TransformMatch = false(nCases,1);
PET_CT_AffineFinite = false(nCases,1);

WorldSpaceStatus = strings(nCases,1);
ResamplingSuitability = strings(nCases,1);
ProblemSummary = strings(nCases,1);

%% ---------------- PROCESS CASES ----------------
fprintf('Auditing world-space geometry...\n');

for c = 1:nCases
    caseFolder = fullfile(rootFolder,CaseID(c),"PSMA");

    PET_Path(c) = fullfile(caseFolder,"PET.nii.gz");
    CT_Path(c)  = fullfile(caseFolder,"CT.nii.gz");
    TTB_Path(c) = fullfile(caseFolder,"TTB.nii.gz");

    problems = strings(0,1);

    if ~isfile(PET_Path(c))
        problems(end+1) = "missing PSMA PET";
    end
    if ~isfile(CT_Path(c))
        problems(end+1) = "missing PSMA CT";
    end
    if ~isfile(TTB_Path(c))
        problems(end+1) = "missing PSMA TTB";
    end

    if isempty(problems)
        try
            petInfo = niftiinfo(char(PET_Path(c)));
            ctInfo  = niftiinfo(char(CT_Path(c)));
            ttbInfo = niftiinfo(char(TTB_Path(c)));

            PET_Size(c) = vector_text(petInfo.ImageSize);
            CT_Size(c) = vector_text(ctInfo.ImageSize);
            PET_Spacing(c) = vector_text(petInfo.PixelDimensions);
            CT_Spacing(c) = vector_text(ctInfo.PixelDimensions);

            A_pet = affine_matrix(petInfo);
            A_ct  = affine_matrix(ctInfo);
            A_ttb = affine_matrix(ttbInfo);

            PET_CT_AffineFinite(c) = all(isfinite(A_pet(:))) && ...
                                     all(isfinite(A_ct(:)));

            PET_TTB_TransformMatch(c) = ...
                isequal(size(A_pet),size(A_ttb)) && ...
                all(abs(A_pet(:)-A_ttb(:))<1e-4);

            if ~PET_TTB_TransformMatch(c)
                problems(end+1) = "PET and TTB affine transforms differ";
            end

            petBox = world_bbox(petInfo,A_pet);
            ctBox  = world_bbox(ctInfo,A_ct);

            PET_Xmin(c)=petBox(1,1); PET_Xmax(c)=petBox(1,2);
            PET_Ymin(c)=petBox(2,1); PET_Ymax(c)=petBox(2,2);
            PET_Zmin(c)=petBox(3,1); PET_Zmax(c)=petBox(3,2);

            CT_Xmin(c)=ctBox(1,1); CT_Xmax(c)=ctBox(1,2);
            CT_Ymin(c)=ctBox(2,1); CT_Ymax(c)=ctBox(2,2);
            CT_Zmin(c)=ctBox(3,1); CT_Zmax(c)=ctBox(3,2);

            [Overlap_X_mm(c),PET_X_Coverage(c)] = ...
                interval_overlap(petBox(1,:),ctBox(1,:));

            [Overlap_Y_mm(c),PET_Y_Coverage(c)] = ...
                interval_overlap(petBox(2,:),ctBox(2,:));

            [Overlap_Z_mm(c),PET_Z_Coverage(c)] = ...
                interval_overlap(petBox(3,:),ctBox(3,:));

            petExtent = petBox(:,2)-petBox(:,1);
            overlapExtent = [Overlap_X_mm(c);Overlap_Y_mm(c);Overlap_Z_mm(c)];

            if all(petExtent>0)
                PET_Volume_Coverage_By_CT(c) = ...
                    prod(max(overlapExtent,0))/prod(petExtent);
            end

            % Tumour-voxel coverage in CT grid.
            ttb = niftiread(char(TTB_Path(c)));
            idx = find(ttb~=0);
            TumourVoxels(c) = numel(idx);

            if isempty(idx)
                TumourVoxelsInsideCT(c) = 0;
                TumourCoverageByCT(c) = 1;
            else
                [r,col,s] = ind2sub(size(ttb),idx);

                % MATLAB voxel indices are 1-based. Convert voxel centres
                % to world coordinates using the NIfTI affine.
                voxPET = [double(col(:)),double(r(:)),double(s(:)), ...
                    ones(numel(idx),1)];

                world = voxPET*A_pet;

                % Map world coordinates into CT voxel coordinates.
                voxCT = world/A_ct;

                xct = voxCT(:,1);
                yct = voxCT(:,2);
                zct = voxCT(:,3);

                ctSize = double(ctInfo.ImageSize);

                inside = xct>=0.5 & xct<=ctSize(2)+0.5 & ...
                         yct>=0.5 & yct<=ctSize(1)+0.5 & ...
                         zct>=0.5 & zct<=ctSize(3)+0.5;

                TumourVoxelsInsideCT(c) = sum(inside);
                TumourCoverageByCT(c) = mean(inside);
            end

            axisPass = PET_X_Coverage(c)>=minAxisOverlap && ...
                       PET_Y_Coverage(c)>=minAxisOverlap && ...
                       PET_Z_Coverage(c)>=minAxisOverlap;

            volumePass = PET_Volume_Coverage_By_CT(c)>=minPETCoverageByCT;
            tumourPass = TumourCoverageByCT(c)>=minTumourCoverageByCT;

            if axisPass && volumePass && tumourPass && ...
                    PET_TTB_TransformMatch(c)
                WorldSpaceStatus(c) = "PASS";
                ResamplingSuitability(c) = "SUITABLE";
            elseif tumourPass && PET_TTB_TransformMatch(c)
                WorldSpaceStatus(c) = "REVIEW";
                ResamplingSuitability(c) = ...
                    "LIKELY_SUITABLE_TUMOUR_FULLY_COVERED";

                if ~axisPass
                    problems(end+1) = ...
                        "PET field-of-view not >=90% covered on every axis";
                end
                if ~volumePass
                    problems(end+1) = ...
                        "PET world-volume coverage by CT below 90%";
                end
            else
                WorldSpaceStatus(c) = "FAIL";
                ResamplingSuitability(c) = "NOT_SUITABLE_WITHOUT_REVIEW";

                if ~tumourPass
                    problems(end+1) = ...
                        "CT does not fully cover PET tumour voxels";
                end
            end

        catch ME
            WorldSpaceStatus(c) = "FAIL";
            ResamplingSuitability(c) = "NOT_ASSESSABLE";
            problems(end+1) = "audit error: " + string(ME.message);
        end
    else
        WorldSpaceStatus(c) = "FAIL";
        ResamplingSuitability(c) = "NOT_ASSESSABLE";
    end

    ProblemSummary(c) = strjoin(problems,"; ");

    if mod(c,10)==0 || c==nCases
        fprintf('  Cases processed: %d/%d\n',c,nCases);
    end
end

%% ---------------- TABLES ----------------
CaseSummary = table( ...
    CaseID,PET_Path,CT_Path,TTB_Path, ...
    PET_Size,CT_Size,PET_Spacing,CT_Spacing, ...
    PET_Xmin,PET_Xmax,PET_Ymin,PET_Ymax,PET_Zmin,PET_Zmax, ...
    CT_Xmin,CT_Xmax,CT_Ymin,CT_Ymax,CT_Zmin,CT_Zmax, ...
    Overlap_X_mm,Overlap_Y_mm,Overlap_Z_mm, ...
    PET_X_Coverage,PET_Y_Coverage,PET_Z_Coverage, ...
    PET_Volume_Coverage_By_CT, ...
    TumourVoxels,TumourVoxelsInsideCT,TumourCoverageByCT, ...
    PET_TTB_TransformMatch,PET_CT_AffineFinite, ...
    WorldSpaceStatus,ResamplingSuitability,ProblemSummary);

Problems = CaseSummary(CaseSummary.WorldSpaceStatus~="PASS",:);

SummaryMetric = [ ...
    "Total cases"; ...
    "World-space PASS"; ...
    "World-space REVIEW"; ...
    "World-space FAIL"; ...
    "PET-TTB affine match"; ...
    "Tumour fully covered by CT (>=99.9%)"; ...
    "PET world volume covered by CT >=90%"; ...
    "Suitable for CT-to-PET resampling"; ...
    "Likely suitable because tumour fully covered"];

SummaryValue = [ ...
    nCases; ...
    sum(WorldSpaceStatus=="PASS"); ...
    sum(WorldSpaceStatus=="REVIEW"); ...
    sum(WorldSpaceStatus=="FAIL"); ...
    sum(PET_TTB_TransformMatch); ...
    sum(TumourCoverageByCT>=minTumourCoverageByCT); ...
    sum(PET_Volume_Coverage_By_CT>=minPETCoverageByCT); ...
    sum(ResamplingSuitability=="SUITABLE"); ...
    sum(ResamplingSuitability=="LIKELY_SUITABLE_TUMOUR_FULLY_COVERED")];

Summary = table(SummaryMetric,SummaryValue);

%% ---------------- SAVE ----------------
fprintf('\nWriting outputs...\n');

writetable(Summary,xlsxPath,'Sheet','Summary');
writetable(CaseSummary,xlsxPath,'Sheet','Case_summary');
writetable(Problems,xlsxPath,'Sheet','Problems');

fid = fopen(txtPath,'w');

if fid>=0
    fprintf(fid,'PSMA Objective 7A2 PET/CT world-space alignment audit\n');
    fprintf(fid,'Generated: %s\n',string(datetime('now')));
    fprintf(fid,'Dataset root: %s\n\n',rootFolder);

    for i = 1:height(Summary)
        fprintf(fid,'%-55s %g\n', ...
            Summary.SummaryMetric(i),Summary.SummaryValue(i));
    end

    fprintf(fid,'\nDecision rules:\n');
    fprintf(fid,'Minimum PET volume coverage by CT: %.1f%%\n', ...
        minPETCoverageByCT*100);
    fprintf(fid,'Minimum tumour voxel coverage by CT: %.2f%%\n', ...
        minTumourCoverageByCT*100);
    fprintf(fid,'Minimum overlap on each axis: %.1f%%\n\n', ...
        minAxisOverlap*100);

    fprintf(fid,'Cases requiring review or failing:\n');

    if isempty(Problems)
        fprintf(fid,'None\n');
    else
        for i = 1:height(Problems)
            fprintf(fid,'%s | %s | %s\n', ...
                Problems.CaseID(i), ...
                Problems.ResamplingSuitability(i), ...
                Problems.ProblemSummary(i));
        end
    end

    fclose(fid);
end

save(matPath,'rootFolder','Summary','CaseSummary','Problems', ...
    'minPETCoverageByCT','minTumourCoverageByCT','minAxisOverlap','-v7.3');

fprintf('\n============================================================\n');
fprintf('WORLD-SPACE ALIGNMENT AUDIT COMPLETED\n');
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
% Return [xmin xmax; ymin ymax; zmin zmax] using voxel-edge corners.

sz = double(info.ImageSize(1:3));

% MATLAB NIfTI affine uses voxel coordinates corresponding to image axes.
% Use voxel-edge coordinates to represent the full physical image extent.
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

function txt = vector_text(v)
v = double(v(:))';
txt = "[" + strjoin(compose("%.9g",v)," ") + "]";
end
