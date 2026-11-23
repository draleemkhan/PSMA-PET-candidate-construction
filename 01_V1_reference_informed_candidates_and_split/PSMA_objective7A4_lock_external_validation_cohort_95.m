function PSMA_objective7A4_lock_external_validation_cohort_95()
% PSMA Objective 7A4 — Lock final 95-case DEEP-PSMA external cohort
%
% Purpose:
%   Create the definitive case-level eligibility file for external
%   validation using only pre-model technical audit findings.
%
% Inclusion rule:
%   Include cases that are:
%     1) World-space PASS in Objective 7A2; OR
%     2) World-space REVIEW in Objective 7A2 and, in Objective 7A3:
%          - all candidates fully covered; OR
%          - only edge non-tumour candidates affected.
%
% Exclusion rule:
%   Exclude cases that are:
%     1) World-space FAIL in Objective 7A2; OR
%     2) Have one or more tumour-overlapping candidates not fully covered
%        by CT in Objective 7A3.
%
% Expected locked result:
%   95 included cases
%   5 excluded cases
%
% Important:
%   - Eligibility is based only on technical image-coverage criteria.
%   - Model predictions and performance are not read or used.
%   - Once generated and verified, this CSV should not be modified after
%     external model evaluation begins.
%
% Outputs:
%   PSMA_objective7_external_validation_eligibility_LOCKED_95.csv
%   Matching XLSX, MAT, TXT and SHA-256 checksum files.

clc;
close all;

fprintf('\n============================================================\n');
fprintf('PSMA Objective 7A4: Lock external validation cohort\n');
fprintf('============================================================\n\n');

%% ---------------- PATHS ----------------
projectRoot = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));

audit7A2Folder = fullfile(projectRoot,"06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7A2_WorldSpace_Alignment_Audit");

audit7A3Folder = fullfile(projectRoot,"06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7A3_Review38_Candidate_CT_Coverage");

outFolder = fullfile(projectRoot,"06_results", ...
    "Objective7_DEEP_PSMA_External_Validation", ...
    "Objective7A4_Locked_Eligibility");

if ~isfolder(outFolder), mkdir(outFolder); end

%% ---------------- FIND LATEST AUDIT FILES ----------------
file7A2 = latest_file(audit7A2Folder, ...
    "PSMA_objective7A2_worldspace_alignment_audit_*.xlsx");

file7A3 = latest_file(audit7A3Folder, ...
    "PSMA_objective7A3_review38_candidate_CT_coverage_*.xlsx");

fprintf('Using Objective 7A2 workbook:\n%s\n\n',file7A2);
fprintf('Using Objective 7A3 workbook:\n%s\n\n',file7A3);

%% ---------------- READ AUDITS ----------------
A2 = readtable(file7A2,'Sheet','Case_summary', ...
    'VariableNamingRule','preserve');

A3 = readtable(file7A3,'Sheet','Case_summary', ...
    'VariableNamingRule','preserve');

A2.Properties.VariableNames = matlab.lang.makeValidName( ...
    A2.Properties.VariableNames);

A3.Properties.VariableNames = matlab.lang.makeValidName( ...
    A3.Properties.VariableNames);

must_have(A2,["CaseID","WorldSpaceStatus", ...
    "ResamplingSuitability","TumourCoverageByCT"],"Objective 7A2");

must_have(A3,["CaseID","Recommendation", ...
    "TumourCandidatesNotFullyCovered", ...
    "NonTumourCandidatesNotFullyCovered"],"Objective 7A3");

A2.CaseID = string(A2.CaseID);
A2.WorldSpaceStatus = string(A2.WorldSpaceStatus);
A2.ResamplingSuitability = string(A2.ResamplingSuitability);

A3.CaseID = string(A3.CaseID);
A3.Recommendation = string(A3.Recommendation);

%% ---------------- DERIVE ELIGIBILITY ----------------
nCases = height(A2);

case_id = A2.CaseID;
include_external_validation = false(nCases,1);
eligibility_status = strings(nCases,1);
exclusion_reason = strings(nCases,1);
technical_basis = strings(nCases,1);

world_space_status = A2.WorldSpaceStatus;
resampling_suitability = A2.ResamplingSuitability;
tumour_CT_coverage_fraction = double(A2.TumourCoverageByCT);

candidate_audit_recommendation = repmat("NOT_APPLICABLE",nCases,1);
tumour_candidates_not_fully_covered = zeros(nCases,1);
non_tumour_candidates_not_fully_covered = zeros(nCases,1);

for i = 1:nCases
    cid = case_id(i);
    idx3 = find(A3.CaseID==cid,1);

    if ~isempty(idx3)
        candidate_audit_recommendation(i) = A3.Recommendation(idx3);
        tumour_candidates_not_fully_covered(i) = ...
            double(A3.TumourCandidatesNotFullyCovered(idx3));
        non_tumour_candidates_not_fully_covered(i) = ...
            double(A3.NonTumourCandidatesNotFullyCovered(idx3));
    end

    if world_space_status(i)=="PASS"
        include_external_validation(i) = true;
        eligibility_status(i) = "INCLUDED";
        technical_basis(i) = ...
            "World-space PASS; PET/CT coverage suitable for locked PET+CT radiomics.";

    elseif world_space_status(i)=="REVIEW"
        if isempty(idx3)
            error('Review case %s is missing from Objective 7A3.',cid);
        end

        rec = candidate_audit_recommendation(i);

        if rec=="PASS_ALL_CANDIDATES_FULLY_COVERED"
            include_external_validation(i) = true;
            eligibility_status(i) = "INCLUDED";
            technical_basis(i) = ...
                "World-space REVIEW resolved: all retained candidates fully covered by CT.";

        elseif rec=="PASS_ONLY_EDGE_NON_TUMOUR_CANDIDATES_AFFECTED"
            include_external_validation(i) = true;
            eligibility_status(i) = "INCLUDED";
            technical_basis(i) = ...
                "World-space REVIEW resolved: only peripheral non-tumour candidate coverage affected.";

        elseif tumour_candidates_not_fully_covered(i)>0 || ...
                rec=="REVIEW_TUMOUR_CANDIDATE_NOT_FULLY_COVERED"
            include_external_validation(i) = false;
            eligibility_status(i) = "EXCLUDED";
            exclusion_reason(i) = ...
                "Incomplete CT coverage for one or more tumour-overlapping SUV-hot candidates.";
            technical_basis(i) = ...
                "Objective 7A3 candidate-level coverage failure.";

        else
            include_external_validation(i) = false;
            eligibility_status(i) = "EXCLUDED";
            exclusion_reason(i) = ...
                "Candidate-level CT coverage did not meet locked technical eligibility criteria.";
            technical_basis(i) = ...
                "Objective 7A3 unresolved technical review.";

        end

    elseif world_space_status(i)=="FAIL"
        include_external_validation(i) = false;
        eligibility_status(i) = "EXCLUDED";
        exclusion_reason(i) = ...
            "Insufficient CT coverage of PET tumour region for valid locked PET+CT radiomics.";
        technical_basis(i) = ...
            "Objective 7A2 world-space FAIL.";

    else
        error('Unexpected WorldSpaceStatus for %s: %s', ...
            cid,world_space_status(i));
    end
end

%% ---------------- BUILD LOCKED TABLE ----------------
lock_date_utc = repmat(string(datetime('now','TimeZone','UTC', ...
    'Format','yyyy-MM-dd''T''HH:mm:ssXXX')),nCases,1);

eligibility_rule_version = repmat( ...
    "DEEP_PSMA_EXTERNAL_ELIGIBILITY_V1_2026-07",nCases,1);

model_status = repmat( ...
    "Eligibility locked before external model application",nCases,1);

LockedEligibility = table( ...
    case_id,include_external_validation,eligibility_status, ...
    exclusion_reason,technical_basis, ...
    world_space_status,resampling_suitability, ...
    tumour_CT_coverage_fraction, ...
    candidate_audit_recommendation, ...
    tumour_candidates_not_fully_covered, ...
    non_tumour_candidates_not_fully_covered, ...
    eligibility_rule_version,lock_date_utc,model_status);

LockedEligibility = sortrows(LockedEligibility,'case_id');

%% ---------------- HARD SAFETY CHECKS ----------------
nIncluded = sum(LockedEligibility.include_external_validation);
nExcluded = sum(~LockedEligibility.include_external_validation);

fprintf('Derived included cases: %d\n',nIncluded);
fprintf('Derived excluded cases: %d\n\n',nExcluded);

if height(LockedEligibility)~=100
    error('Safety stop: expected 100 total cases, found %d.', ...
        height(LockedEligibility));
end

if nIncluded~=95 || nExcluded~=5
    fprintf('Excluded cases detected:\n');
    disp(LockedEligibility(~LockedEligibility.include_external_validation, ...
        {'case_id','exclusion_reason','world_space_status', ...
         'candidate_audit_recommendation'}));

    error(['Safety stop: expected exactly 95 included and 5 excluded, ' ...
        'but derived %d included and %d excluded.'],nIncluded,nExcluded);
end

if any(strlength(LockedEligibility.exclusion_reason( ...
        ~LockedEligibility.include_external_validation))==0)
    error('Safety stop: at least one excluded case lacks an exclusion reason.');
end

if any(strlength(LockedEligibility.exclusion_reason( ...
        LockedEligibility.include_external_validation))>0)
    error('Safety stop: an included case has an exclusion reason.');
end

%% ---------------- WRITE LOCKED FILES ----------------
csvPath = fullfile(outFolder, ...
    "PSMA_objective7_external_validation_eligibility_LOCKED_95.csv");

xlsxPath = fullfile(outFolder, ...
    "PSMA_objective7_external_validation_eligibility_LOCKED_95.xlsx");

matPath = fullfile(outFolder, ...
    "PSMA_objective7_external_validation_eligibility_LOCKED_95.mat");

txtPath = fullfile(outFolder, ...
    "PSMA_objective7_external_validation_eligibility_LOCKED_95_README.txt");

checksumPath = fullfile(outFolder, ...
    "PSMA_objective7_external_validation_eligibility_LOCKED_95_SHA256.txt");

% Refuse silent overwrite of an existing locked CSV.
if isfile(csvPath)
    existing = readtable(csvPath,'VariableNamingRule','preserve');

    if ~isequaln(existing,LockedEligibility)
        backupName = fullfile(outFolder, ...
            "PREVIOUS_LOCKED_FILE_DO_NOT_USE_" + ...
            string(datetime('now','Format','yyyyMMdd_HHmmss')) + ".csv");
        copyfile(csvPath,backupName);

        error(['An existing locked CSV differs from the newly derived ' ...
            'eligibility table. It was backed up to:\n%s\n' ...
            'No new locked file was written.'],backupName);
    else
        fprintf('Existing locked CSV is identical; rewriting companion files only.\n');
    end
else
    writetable(LockedEligibility,csvPath);
end

ExcludedCases = LockedEligibility( ...
    ~LockedEligibility.include_external_validation,:);

IncludedCases = LockedEligibility( ...
    LockedEligibility.include_external_validation,:);

writetable(LockedEligibility,xlsxPath,'Sheet','Locked_eligibility');
writetable(IncludedCases,xlsxPath,'Sheet','Included_95');
writetable(ExcludedCases,xlsxPath,'Sheet','Excluded_5');

SourceAudit = table( ...
    ["Objective 7A2 workbook";"Objective 7A3 workbook"], ...
    [string(file7A2);string(file7A3)], ...
    'VariableNames',{'Source','Path'});

writetable(SourceAudit,xlsxPath,'Sheet','Source_audits');

save(matPath,'LockedEligibility','IncludedCases','ExcludedCases', ...
    'file7A2','file7A3','-v7.3');

%% ---------------- README ----------------
fid = fopen(txtPath,'w');

if fid<0
    error('Could not create README file.');
end

fprintf(fid,'DEEP-PSMA EXTERNAL VALIDATION COHORT LOCK\n');
fprintf(fid,'=========================================\n\n');
fprintf(fid,'Total available cases: 100\n');
fprintf(fid,'Included cases: 95\n');
fprintf(fid,'Excluded cases: 5\n\n');
fprintf(fid,'Lock timestamp (UTC): %s\n',lock_date_utc(1));
fprintf(fid,'Eligibility rule version: %s\n\n', ...
    eligibility_rule_version(1));

fprintf(fid,'Inclusion criteria:\n');
fprintf(fid,['1. Objective 7A2 world-space PASS; or\n' ...
    '2. Objective 7A2 REVIEW resolved by Objective 7A3 showing either:\n' ...
    '   a. all retained SUV-hot candidates fully covered by CT; or\n' ...
    '   b. only peripheral non-tumour candidates affected.\n\n']);

fprintf(fid,'Exclusion criteria:\n');
fprintf(fid,['1. Objective 7A2 world-space FAIL; or\n' ...
    '2. one or more tumour-overlapping candidates not fully covered by CT.\n\n']);

fprintf(fid,['Eligibility was determined before external model application and ' ...
    'without using model predictions or performance.\n\n']);

fprintf(fid,'Excluded cases:\n');

for i = 1:height(ExcludedCases)
    fprintf(fid,'%s | %s\n', ...
        ExcludedCases.case_id(i),ExcludedCases.exclusion_reason(i));
end

fprintf(fid,'\nSource Objective 7A2 workbook:\n%s\n',file7A2);
fprintf(fid,'\nSource Objective 7A3 workbook:\n%s\n',file7A3);

fclose(fid);

%% ---------------- SHA-256 CHECKSUM ----------------
sha256 = file_sha256(csvPath);

fid = fopen(checksumPath,'w');

if fid>=0
    fprintf(fid,'%s  %s\n',sha256, ...
        "PSMA_objective7_external_validation_eligibility_LOCKED_95.csv");
    fclose(fid);
end

%% ---------------- FINAL DISPLAY ----------------
fprintf('\n============================================================\n');
fprintf('EXTERNAL VALIDATION COHORT LOCKED SUCCESSFULLY\n');
fprintf('============================================================\n');
fprintf('Included: %d\n',nIncluded);
fprintf('Excluded: %d\n',nExcluded);
fprintf('SHA-256: %s\n',sha256);

fprintf('\nLocked CSV:\n%s\n',csvPath);
fprintf('\nLocked XLSX:\n%s\n',xlsxPath);
fprintf('\nMAT workspace:\n%s\n',matPath);
fprintf('\nREADME:\n%s\n',txtPath);
fprintf('\nChecksum:\n%s\n\n',checksumPath);

fprintf('Excluded cases:\n');
disp(ExcludedCases(:,{'case_id','exclusion_reason'}));

end

%% ========================================================================
function filePath = latest_file(folderPath,pattern)
if ~isfolder(folderPath)
    error('Required audit folder does not exist:\n%s',folderPath);
end

D = dir(fullfile(folderPath,pattern));

if isempty(D)
    error('No files matching %s were found in:\n%s',pattern,folderPath);
end

[~,idx] = max([D.datenum]);
filePath = string(fullfile(D(idx).folder,D(idx).name));
end

function must_have(T,names,source)
available = string(T.Properties.VariableNames);
missing = names(~ismember(names,available));

if ~isempty(missing)
    fprintf('\nColumns found in %s:\n',source);
    disp(available');

    error('Missing required columns in %s: %s', ...
        source,strjoin(missing,', '));
end
end

function hash = file_sha256(filePath)
md = java.security.MessageDigest.getInstance('SHA-256');
fis = java.io.FileInputStream(java.io.File(char(filePath)));
cleanup = onCleanup(@() fis.close());

buffer = zeros(1,8192,'int8');

while true
    n = fis.read(buffer,0,numel(buffer));

    if n==-1
        break;
    end

    md.update(buffer(1:n));
end

digest = typecast(md.digest(),'uint8');
hash = lower(join(compose('%02x',digest),''));
hash = string(hash);
end
