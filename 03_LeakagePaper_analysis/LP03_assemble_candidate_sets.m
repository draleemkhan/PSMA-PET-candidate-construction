function LP03_assemble_candidate_sets()
%LP03_ASSEMBLE_CANDIDATE_SETS  Leakage paper - build the candidate sets.
%
% Locked methods: LeakagePaper_Methods_LOCKED_v1.0_2026-11-09.md
%   (SHA256 d9ca33362c8459f29bd8d870e5df55c1e63f1d8cfcfe88edf427e4b3695c0fe4)
%
% Builds, with one common predictor schema:
%   internal_B      Pipeline B (PET-only) candidates, all splits      (17,404)
%   internal_Ahist  reference-derived positives + V1 carved negatives (22,794)
%   internal_Astd   reference-derived positives + B negatives with no voxel
%                   within 5 mm (physical ellipsoid) of GT, dropped whole
%   external_B      Pipeline B candidates, locked 95 cohort           (4,173)
%   external_Ahist  reference-derived candidates, locked 95 cohort    (4,658)
%
% Also recomputes the A-std 5 mm exclusion flag per candidate (the
% pre-analysis audit stored counts only). Expected: 30 internal, 0 external.
%
% Predictors (Methods 5.1, 5.2, 6): 15 routine numeric + 3 location zones
% + 107 PET_original_* + 93 CT_original_*. Excluded: scanner, radionuclide,
% contrast, voxel_volume_mL, candidate_voxels(*), *_label_value_verified,
% centroid_*_voxel, CT_resampled_to_PET_grid, and the processing duplicates
% in the V2 master table (PET_rad__/CT_rad__ *_QC copies, SUV_threshold,
% minimum_volume_mL, connectivity).
%
% READ-ONLY on all existing files. Writes only to a new folder:
%   06_results\LeakagePaper_ModelData\LP03_<timestamp>\
% No model is trained and no outcome is evaluated here.

clc;
t0 = tic;
V2 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project_V2_PET_Derived"));
V1 = string(fullfile(psma_root("PSMA_PROJECT_ROOT"),"PSMA_AI_Project"));
R9 = fullfile(V2,"06_results","Objective9_PET_Derived_Candidate_Pipeline");
RA = fullfile(V2,"06_results","Objective9_LockedExternal_Addendum_train0096_0100");
FE = fullfile(V2,"06_results","LeakagePaper_FeatureExtraction");
C5 = fullfile(V1,"06_results","Objective5C FINAL-FINAL-FINALL routine plus_Radiomics");

src = struct();
src.master9F   = fullfile(R9,"Objective9F_MASTER_Candidate_Dataset","internal","PSMA_objective9F_internal_MASTER_candidate_dataset.csv");
src.ext9B      = fullfile(R9,"Objective9B_GT_Overlap_and_Label_Assignment","external","external","PSMA_objective9B_candidate_overlap_table_20260716_090139.csv");
src.ext9B_add  = fullfile(RA,"Objective9B_GT_Overlap_and_Label_Assignment","external","external","PSMA_objective9B_candidate_overlap_table_20261109_150612.csv");
src.ext9C      = fullfile(R9,"Objective9C_Routine_Feature_Extraction","external","external","PSMA_objective9C_external_external_routine_features_checkpoint.csv");
src.ext9C_add  = fullfile(RA,"Objective9C_Routine_Feature_Extraction","external","external","PSMA_objective9C_external_external_routine_features_checkpoint.csv");
src.extPET     = fullfile(R9,"Objective9D9E_Direct_PET_CT_Radiomics","external","PET","PSMA_objective9PET_external_PET_radiomics_checkpoint.csv");
src.extPET_add = fullfile(RA,"Objective9D9E_Direct_PET_CT_Radiomics","external","PET","PSMA_objective9PET_external_PET_radiomics_checkpoint.csv");
src.extCT      = fullfile(R9,"Objective9D9E_Direct_PET_CT_Radiomics","external","CT","PSMA_objective9CT_external_CT_radiomics_checkpoint.csv");
src.extCT_add  = fullfile(RA,"Objective9D9E_Direct_PET_CT_Radiomics","external","CT","PSMA_objective9CT_external_CT_radiomics_checkpoint.csv");
src.lp02_int   = fullfile(FE,"routine_internal","LP02_routine_internal_checkpoint.csv");
src.lp02_ext   = fullfile(FE,"routine_external","LP02_routine_external_checkpoint.csv");
src.lp01_int   = fullfile(FE,"internal_ct","LP01_internal_ct_CT_radiomics_checkpoint.csv");
src.lp01_ext   = fullfile(FE,"external_ct","LP01_external_ct_CT_radiomics_checkpoint.csv");
src.v1PET_int  = fullfile(C5,"PSMA_objective5C_PET_radiomics_FULL.csv");
src.v1CT_int   = fullfile(C5,"PSMA_objective5C_CT_radiomics_FULL.csv");
src.v1PET_ext  = fullfile(V1,"06_results","Objective8_Corrected_External_Validation","Objective8D8E_Resumable_PET_CT_Radiomics","PSMA_objective8PET_pet_radiomics_checkpoint.csv");
src.lock95     = fullfile(V1,"06_results","Objective7_DEEP_PSMA_External_Validation","Objective7A4_Locked_Eligibility","PSMA_objective7_external_validation_eligibility_LOCKED_95.csv");

fn = string(fieldnames(src));
for i = 1:numel(fn)
    if ~isfile(src.(fn(i))), error("Missing input (%s):\n%s", fn(i), src.(fn(i))); end
end

outDir = fullfile(V2,"06_results","LeakagePaper_ModelData","LP03_" + string(datetime("now","Format","yyyyMMdd_HHmmss")));
mkdir(outDir);
logFile = fullfile(outDir,"LP03_log.txt");
diary(char(logFile)); cleanupDiary = onCleanup(@() diary("off"));

fprintf("LP03 assembly started %s\nMATLAB %s\nOutput: %s\n\n", string(datetime("now")), version, outDir);

routineNum = ["candidate_volume_mL","SUVmax","SUVmean","SUVmedian","SUVstd", ...
    "CTmean","CTmedian","CTstd","centroid_x_norm","centroid_y_norm","centroid_z_norm", ...
    "distance_from_center_mm","bbox_x_mm","bbox_y_mm","bbox_z_mm"];
zoneVars = ["lateral_zone","AP_zone","craniocaudal_zone"];
marginMm = 5.0;

%% ---------------- Internal B (V2 9F master) ----------------------------
fprintf("Reading V2 9F internal master...\n");
M = readCsvLP(src.master9F);
petNames = extractPrefixed(M, "PET_rad__PET_original_", "PET_rad__");
ctNames  = extractPrefixed(M, "CT_rad__CT_original_",  "CT_rad__");
checkCount("PET radiomics features (9F)", numel(petNames), 107);
checkCount("CT radiomics features (9F)",  numel(ctNames), 93);
predNumeric = [routineNum, petNames, ctNames];

F = table();
for v = routineNum, F.(v) = double(M.(v)); end
for v = zoneVars,   F.(v) = string(M.(v)); end
for v = petNames,   F.(v) = double(M.("PET_rad__" + v)); end
for v = ctNames,    F.(v) = double(M.("CT_rad__" + v)); end

iB = makeMetaB(M, "internal");
iB = [iB, F];
iB.routine_source = repmat("V2_9C", height(iB),1);
iB.PET_source = repmat("V2_9D", height(iB),1);
iB.CT_source = repmat("V2_9E", height(iB),1);
reportBMismatch(M, iB, "internal");

checkCount("internal B rows", height(iB), 17404);
checkCount("internal B positives", sum(iB.y==1), 8200);
checkCount("internal B negatives", sum(iB.y==0), 9204);
checkCount("internal B test rows", sum(iB.split=="test"), 2673);
checkCount("internal B test positives", sum(iB.split=="test" & iB.y==1), 1263);
checkCount("internal test patients", numel(unique(iB.patient_id(iB.split=="test"))), 58);
checkCount("internal test studies", numel(unique(iB.case_id(iB.split=="test"))), 87);
clear M F;

%% ---------------- A-std 5 mm flags (internal B negatives) --------------
fprintf("\nComputing 5 mm exclusion flags for internal B negatives (%d cases)...\n", numel(unique(iB.case_id)));
[iB.Astd_excluded_5mm, touchRawI, caseSumI] = fiveMmFlags(iB, marginMm);
checkCount("internal B negatives touching raw GT (must be 0)", sum(touchRawI), 0);
checkCount("internal B negatives excluded by 5 mm rule", sum(iB.Astd_excluded_5mm), 30);
writetable(caseSumI, fullfile(outDir,"LP03_Astd_5mm_case_summary_internal.csv"));

%% ---------------- Internal reference-derived (A) -----------------------
fprintf("\nReading internal reference-derived features (LP02, V1 5C PET, V1 5C CT, LP01 CT)...\n");
R  = readCsvLP(src.lp02_int);
P1 = readCsvLP(src.v1PET_int);
C1 = readCsvLP(src.v1CT_int);
CL = readCsvLP(src.lp01_int);
checkCount("LP02 internal rows", height(R), 22794);
checkCount("LP01 internal CT rows", height(CL), 1970);
checkCount("LP01 internal CT studies", numel(unique(CL.case_id)), 107);
mustHaveAll(P1, petNames, "V1 5C PET"); mustHaveAll(C1, ctNames, "V1 5C CT"); mustHaveAll(CL, ctNames, "LP01 internal CT");

key = string(R.candidate_global_id);
P1 = alignByKey(P1, key, R.case_id, "V1 5C PET");
C1 = alignByKey(C1, key, R.case_id, "V1 5C CT");
if any(double(P1.candidate_label_binary) ~= double(R.v1_candidate_label_binary))
    error("V1 label mismatch between LP02 and V1 5C PET table.");
end
useLP01 = ismember(key, string(CL.candidate_global_id));
checkCount("internal A rows using LP01 CT", sum(useLP01), 1970);
% every candidate of a re-extracted study must come from LP01
if any(ismember(R.case_id(~useLP01), unique(CL.case_id)))
    error("Some candidates of LP01 studies are missing from LP01 output.");
end
CLa = alignByKey(CL, key(useLP01), R.case_id(useLP01), "LP01 internal CT");

F = table();
for v = routineNum, F.(v) = double(R.(v)); end
for v = zoneVars,   F.(v) = string(R.(v)); end
for v = petNames,   F.(v) = double(P1.(v)); end
for v = ctNames
    x = double(C1.(v));
    x(useLP01) = double(CLa.(v));
    F.(v) = x;
end

iA = table();
iA.row_id = "A_INT|" + key;
iA.source_candidate_id = key;
iA.cohort = repmat("internal", height(R),1);
iA.split = normSplit(R.split);
iA.case_id = string(R.case_id);
iA.patient_id = patientId(iA.case_id, "internal");
iA.candidate_label_value = double(R.candidate_label_value);
iA.candidate_labelmap_path = string(R.candidate_labelmap_path);
iA.GT_path = strings(height(R),1);
iA.y = double(R.v1_candidate_label_binary);
iA.f_overlap = NaN(height(R),1);
iA.overlap_volume_mL = NaN(height(R),1);
iA.label_any = NaN(height(R),1);
iA.label_25pct = NaN(height(R),1);
iA.label_highpurity = NaN(height(R),1);
iA.occupancy_category = strings(height(R),1);
iA.Astd_excluded_5mm = false(height(R),1);
iA = [iA, F];
iA.routine_source = repmat("LP02_V2_9C_on_V1_mask", height(iA),1);
iA.PET_source = repmat("V1_5C_reused", height(iA),1);
iA.CT_source = repmat("V1_5C_reused", height(iA),1);
iA.CT_source(useLP01) = "LP01_V2_rule";
clear R P1 C1 CL CLa F;

checkCount("internal A-hist rows", height(iA), 22794);
checkCount("internal A positives", sum(iA.y==1), 11380);
checkCount("internal A-hist negatives", sum(iA.y==0), 11414);
checkSplitConsistency(iB, iA);

% A-hist = all reference-derived candidates
iAhist = iA;
iAhist.pipeline = repmat("Ahist", height(iAhist),1);

% A-std = reference-derived positives + B negatives retained by 5 mm rule
pos = iA(iA.y==1,:);
negB = iB(iB.y==0 & ~iB.Astd_excluded_5mm,:);
negB = negB(:, pos.Properties.VariableNames);
iAstd = [pos; negB];
iAstd.pipeline = repmat("Astd", height(iAstd),1);
checkCount("internal A-std negatives", sum(iAstd.y==0), 9174);

iB.pipeline = repmat("B", height(iB),1);

%% ---------------- External locked 95 -----------------------------------
L95 = readCsvLP(src.lock95);
lockCases = string(L95.case_id(double(L95.include_external_validation)==1));
checkCount("locked external cases", numel(lockCases), 95);

fprintf("\nReading external Pipeline B tables (original 0001-0095 + addendum 0096-0100)...\n");
B9 = stackTables(readCsvLP(src.ext9B), readCsvLP(src.ext9B_add));
C9 = stackTables(readCsvLP(src.ext9C), readCsvLP(src.ext9C_add));
PE = stackTables(readCsvLP(src.extPET), readCsvLP(src.extPET_add));
CE = stackTables(readCsvLP(src.extCT), readCsvLP(src.extCT_add));
mustHaveAll(PE, petNames, "V2 external PET"); mustHaveAll(CE, ctNames, "V2 external CT");

B9 = B9(ismember(string(B9.case_id), lockCases),:);
ekey = string(B9.candidate_global_id);
C9 = alignByKey(C9, ekey, B9.case_id, "V2 external 9C");
PE = alignByKey(PE, ekey, B9.case_id, "V2 external PET");
CE = alignByKey(CE, ekey, B9.case_id, "V2 external CT");

F = table();
for v = routineNum, F.(v) = double(C9.(v)); end
for v = zoneVars,   F.(v) = string(C9.(v)); end
for v = petNames,   F.(v) = double(PE.(v)); end
for v = ctNames,    F.(v) = double(CE.(v)); end
eB = makeMetaB(B9, "external");
eB = [eB, F];
eB.routine_source = repmat("V2_9C", height(eB),1);
eB.PET_source = repmat("V2_9D", height(eB),1);
eB.CT_source = repmat("V2_9E", height(eB),1);
reportBMismatch(B9, eB, "external");
checkCount("external B rows", height(eB), 4173);
checkCount("external B positives", sum(eB.y==1), 3026);
checkCount("external B studies", numel(unique(eB.case_id)), 95);
clear B9 C9 PE CE F;

fprintf("\nComputing 5 mm exclusion flags for external B negatives (audit only)...\n");
[eB.Astd_excluded_5mm, touchRawE, caseSumE] = fiveMmFlags(eB, marginMm);
checkCount("external B negatives touching raw GT (must be 0)", sum(touchRawE), 0);
checkCount("external B negatives excluded by 5 mm rule", sum(eB.Astd_excluded_5mm), 0);
writetable(caseSumE, fullfile(outDir,"LP03_Astd_5mm_case_summary_external.csv"));
eB.pipeline = repmat("B", height(eB),1);

fprintf("\nReading external reference-derived features (LP02, V1 8 PET, LP01 CT)...\n");
R  = readCsvLP(src.lp02_ext);
P8 = readCsvLP(src.v1PET_ext);
CL = readCsvLP(src.lp01_ext);
P8 = renamePrefix(P8, "original_", "PET_original_");
mustHaveAll(P8, petNames, "V1 8 PET"); mustHaveAll(CL, ctNames, "LP01 external CT");
R = R(ismember(string(R.case_id), lockCases),:);
key = string(R.candidate_global_id);
P8 = alignByKey(P8, key, R.case_id, "V1 8 PET");
CL = alignByKey(CL, key, R.case_id, "LP01 external CT");
if any(double(P8.candidate_label_binary) ~= double(R.v1_candidate_label_binary))
    error("V1 label mismatch between LP02 external and V1 8 PET table.");
end
F = table();
for v = routineNum, F.(v) = double(R.(v)); end
for v = zoneVars,   F.(v) = string(R.(v)); end
for v = petNames,   F.(v) = double(P8.(v)); end
for v = ctNames,    F.(v) = double(CL.(v)); end
eA = table();
eA.row_id = "A_EXT|" + key;
eA.source_candidate_id = key;
eA.cohort = repmat("external", height(R),1);
eA.split = repmat("external", height(R),1);
eA.case_id = string(R.case_id);
eA.patient_id = patientId(eA.case_id, "external");
eA.candidate_label_value = double(R.candidate_label_value);
eA.candidate_labelmap_path = string(R.candidate_labelmap_path);
eA.GT_path = strings(height(R),1);
eA.y = double(R.v1_candidate_label_binary);
eA.f_overlap = NaN(height(R),1);
eA.overlap_volume_mL = NaN(height(R),1);
eA.label_any = NaN(height(R),1);
eA.label_25pct = NaN(height(R),1);
eA.label_highpurity = NaN(height(R),1);
eA.occupancy_category = strings(height(R),1);
eA.Astd_excluded_5mm = false(height(R),1);
eA = [eA, F];
eA.routine_source = repmat("LP02_V2_9C_on_V1_mask", height(eA),1);
eA.PET_source = repmat("V1_8_reused", height(eA),1);
eA.CT_source = repmat("LP01_V2_rule", height(eA),1);
eA.pipeline = repmat("Ahist", height(eA),1);
checkCount("external A-hist rows", height(eA), 4658);
checkCount("external A positives", sum(eA.y==1), 3250);
checkCount("external A-hist studies", numel(unique(eA.case_id)), 95);
clear R P8 CL F;

%% ---------------- Common column order and save ------------------------
metaOrder = ["row_id","source_candidate_id","pipeline","cohort","split","case_id","patient_id", ...
    "candidate_label_value","candidate_labelmap_path","GT_path","y","f_overlap","overlap_volume_mL", ...
    "label_any","label_25pct","label_highpurity","occupancy_category","Astd_excluded_5mm", ...
    "routine_source","PET_source","CT_source"];
colOrder = [metaOrder, routineNum, zoneVars, petNames, ctNames];
data = struct();
data.internal_B     = iB(:, colOrder);
data.internal_Astd  = iAstd(:, colOrder);
data.internal_Ahist = iAhist(:, colOrder);
data.external_B     = eB(:, colOrder);
data.external_Ahist = eA(:, colOrder);

predictors = struct();
predictors.numeric = predNumeric;
predictors.categorical = zoneVars;
predictors.routine_numeric = routineNum;
predictors.PET = petNames;
predictors.CT = ctNames;
predictors.ablation_removed = ["PET_original_firstorder_Minimum","PET_original_firstorder_10Percentile","PET_original_firstorder_Range"];
if ~all(ismember(predictors.ablation_removed, petNames)), error("Ablation predictor names not found."); end

info = struct();
info.created = string(datetime("now"));
info.matlab_version = string(version);
info.sources = src;
info.locked_methods_sha256 = "d9ca33362c8459f29bd8d870e5df55c1e63f1d8cfcfe88edf427e4b3695c0fe4";
info.margin_mm = marginMm;

% counts table
sets = string(fieldnames(data));
cnt = table();
for s = sets'
    T = data.(s);
    for sp = unique(T.split)'
        m = T.split == sp;
        cnt = [cnt; table(s, sp, sum(m), sum(T.y(m)==1), sum(T.y(m)==0), numel(unique(T.case_id(m))), numel(unique(T.patient_id(m))), ...
            'VariableNames',{'set','split','n','n_pos','n_neg','n_studies','n_patients'})]; %#ok<AGROW>
    end
end
disp(cnt);
writetable(cnt, fullfile(outDir,"LP03_counts.csv"));

fprintf("\nWriting CSV copies...\n");
for s = sets'
    writetable(data.(s), fullfile(outDir, "LP03_" + s + ".csv"));
end
flags = iB(iB.y==0, ["row_id","case_id","split","candidate_label_value","Astd_excluded_5mm"]);
writetable(flags, fullfile(outDir,"LP03_Astd_5mm_flags_internal_B_negatives.csv"));
flagsE = eB(eB.y==0, ["row_id","case_id","split","candidate_label_value","Astd_excluded_5mm"]);
writetable(flagsE, fullfile(outDir,"LP03_Astd_5mm_flags_external_B_negatives.csv"));

matFile = fullfile(outDir,"LP03_candidate_sets.mat");
save(matFile, "data","predictors","info","cnt","-v7.3");

fprintf("\n============================================================\n");
fprintf("LP03 COMPLETED (%.1f min). All count checks passed.\n", toc(t0)/60);
fprintf("MAT for LP04: %s\n", matFile);
fprintf("============================================================\n");
end

%% ======================================================================
function T = readCsvLP(f)
    opts = detectImportOptions(f, "FileType","text","Delimiter",",","VariableNamingRule","preserve");
    names = string(opts.VariableNames);
    names = erase(names, char(65279));
    names = erase(names, string(char([239 187 191])));
    names = strtrim(erase(names, '"'));
    opts.VariableNames = cellstr(names);
    strVars = ["candidate_global_id","case_id","split","dataset_name","cohort","candidate_labelmap_path", ...
        "candidate_label_map_path","GT_path","PET_path","CT_path","lateral_zone","AP_zone","craniocaudal_zone", ...
        "overlap_category","source_table","processing_status","error_message","candidate_source","split_set", ...
        "Image","Mask","modality","radiomics_status","radiomics_modality","radiomics_image_path", ...
        "scanner","radionuclide","contrast","eligibility_status","exclusion_reason"];
    sv = intersect(strVars, names, "stable");
    if ~isempty(sv), opts = setvartype(opts, cellstr(sv), "string"); end
    numKnown = ["candidate_volume_mL","SUVmax","SUVmean","SUVmedian","SUVstd","CTmean","CTmedian","CTstd", ...
        "centroid_x_norm","centroid_y_norm","centroid_z_norm","distance_from_center_mm","bbox_x_mm","bbox_y_mm","bbox_z_mm", ...
        "candidate_label_value","candidate_overlap_fraction","overlap_volume_mL","candidate_label_any_overlap", ...
        "candidate_label_overlap_25pct","candidate_label_overlap_50pct","v1_candidate_label_binary", ...
        "candidate_label_binary","include_external_validation"];
    nv = names(startsWith(names, ["PET_original_","CT_original_","original_","PET_rad__PET_original_","CT_rad__CT_original_"]) ...
        | ismember(names, numKnown));
    if ~isempty(nv), opts = setvartype(opts, cellstr(nv), "double"); end
    T = readtable(f, opts);
    fprintf("  read %-60s %6d rows x %3d cols\n", shortName(f), height(T), width(T));
end

function s = shortName(f)
    [~,a,b] = fileparts(f); s = string(a) + string(b);
    if strlength(s) > 60, s = extractBefore(s, 58) + ".."; end
end

function names = extractPrefixed(T, prefix, strip)
    v = string(T.Properties.VariableNames);
    names = erase(v(startsWith(v, prefix)), strip);
end

function mustHaveAll(T, names, label)
    miss = setdiff(names, string(T.Properties.VariableNames));
    if ~isempty(miss)
        error("%s is missing %d feature columns, e.g. %s", label, numel(miss), miss(1));
    end
end

function T = renamePrefix(T, oldP, newP)
    v = string(T.Properties.VariableNames);
    m = startsWith(v, oldP);
    v(m) = newP + extractAfter(v(m), strlength(oldP));
    T.Properties.VariableNames = cellstr(v);
end

function S = stackTables(A, B)
    va = string(A.Properties.VariableNames); vb = string(B.Properties.VariableNames);
    common = intersect(va, vb, "stable");
    S = [A(:, common); B(:, common)];
    if numel(unique(string(S.candidate_global_id))) ~= height(S)
        error("Duplicate candidate IDs after stacking original and addendum tables.");
    end
end

function T = alignByKey(T, key, caseIds, label)
    k = string(T.candidate_global_id);
    if numel(unique(k)) ~= numel(k), error("%s: duplicate candidate IDs.", label); end
    [tf, loc] = ismember(key, k);
    if ~all(tf), error("%s: %d of %d candidates missing.", label, sum(~tf), numel(key)); end
    T = T(loc,:);
    if ismember("case_id", string(T.Properties.VariableNames))
        if any(string(T.case_id) ~= string(caseIds))
            error("%s: case_id mismatch after alignment.", label);
        end
    end
end

function sp = normSplit(x)
    sp = lower(strtrim(string(x)));
    sp(sp=="val") = "validation";
end

function p = patientId(caseIds, cohort)
    caseIds = string(caseIds);
    if cohort == "internal"
        ok = ~cellfun(@isempty, regexp(cellstr(caseIds), '_\d{4}-\d{2}-\d{2}$', 'once'));
        if ~all(ok), error("Unexpected internal case_id format, e.g. %s", caseIds(find(~ok,1))); end
        p = string(regexprep(cellstr(caseIds), '_\d{4}-\d{2}-\d{2}$', ''));
    else
        p = caseIds;   % one study per DEEP-PSMA case
    end
end

function B = makeMetaB(M, cohort)
    n = height(M);
    id = string(M.candidate_global_id);
    B = table();
    B.row_id = "B|" + id;
    B.source_candidate_id = id;
    B.cohort = repmat(cohort, n,1);
    if cohort == "internal", B.split = normSplit(M.split); else, B.split = repmat("external", n,1); end
    B.case_id = string(M.case_id);
    B.patient_id = patientId(B.case_id, cohort);
    B.candidate_label_value = double(M.candidate_label_value);
    B.candidate_labelmap_path = string(M.candidate_labelmap_path);
    B.GT_path = string(M.GT_path);
    f = double(M.candidate_overlap_fraction);
    if any(~isfinite(f)) || any(f < 0 | f > 1), error("Invalid overlap fraction in %s B.", cohort); end
    B.y = double(f > 0);
    B.f_overlap = f;
    B.overlap_volume_mL = double(M.overlap_volume_mL);
    B.label_any = double(f > 0);
    B.label_25pct = double(f >= 0.25);
    hp = NaN(n,1); hp(f == 0) = 0; hp(f >= 0.5) = 1;
    B.label_highpurity = hp;
    oc = repmat("C0", n,1);
    oc(f > 0 & f < 0.25) = "C1"; oc(f >= 0.25 & f < 0.5) = "C2"; oc(f >= 0.5) = "C3";
    B.occupancy_category = oc;
    B.Astd_excluded_5mm = false(n,1);
end

function reportBMismatch(M, B, cohort)
    a = double(M.candidate_label_any_overlap);
    fprintf("  %s B: any-overlap label vs f>0 disagreements: %d\n", cohort, sum(a ~= B.label_any));
    if any(a ~= B.label_any), error("%s B: stored any-overlap label disagrees with f > 0.", cohort); end
    v = string(M.Properties.VariableNames);
    if ismember("candidate_label_overlap_25pct", v)
        fprintf("  %s B: stored 25%% label vs f>=0.25 disagreements: %d (methods definition f>=0.25 is used)\n", ...
            cohort, sum(double(M.candidate_label_overlap_25pct) ~= B.label_25pct));
    end
    if ismember("candidate_label_overlap_50pct", v)
        fprintf("  %s B: stored 50%% label vs f>=0.50 disagreements: %d\n", ...
            cohort, sum(double(M.candidate_label_overlap_50pct) ~= double(B.f_overlap >= 0.5)));
    end
end

function checkSplitConsistency(B, A)
    [cb, ib] = unique(B.case_id); sb = B.split(ib);
    [ca, ia] = unique(A.case_id); sa = A.split(ia);
    [tf, loc] = ismember(ca, cb);
    bad = tf & sa ~= sb(max(loc,1));
    if any(bad), error("Split mismatch between B and A for %d studies (e.g. %s).", sum(bad), ca(find(bad,1))); end
    fprintf("  Split consistency B vs A: %d shared studies, 0 mismatches; A-only studies: %d; B-only studies: %d\n", ...
        sum(tf), sum(~tf), numel(setdiff(cb, ca)));
    pb = unique([B.patient_id, B.split], "rows"); pa = unique([A.patient_id, A.split], "rows");
    pp = unique([pb; pa], "rows");
    if numel(unique(pp(:,1))) ~= size(pp,1), error("A patient appears in more than one split."); end
end

function checkCount(label, got, expected)
    if got ~= expected
        error("CHECK FAILED: %s = %d (expected %d).", label, got, expected);
    end
    fprintf("  check OK: %-55s %d\n", label, got);
end

function [excluded, touchRaw, caseSum] = fiveMmFlags(B, marginMm)
% Pipeline B negatives with any voxel within marginMm (physical ellipsoid)
% of the reference mask; flagged whole. Mirrors negative_5mm_retention_audit.py.
    n = height(B);
    excluded = false(n,1);
    touchRaw = false(n,1);
    cases = unique(B.case_id, "stable");
    nc = numel(cases);
    cid = strings(nc,1); nNeg = zeros(nc,1); nExc = zeros(nc,1); spx = NaN(nc,3);
    for c = 1:nc
        idx = find(B.case_id == cases(c));
        neg = idx(B.y(idx) == 0);
        cid(c) = cases(c); nNeg(c) = numel(neg);
        if isempty(neg), continue; end
        lp = B.candidate_labelmap_path(idx(1)); gp = B.GT_path(idx(1));
        if ~isfile(lp), error("Label map not found: %s", lp); end
        if ~isfile(gp), error("GT not found: %s", gp); end
        iL = niftiinfo(lp); iG = niftiinfo(gp);
        L = niftiread(iL); G = niftiread(iG) > 0;
        sp = double(iL.PixelDimensions(1:3)); spG = double(iG.PixelDimensions(1:3));
        if ~isequal(size(L), size(G)) || any(abs(sp - spG) > 1e-3)
            error("Label map / GT geometry mismatch for %s", cases(c));
        end
        spx(c,:) = sp;
        hit = []; rawHit = [];
        if any(G(:))
            rad = floor(marginMm ./ sp);
            [ii,jj,kk] = ndgrid(-rad(1):rad(1), -rad(2):rad(2), -rad(3):rad(3));
            nh = (ii*sp(1)).^2 + (jj*sp(2)).^2 + (kk*sp(3)).^2 <= marginMm^2 + 1e-9;
            [x,y,z] = ind2sub(size(G), find(G));
            lo = max([min(x) min(y) min(z)] - rad - 1, 1);
            hi = min([max(x) max(y) max(z)] + rad + 1, size(G));
            Gs = G(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
            Ls = L(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
            D = imdilate(Gs, nh);
            hit = unique(double(Ls(D & Ls > 0)));
            rawHit = unique(double(L(G & L > 0)));
        end
        lv = B.candidate_label_value(neg);
        excluded(neg) = ismember(lv, hit);
        touchRaw(neg) = ismember(lv, rawHit);
        nExc(c) = sum(excluded(neg));
        if mod(c, 50) == 0 || c == nc
            fprintf("  5 mm flags: %d / %d studies, excluded so far %d\n", c, nc, sum(excluded));
        end
    end
    caseSum = table(cid, nNeg, nExc, spx(:,1), spx(:,2), spx(:,3), ...
        'VariableNames', {'case_id','n_neg','n_neg_excluded_5mm','sx','sy','sz'});
end
