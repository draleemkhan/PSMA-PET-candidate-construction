%% PSMA_00A_build_manifest_from_public_dataset.m
% Build the internal case manifest from the public PSMA-PET-CT-Lesions dataset
% (TCIA, doi:10.7937/R7EP-3X37, version used: v3 file layout).
%
% Input  (PSMA_INTERNAL_DATA):
%   imagesTr\<case_id>_0000.nii.gz   CT
%   imagesTr\<case_id>_0001.nii.gz   PET (SUV)
%   labelsTr\<case_id>.nii.gz        reference tumour mask
%   psma_metadata.csv                Subject ID, Study Date, age, scanner,
%                                    radionuclide, CT contrast
% Output (PSMA_PROJECT_ROOT\PSMA_AI_Project\01_manifest):
%   PSMA_final_manifest.xlsx / .csv            read by PSMA_objective4_...
%   PSMA_final_manifest_D_DRIVE.xlsx / .csv    read by PSMA_objective5C_01_...
%   (identical content; both names are kept because the downstream scripts
%    were written for them)
%
% case_id has the form psma_<subject hash>_<YYYY-MM-DD>; the subject key is
% "PSMA_<hash>" and the study date is the date part. Rows are sorted by case_id.
% Run config/set_psma_paths.m first.

clear; clc;

dataRoot = psma_root("PSMA_INTERNAL_DATA");
outDir   = fullfile(psma_root("PSMA_PROJECT_ROOT"), "PSMA_AI_Project", "01_manifest");
if ~isfolder(outDir), mkdir(outDir); end

imgDir = fullfile(dataRoot, "imagesTr");
labDir = fullfile(dataRoot, "labelsTr");
metaFile = fullfile(dataRoot, "psma_metadata.csv");
assert(isfolder(imgDir) && isfolder(labDir) && isfile(metaFile), ...
    "PSMA_INTERNAL_DATA must contain imagesTr, labelsTr and psma_metadata.csv.");

% ---- cases from CT files ------------------------------------------------
d = dir(fullfile(imgDir, "*_0000.nii.gz"));
case_id = sort(string(erase({d.name}', "_0000.nii.gz")));
n = numel(case_id);
fprintf("Cases found: %d\n", n);

CT_path    = fullfile(imgDir, case_id + "_0000.nii.gz");
PET_path   = fullfile(imgDir, case_id + "_0001.nii.gz");
Label_path = fullfile(labDir, case_id + ".nii.gz");
CT_exists    = double(arrayfun(@isfile, CT_path));
PET_exists   = double(arrayfun(@isfile, PET_path));
Label_exists = double(arrayfun(@isfile, Label_path));
has_CT = CT_exists; has_PET = PET_exists; has_Label = Label_exists;
all_images_exist = double(CT_exists & PET_exists & Label_exists);

tok = regexp(case_id, "^psma_([0-9a-f]+)_(\d{4}-\d{2}-\d{2})$", "tokens", "once");
bad = cellfun(@isempty, tok);
if any(bad)
    error("Unexpected case_id format: %s", strjoin(case_id(bad), ", "));
end
tok = vertcat(tok{:});
SubjectID_key = "psma_" + string(tok(:,1));
StudyDate_key = string(tok(:,2));
SubjectID     = "PSMA_" + string(tok(:,1));
StudyDate     = StudyDate_key;

% ---- metadata join --------------------------------------------------------
opts = detectImportOptions(metaFile, "VariableNamingRule", "preserve");
opts = setvartype(opts, opts.VariableNames, "string");
Mt = readtable(metaFile, opts);
Mt.SubjectID_key = lower(strtrim(Mt.("Subject ID")));
Mt.StudyDate_key = string(datetime(strtrim(Mt.("Study Date")), ...
    "InputFormat", "yyyy-MM-dd"), "yyyy-MM-dd");
key  = SubjectID_key + "|" + StudyDate_key;
mkey = Mt.SubjectID_key + "|" + Mt.StudyDate_key;
assert(numel(unique(mkey)) == numel(mkey), "Duplicate subject/date rows in psma_metadata.csv.");
[tf, loc] = ismember(key, mkey);
has_metadata = double(tf);

age = nan(n,1);
manufacturer_model_name = strings(n,1);
pet_radionuclide = strings(n,1);
ct_contrast_agent = strings(n,1);
age(tf) = str2double(Mt.age(loc(tf)));
manufacturer_model_name(tf) = Mt.manufacturer_model_name(loc(tf));
pet_radionuclide(tf) = Mt.pet_radionuclide(loc(tf));
ct_contrast_agent(tf) = Mt.ct_contrast_agent(loc(tf));

T = table(case_id, CT_path, PET_path, Label_path, has_CT, has_PET, has_Label, ...
    SubjectID, StudyDate, age, manufacturer_model_name, pet_radionuclide, ...
    ct_contrast_agent, SubjectID_key, StudyDate_key, has_metadata, ...
    CT_exists, PET_exists, Label_exists, all_images_exist);

fprintf("All images present: %d of %d\n", sum(all_images_exist), n);
fprintf("Metadata matched:   %d of %d\n", sum(has_metadata), n);
if sum(all_images_exist) ~= n || sum(has_metadata) ~= n
    warning("Some cases lack images or metadata; check before continuing.");
end
if n ~= 597
    warning("Expected 597 studies (as in the paper); found %d.", n);
end

for base = ["PSMA_final_manifest", "PSMA_final_manifest_D_DRIVE"]
    writetable(T, fullfile(outDir, base + ".csv"));
    xl = fullfile(outDir, base + ".xlsx");
    if isfile(xl), delete(xl); end
    writetable(T, xl);
end
fprintf("Manifest written to:\n%s\n", outDir);
