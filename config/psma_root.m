function p = psma_root(name)
%PSMA_ROOT Return a root folder defined by an environment variable.
%   p = psma_root("PSMA_PROJECT_ROOT")   output root; scripts create and use
%                                        <root>\PSMA_AI_Project and
%                                        <root>\PSMA_AI_Project_V2_PET_Derived
%   p = psma_root("PSMA_INTERNAL_DATA")  PSMA-PET-CT-Lesions root
%                                        (imagesTr, labelsTr, psma_metadata.csv)
%   p = psma_root("PSMA_EXTERNAL_DATA")  extracted DEEP-PSMA root
%
%   Set the variables once per MATLAB session with config/set_psma_paths.m
%   (copy set_psma_paths_TEMPLATE.m and edit it).
p = string(getenv(name));
if strlength(p) == 0
    error("psma_root:unset", ...
        "Environment variable %s is not set. Run config/set_psma_paths.m first.", name);
end
if ~isfolder(p)
    error("psma_root:missing", "%s points to a folder that does not exist: %s", name, p);
end
end
