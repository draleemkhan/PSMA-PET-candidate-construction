# Reference-informed versus PET-only candidate construction in PSMA PET/CT — analysis code

Code for the study *"Reference-informed versus PET-only candidate construction in PSMA PET/CT: effects on classification, tumour segmentation and burden estimation"* (manuscript under review; citation to be added on publication).

The repository contains **code only**. It contains no images, masks, feature tables, model files, results or logs. The `split/` folder holds two identifier lists from the public datasets, which are needed to reproduce the patient-level split and the locked external cohort.

## Data (download separately)

| Cohort | Dataset | Access |
|---|---|---|
| Internal (development and internal test) | PSMA-PET-CT-Lesions, LMU Munich — 597 [¹⁸F]PSMA-1007 / [⁶⁸Ga]PSMA-11 PET/CT studies, 378 patients | The Cancer Imaging Archive, doi:10.7937/R7EP-3X37 (release cited in the paper); data descriptor: Jeblick et al., *Scientific Data* 13, 1023 (2026), doi:10.1038/s41597-026-07821-z |
| External | DEEP-PSMA (Peter MacCallum Cancer Centre) | Zenodo, doi:10.5281/zenodo.15281784 (CC BY-NC 4.0) |

Expected layout: `PSMA_INTERNAL_DATA` must contain `imagesTr/` (`<case>_0000.nii.gz` CT, `<case>_0001.nii.gz` PET in SUV), `labelsTr/` (`<case>.nii.gz`) and `psma_metadata.csv`. `PSMA_EXTERNAL_DATA` is the extracted DEEP-PSMA folder with one sub-folder per case (`train_0001` … `train_0100`).

## Software

- MATLAB R2025a (25.1.0.2943329) with Image Processing Toolbox and Statistics and Machine Learning Toolbox.
- Python 3.10.11 with the packages in `requirements.txt` (PyRadiomics 3.0.1, SimpleITK 2.5.5, NumPy 1.26.4, pandas 2.3.3).
- Windows with PowerShell. Two Python steps (`…9C_…` routine features) map a free drive letter (`O:`, or `Q:` for the addendum) with `subst` to avoid long-path limits; these letters must be free.

## Setup

1. Copy `config/set_psma_paths_TEMPLATE.m` to `config/set_psma_paths.m` and `config/set_psma_paths_TEMPLATE.ps1` to `config/set_psma_paths.ps1`; set the three folders:
   - `PSMA_PROJECT_ROOT` — empty output folder. The scripts create `PSMA_AI_Project/` (stage 1 outputs) and `PSMA_AI_Project_V2_PET_Derived/` (stage 2–3 outputs) inside it.
   - `PSMA_INTERNAL_DATA`, `PSMA_EXTERNAL_DATA` — the two datasets.
2. In MATLAB run `config/set_psma_paths.m` once per session (it also adds `config/` to the path). In PowerShell dot-source `config/set_psma_paths.ps1` before running the Python steps.
3. Radiomics settings are in `config/PSMA_objective5C_PET_radiomics_params.yaml` and `config/PSMA_objective5C_CT_radiomics_params.yaml`; all extraction steps use these two files.

## Run order

Script names are kept as used in the study so that they match the analysis registers.

### Stage 1 — reference-informed candidates, split and external cohort (`01_V1_reference_informed_candidates_and_split/`)

| # | Script | Purpose |
|---|---|---|
| 1 | `PSMA_00A_build_manifest_from_public_dataset.m` | Case manifest from the downloaded dataset (written for this release) |
| 2 | `PSMA_objective4_create_TRUE_vs_FALSE_candidate_table.m` | Historical reference-informed candidates (A-hist): reference components ≥ 1 mL; SUV ≥ 3 uptake outside the reference dilated by one voxel |
| 3 | `PSMA_objective4B_patient_based_train_val_test_split.m` | Patient-level 70/15/15 split, seed 2026 (compare with `split/`) |
| 4 | `PSMA_00B_create_physical_train_val_test_folders.m` | Split folders read by stage 2 |
| 5 | `PSMA_objective5C_01_create_candidate_labelmaps_for_radiomics.m` | Label maps of reference-derived candidates |
| 6 | `PSMA_objective5C_02_run_pyradiomics_FIXED_VERBOSE.py` | PET radiomics of reference-derived internal candidates |
| 7 | `PSMA_objective7A2_…`, `PSMA_objective7A3_…` | DEEP-PSMA world-space alignment and CT-coverage audits |
| 8 | `PSMA_objective7A4_lock_external_validation_cohort_95.m` | Locked external cohort (95 of 100; compare with `split/`) |
| 9 | `PSMA_objective7B1_…_CORRECTED.m` | External reference-mask analysis |
| 10 | `PSMA_objective8B_…` | External reference-derived candidates |
| 11 | `Run-PSMA-Objective8D8E-SEPARATE-YAML.ps1` | PET radiomics of external reference-derived candidates |

### Stage 2 — PET-only (Pipeline B) candidates (`02_V2_PET_only_candidate_pipeline/`)

| # | Script | Purpose |
|---|---|---|
| 1 | `PSMA_objective9A_generate_PET_only_candidates_CORRECTED.m`, `PSMA_objective9A_external_PSMA_only_candidates.m` | SUV ≥ 3 connected components ≥ 1 mL (26-connectivity), reference not read |
| 2 | `PSMA_objective9B_measure_GT_overlap_and_assign_labels_DELIMITER_FIXED.m`, `PSMA_objective9B_external_PSMA_overlap_and_labels.m` | Reference overlap, occupancy and labels |
| 3 | `Run-PSMA-Objective9C-FINAL.ps1` | Routine features |
| 4 | `Run-PSMA-Objective9D9E-Direct-Radiomics.ps1` (internal), `Run-PSMA-Objective9-External-PET-Radiomics.ps1`, `Run-PSMA-Objective9-External-CT-Radiomics.ps1` | PET and CT radiomics |
| 5 | `PSMA_objective9F_build_internal_MASTER_candidate_dataset_REVISED_v5.m` | Internal master candidate table |
| 6 | `ADDENDUM5_external_train0096_0100/` | Same steps for external cases train_0096–0100 (documented addendum run, so that the evaluated set equals the locked 95-case cohort); includes a check of the LogitBoost probability transform |

### Stage 3 — study analysis (`03_LeakagePaper_analysis/`)

| # | Script | Purpose |
|---|---|---|
| 1 | `LP01_extract_CT_radiomics_reference_derived.py` | CT radiomics of reference-derived candidates (current rule) and PET spot check |
| 2 | `LP02_extract_routine_features_reference_derived.py` | Routine features of reference-derived candidates |
| 3 | `LP03_assemble_candidate_sets.m` | Pipelines B, A-std and A-hist |
| 4 | `LP04_train_models_and_thresholds.m` | LogitBoost models and locked validation thresholds |
| 5 | `LP05a_candidate_level_evaluation.m` | Primary and secondary candidate-level endpoints, bootstrap |
| 6 | `LP05aD_mechanism_diagnostic.m`, `LP05aD2_SUVmin_band_diagnostic.m` | Post hoc diagnostics |
| 7 | `LP05b_segmentation_burden_detection.m` | Segmentation, burden and lesion detection |
| 8 | `LP06_main_display_items.m`, `LP06_main_Fig1_Fig4_revision2.m`, `LP06_supplementary_items.m` | Tables and figures |
| 9 | `LP06_Fig5_illustrative_case.m` | Figure 5 (representative case; the study identifier is in the script and in `split/`) |

## Notes on re-running

- **Frozen run identifiers.** The study froze each stage. Stage 3 scripts therefore read earlier outputs by their run folder or time-stamped file name (for example `LP03_20261109_203442`, `PSMA_objective9B_candidate_overlap_table_20260716_090139.csv`). After a new run, edit these names near the top of each script to point to your own outputs.
- **Built-in checks.** Several scripts stop if counts differ from the published ones (for example 2,673 internal test candidates) or if input SHA256 hashes differ from the recorded ones. Count checks are intended as reproduction checks; hash checks identify our frozen files and will differ for a new run.
- **Split.** `split/internal_PSMA-PET-CT-Lesions_patient_split.csv` (study_id, patient_id, split) is the split used in the paper; `split/external_DEEP-PSMA_locked95_eligibility.csv` lists the 100 DEEP-PSMA cases with the locked inclusion flag.

## Licence

PolyForm Noncommercial License 1.0.0 — free for research and other noncommercial use; commercial use is not permitted. See `LICENSE.md`. The datasets have their own licences.

## Citation

See `CITATION.cff`. Please also cite the two datasets.
