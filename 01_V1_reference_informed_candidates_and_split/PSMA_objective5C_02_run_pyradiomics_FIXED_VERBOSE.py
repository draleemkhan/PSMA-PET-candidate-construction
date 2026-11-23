"""
Objective 5C Step 02: Extract PET and CT radiomics for PSMA candidate labelmaps.

Run first in pilot mode. After successful pilot, change RUN_MODE = "full".

Required packages:
    pip install pyradiomics pandas openpyxl pyyaml tqdm

Input batch CSV files are created by:
    PSMA_objective5C_01_create_candidate_labelmaps_for_radiomics.m
"""

from __future__ import annotations

import os
import sys
import time
import logging
from pathlib import Path
from typing import Dict, Any, List

import pandas as pd
from tqdm import tqdm
from radiomics import featureextractor


# ============================ USER SETTINGS =============================
PROJECT_ROOT = (Path(os.environ["PSMA_PROJECT_ROOT"]) / "PSMA_AI_Project")
RAD_ROOT = PROJECT_ROOT / "06_results" / "Objective5C_Radiomics"
FEATURE_DIR = RAD_ROOT / "radiomics_features"
LOG_DIR = RAD_ROOT / "logs"

PET_BATCH_CSV = FEATURE_DIR / "pyradiomics_batch_PET_candidates.csv"
CT_BATCH_CSV = FEATURE_DIR / "pyradiomics_batch_CT_candidates.csv"

PET_PARAM_YAML = Path(__file__).resolve().parents[1] / "config" / "PSMA_objective5C_PET_radiomics_params.yaml"
CT_PARAM_YAML = Path(__file__).resolve().parents[1] / "config" / "PSMA_objective5C_CT_radiomics_params.yaml"

# Start safely with pilot. After it succeeds, change to "full".
RUN_MODE = "full"      # "pilot" or "full"
PILOT_N_ROWS = 100       # used only when RUN_MODE = "pilot"

# If extraction stops midway, set RESUME = True to skip rows already present in output CSV.
RESUME = True

# Save every N successful/error rows so progress is not lost.
SAVE_EVERY = 25
# ========================================================================


def setup_logging() -> None:
    print("Starting Objective 5C Step 02 script...", flush=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / f"PSMA_objective5C_02_pyradiomics_{RUN_MODE}.log"

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )
    print(f"Log file: {log_file}", flush=True)

    # Reduce PyRadiomics internal verbosity.
    logging.getLogger("radiomics").setLevel(logging.ERROR)


def validate_inputs() -> None:
    required = [PET_BATCH_CSV, CT_BATCH_CSV, PET_PARAM_YAML, CT_PARAM_YAML]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required input files:\n" + "\n".join(missing))


def read_batch(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required_cols = {"Image", "Mask", "Label", "case_id", "candidate_global_id"}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {sorted(missing)}")

    if RUN_MODE.lower() == "pilot":
        df = df.head(PILOT_N_ROWS).copy()

    df["Label"] = df["Label"].astype(int)
    df["candidate_global_id"] = df["candidate_global_id"].astype(int)
    return df


def make_extractor(param_yaml: Path) -> featureextractor.RadiomicsFeatureExtractor:
    ext = featureextractor.RadiomicsFeatureExtractor(str(param_yaml))
    return ext


def clean_feature_dict(raw: Dict[str, Any], prefix: str) -> Dict[str, Any]:
    """Keep radiomic features and remove PyRadiomics diagnostics fields."""
    out: Dict[str, Any] = {}
    for k, v in raw.items():
        key = str(k)
        if key.startswith("diagnostics_"):
            continue
        # Convert numpy scalar-like objects when possible.
        try:
            if hasattr(v, "item"):
                v = v.item()
        except Exception:
            pass
        out[f"{prefix}_{key}"] = v
    return out


def load_existing_ids(out_csv: Path) -> set[int]:
    if not out_csv.exists() or not RESUME:
        return set()
    try:
        existing = pd.read_csv(out_csv, usecols=["candidate_global_id"])
        return set(existing["candidate_global_id"].dropna().astype(int).tolist())
    except Exception:
        return set()


def append_csv(rows: List[Dict[str, Any]], out_csv: Path) -> None:
    if not rows:
        return
    df = pd.DataFrame(rows)
    write_header = not out_csv.exists()
    df.to_csv(out_csv, mode="a", header=write_header, index=False)


def run_one_modality(
    modality_name: str,
    batch_csv: Path,
    param_yaml: Path,
    prefix: str,
    out_csv: Path,
    error_csv: Path,
) -> None:
    logging.info("============================================================")
    logging.info("Starting %s radiomics", modality_name)
    logging.info("Batch CSV: %s", batch_csv)
    logging.info("Params:    %s", param_yaml)
    logging.info("Output:    %s", out_csv)

    df = read_batch(batch_csv)
    extractor = make_extractor(param_yaml)

    done_ids = load_existing_ids(out_csv)
    if done_ids:
        logging.info("Resume enabled: %d candidates already in %s", len(done_ids), out_csv.name)

    rows_buffer: List[Dict[str, Any]] = []
    errors_buffer: List[Dict[str, Any]] = []
    n_success = 0
    n_error = 0
    n_skipped = 0

    t0 = time.time()

    for _, row in tqdm(df.iterrows(), total=len(df), desc=f"{modality_name} radiomics"):
        cid = int(row["candidate_global_id"])
        if cid in done_ids:
            n_skipped += 1
            continue

        image_path = str(row["Image"])
        mask_path = str(row["Mask"])
        label = int(row["Label"])

        base = {
            "candidate_global_id": cid,
            "case_id": row.get("case_id", ""),
            "candidate_label_binary": row.get("candidate_label_binary", ""),
            "candidate_source": row.get("candidate_source", ""),
            "split_set": row.get("split_set", ""),
            "label_value_in_mask": label,
            "Image": image_path,
            "Mask": mask_path,
        }

        try:
            if not os.path.exists(image_path):
                raise FileNotFoundError(f"Image not found: {image_path}")
            if not os.path.exists(mask_path):
                raise FileNotFoundError(f"Mask not found: {mask_path}")

            result = extractor.execute(image_path, mask_path, label=label)
            feats = clean_feature_dict(result, prefix=prefix)
            rows_buffer.append({**base, **feats})
            n_success += 1

        except Exception as e:
            errors_buffer.append({**base, "error_message": str(e)})
            n_error += 1

        if (n_success + n_error) % SAVE_EVERY == 0:
            append_csv(rows_buffer, out_csv)
            append_csv(errors_buffer, error_csv)
            rows_buffer.clear()
            errors_buffer.clear()

    append_csv(rows_buffer, out_csv)
    append_csv(errors_buffer, error_csv)

    elapsed_min = (time.time() - t0) / 60.0
    logging.info("DONE %s radiomics", modality_name)
    logging.info("Rows requested: %d", len(df))
    logging.info("Success:        %d", n_success)
    logging.info("Errors:         %d", n_error)
    logging.info("Skipped:        %d", n_skipped)
    logging.info("Elapsed min:    %.2f", elapsed_min)
    logging.info("Output CSV:     %s", out_csv)
    logging.info("Error CSV:      %s", error_csv)


def main() -> None:
    setup_logging()
    FEATURE_DIR.mkdir(parents=True, exist_ok=True)
    validate_inputs()

    suffix = "PILOT" if RUN_MODE.lower() == "pilot" else "FULL"

    pet_out = FEATURE_DIR / f"PSMA_objective5C_PET_radiomics_{suffix}.csv"
    pet_err = FEATURE_DIR / f"PSMA_objective5C_PET_radiomics_errors_{suffix}.csv"
    ct_out = FEATURE_DIR / f"PSMA_objective5C_CT_radiomics_{suffix}.csv"
    ct_err = FEATURE_DIR / f"PSMA_objective5C_CT_radiomics_errors_{suffix}.csv"

    print("Objective 5C Step 02 PyRadiomics extraction", flush=True)
    logging.info("Objective 5C Step 02 PyRadiomics extraction")
    logging.info("RUN_MODE: %s", RUN_MODE)
    logging.info("PROJECT_ROOT: %s", PROJECT_ROOT)

    print("Starting PET radiomics...", flush=True)
    run_one_modality("PET", PET_BATCH_CSV, PET_PARAM_YAML, "PET", pet_out, pet_err)
    print("Starting CT radiomics...", flush=True)
    run_one_modality("CT", CT_BATCH_CSV, CT_PARAM_YAML, "CT", ct_out, ct_err)

    print("DONE Objective 5C Step 02.", flush=True)
    logging.info("============================================================")
    logging.info("DONE Objective 5C Step 02.")
    logging.info("Next step: merge PET+CT radiomics with Objective 4 routine candidate table.")
    logging.info("============================================================")


if __name__ == "__main__":
    main()
