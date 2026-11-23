#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Leakage paper - LP02: routine features for reference-derived candidates
=======================================================================
Locked methods: LeakagePaper_Methods_LOCKED_v1.0_2026-11-09.md, Section 5.1 and 5.3.

Purpose
-------
Compute routine candidate features with the V2 Objective 9C code for all reference-derived
candidates (Pipeline A positives = reference lesions >= 1 mL; A-hist negatives = V1 carved
components), on the existing V1 masks:

  --set internal : V1 Objective 5C label maps, all 597 studies (expected 22,794 candidates:
                   11,380 positives + 11,414 A-hist negatives)
  --set external : V1 Objective 8B label maps, locked 95 cohort (expected 4,658 candidates:
                   3,250 positives + 1,408 A-hist negatives)

Identity with V2
----------------
The per-case computation is NOT re-implemented. The V2 script
PSMA_objective9C_extract_routine_features_FINAL.py is loaded and its extract_case() is called
unchanged (SUV/CT statistics with SD N-1, linear CT-to-PET resampling, centroid and bounding-box
definitions, location zones, strict label-set check, voxel/volume reconciliation).
One line is changed IN MEMORY ONLY when loading it: the module-level call that maps a short
"O:" drive alias (ensure_short_drive) is replaced by the long project path, so this script does
not create or change any drive mapping. The V2 file itself is not modified. All inputs used here
lie outside the folder that alias refers to, so the alias is never needed.

Safety
------
Reads only. Writes only into a NEW folder:
  ...\PSMA_AI_Project_V2_PET_Derived\06_results\LeakagePaper_FeatureExtraction\routine_<set>\
Resumable: completed studies are skipped on rerun.

Run (PowerShell, same Python environment as the V2 pipeline):
  python LP02_extract_routine_features_reference_derived.py --set internal --workers 4
  python LP02_extract_routine_features_reference_derived.py --set external --workers 4
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import sys
import time
import traceback
import types
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os

import pandas as pd

# -------------------------------------------------------------------------
# Locked paths
# -------------------------------------------------------------------------
PROJECT = Path(os.environ["PSMA_PROJECT_ROOT"])
V1 = PROJECT / "PSMA_AI_Project"
V2 = PROJECT / "PSMA_AI_Project_V2_PET_Derived"
V2_9C_FILE = Path(__file__).resolve().parents[1] / "02_V2_PET_only_candidate_pipeline" / "PSMA_objective9C_extract_routine_features_FINAL.py"
DEEP_PSMA = Path(os.environ["PSMA_EXTERNAL_DATA"])

V1_5C_DIR = V1 / "06_results" / "Objective5C FINAL-FINAL-FINALL routine plus_Radiomics"
V1_5C_INDEX = V1_5C_DIR / "PSMA_objective5C_candidate_labelmap_index.csv"
V1_5C_MASK_DIR = V1_5C_DIR / "candidate_masks"
V1_8B_TABLE = (V1 / "06_results" / "Objective8_Corrected_External_Validation"
               / "Objective8B_Exact_InternalStyle_Candidate_Dataset"
               / "PSMA_objective8B_exact_internal_style_candidates_20260714_125116.csv")
LOCK_TABLE = (V1 / "06_results" / "Objective7_DEEP_PSMA_External_Validation"
              / "Objective7A4_Locked_Eligibility"
              / "PSMA_objective7_external_validation_eligibility_LOCKED_95.csv")
V2_9B_INTERNAL = (V2 / "06_results" / "Objective9_PET_Derived_Candidate_Pipeline"
                  / "Objective9B_GT_Overlap_and_Label_Assignment" / "internal")
OUT_ROOT = V2 / "06_results" / "LeakagePaper_FeatureExtraction"

EXPECTED = {"internal": (597, 22794, 11380, 11414), "external": (95, 4658, 3250, 1408)}

# -------------------------------------------------------------------------
# Load the V2 9C module unchanged except the drive-alias line (in memory only)
# -------------------------------------------------------------------------
_MOD_NAME = "psma_v2_9c"
_src = V2_9C_FILE.read_text(encoding="utf-8")
_OLD = "OBJ9_ROOT = ensure_short_drive()"
if _src.count(_OLD) != 1:
    raise RuntimeError("Unexpected V2 9C source: drive-alias line not found exactly once")
_src = _src.replace(_OLD, "OBJ9_ROOT = LONG_OBJ9_ROOT  # LP02: no drive alias (in-memory change)")
V2_9C = types.ModuleType(_MOD_NAME)
V2_9C.__file__ = str(V2_9C_FILE)
sys.modules[_MOD_NAME] = V2_9C
exec(compile(_src, str(V2_9C_FILE), "exec"), V2_9C.__dict__)


# -------------------------------------------------------------------------
# Candidate tables
# -------------------------------------------------------------------------
def internal_table() -> pd.DataFrame:
    x = pd.read_csv(V1_5C_INDEX, low_memory=False)
    files = list(V2_9B_INTERNAL.glob("*/checkpoints/PSMA_objective9B_candidate_overlap_checkpoint.csv"))
    if len(files) != 3:
        raise RuntimeError(f"Expected 3 internal 9B checkpoints, found {len(files)}")
    pet = pd.concat([pd.read_csv(f, usecols=["case_id", "PET_path", "split"]) for f in files]).drop_duplicates("case_id").set_index("case_id")
    if x.case_id.map(pet.PET_path).isna().any():
        raise RuntimeError("Some V1 internal studies have no V2 PET path")
    return pd.DataFrame({
        "candidate_global_id": x.candidate_global_id.astype(str),
        "dataset_name": "PSMA-PET-CT-Lesions",
        "cohort": "internal",
        "split": x.case_id.map(pet.split).astype(str),
        "case_id": x.case_id.astype(str),
        "candidate_label_value": x.candidate_label_value.astype(int),
        "candidate_voxels": x.candidate_voxels.astype(int),
        "candidate_volume_mL": x.candidate_volume_mL.astype(float),
        "PET_path": x.case_id.map(pet.PET_path).astype(str),
        "candidate_labelmap_path": [str(V1_5C_MASK_DIR / Path(str(p).replace("\\", "/")).name) for p in x.candidate_labelmap_path],
        "v1_candidate_label_binary": x.candidate_label_binary.astype(int),
        "source_table": "V1_Objective5C_labelmap_index",
    })


def external_table() -> pd.DataFrame:
    lock = pd.read_csv(LOCK_TABLE)
    inc = set(lock.loc[lock.include_external_validation.astype(str).str.lower().isin(["1", "true"]), "case_id"])
    b = pd.read_csv(V1_8B_TABLE, low_memory=False)
    if set(b.case_id) != inc:
        raise RuntimeError("V1 8B case set differs from the locked 95-case cohort")
    return pd.DataFrame({
        "candidate_global_id": b.candidate_global_id.astype(str),
        "dataset_name": "DEEP-PSMA",
        "cohort": "external",
        "split": "external",
        "case_id": b.case_id.astype(str),
        "candidate_label_value": b.candidate_label_value.astype(int),
        "candidate_voxels": b.candidate_voxels.astype(int),
        "candidate_volume_mL": b.candidate_volume_mL.astype(float),
        "PET_path": [str(DEEP_PSMA / c / "PSMA" / "PET.nii.gz") for c in b.case_id],
        "candidate_labelmap_path": b.candidate_label_map_path.astype(str),
        "v1_candidate_label_binary": b.candidate_label_binary.astype(int),
        "source_table": "V1_Objective8B",
    })


def validate(t: pd.DataFrame, which: str) -> None:
    n_studies, n_rows, n_pos, n_neg = EXPECTED[which]
    got = (t.case_id.nunique(), len(t), int((t.v1_candidate_label_binary == 1).sum()),
           int((t.v1_candidate_label_binary == 0).sum()))
    if got != (n_studies, n_rows, n_pos, n_neg):
        raise RuntimeError(f"{which}: expected (studies, rows, pos, neg) = {(n_studies, n_rows, n_pos, n_neg)}, got {got}")
    if t.candidate_global_id.duplicated().any():
        raise RuntimeError(f"{which}: duplicate candidate IDs")
    if t.duplicated(["case_id", "candidate_label_value"]).any():
        raise RuntimeError(f"{which}: duplicate label values within a study")
    if (t.groupby("case_id").candidate_labelmap_path.nunique() != 1).any():
        raise RuntimeError(f"{which}: more than one label map per study")
    missing = [p for p in set(t.candidate_labelmap_path) | set(t.PET_path) if not Path(p).is_file()]
    if missing:
        raise RuntimeError(f"{which}: {len(missing)} input files missing, e.g. {missing[:3]}")


# -------------------------------------------------------------------------
# Worker (one study)
# -------------------------------------------------------------------------
def run_case(case_records: list, scope: str, split: str):
    case_df = pd.DataFrame(case_records)
    case_id = str(case_df.case_id.iloc[0])
    t0 = time.time()
    try:
        rows = V2_9C.extract_case(case_df, scope, split)
        extra = case_df.set_index("candidate_global_id")
        for r in rows:
            cid = r["candidate_global_id"]
            r["v1_candidate_label_binary"] = int(extra.loc[cid, "v1_candidate_label_binary"])
            r["source_table"] = str(extra.loc[cid, "source_table"])
        return case_id, rows, [], time.time() - t0
    except Exception as exc:
        err = [{"case_id": case_id, "n_candidates": len(case_df), "error_message": str(exc),
                "traceback": traceback.format_exc(limit=4)}]
        return case_id, [], err, time.time() - t0


def versions() -> dict:
    v = {"python": sys.version.split()[0], "platform": platform.platform(), "V2_9C": str(V2_9C_FILE)}
    for name in ["SimpleITK", "numpy", "pandas"]:
        try:
            v[name] = __import__(name).__version__
        except Exception as exc:
            v[name] = f"unavailable ({exc})"
    return v


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["internal", "external"])
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    which = args.set

    out_dir = OUT_ROOT / f"routine_{which}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"LP02_routine_{which}_checkpoint.csv"
    fail_csv = out_dir / f"LP02_routine_{which}_failures.csv"
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(out_dir / f"LP02_routine_{which}.log", encoding="utf-8")],
                        force=True)
    ver = versions()
    logging.info("LP02 routine features | set %s | versions %s", which, json.dumps(ver))

    t = internal_table() if which == "internal" else external_table()
    validate(t, which)
    manifest = out_dir / f"LP02_routine_{which}_input_manifest.csv"
    if not manifest.exists():
        t.to_csv(manifest, index=False)
    logging.info("Candidates: %d in %d studies (positives %d, A-hist negatives %d)", len(t), t.case_id.nunique(),
                 int((t.v1_candidate_label_binary == 1).sum()), int((t.v1_candidate_label_binary == 0).sum()))

    done = set()
    if out_csv.exists():
        done = set(pd.read_csv(out_csv, usecols=["case_id"]).case_id.astype(str))
    todo = t[~t.case_id.isin(done)]
    logging.info("Studies already completed: %d | to process now: %d", len(done), todo.case_id.nunique())

    cols = V2_9C.OUTPUT_COLUMNS + ["v1_candidate_label_binary", "source_table"]
    t0 = time.time(); n_ok = n_err = 0
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = []
        for cid, g in todo.groupby("case_id", sort=True):
            futs.append(ex.submit(run_case, g.to_dict("records"), which, str(g.split.iloc[0])))
        for i, fu in enumerate(as_completed(futs), 1):
            case_id, rows, err, sec = fu.result()
            if rows:
                pd.DataFrame(rows)[cols].to_csv(out_csv, mode="a", header=not out_csv.exists(), index=False)
            if err:
                pd.DataFrame(err).to_csv(fail_csv, mode="a", header=not fail_csv.exists(), index=False)
            n_ok += len(rows); n_err += len(err)
            logging.info("[%d/%d] %s | rows=%d failed=%s | %.1fs", i, len(futs), case_id, len(rows), bool(err), sec)

    total = len(pd.read_csv(out_csv, usecols=["candidate_global_id"])) if out_csv.exists() else 0
    state = {"set": which, "expected_candidates": len(t), "success_total": total, "new_rows": n_ok,
             "new_failed_studies": n_err, "minutes": round((time.time() - t0) / 60, 2), "versions": ver,
             "updated": time.strftime("%Y-%m-%d %H:%M:%S")}
    (out_dir / f"LP02_routine_{which}_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")
    logging.info("DONE: %d / %d candidates | failed studies this run: %d", total, len(t), n_err)
    return 0 if total == len(t) else 1


if __name__ == "__main__":
    sys.exit(main())
