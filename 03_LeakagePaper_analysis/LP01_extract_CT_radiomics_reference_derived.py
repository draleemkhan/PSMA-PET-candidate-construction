#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
Leakage paper - LP01: CT radiomics (and PET spot-check) for reference-derived candidates
==========================================================================================
Locked methods: LeakagePaper_Methods_LOCKED_v1.0_2026-11-09.md, Section 5.2-5.3.

Purpose
-------
Extract CT radiomics with the V2 Objective 9D/9E rule for reference-derived candidates
(Pipeline A positives = reference lesions >= 1 mL, and A-hist negatives = V1 carved
components) where V1 CT values cannot be reused:

  --set internal_ct  : internal studies in which V2 resamples CT to the PET grid (107 studies;
                       expected 662 positives + 1,308 A-hist negatives = 1,970 regions)
  --set external_ct  : external locked 95 cohort (expected 3,250 positives + 1,408 A-hist
                       negatives = 4,658 regions)
  --set pet_spotcheck: PET radiomics for 20 internal reference-derived positives (random,
                       seed 20261109), compared with the reused V1 PET values

Identity with V2
----------------
The per-case extraction is NOT re-implemented. This script imports the V2 module
PSMA_objective9D9E_extract_radiomics_DIRECT_RESUMABLE_LONG_PATH_FIXED.py and calls its
process_case() unchanged: same PyRadiomics YAML files, same CT reading, same geometry
tolerance (1e-5), same B-spline CT resampling (default -1024 HU), same 3-voxel crop,
same label handling and same feature names (CT_original_..., PET_original_...).
Masks are the existing V1 label maps (internal: Objective 5C; external: Objective 8B).

Safety
------
Reads only. Writes only into a NEW output folder:
  ...\PSMA_AI_Project_V2_PET_Derived\06_results\LeakagePaper_FeatureExtraction\<set>\
No existing file is modified. Resumable: completed candidate IDs are skipped on rerun.

Run (PowerShell, same Python environment as the V2 radiomics):
  python LP01_extract_CT_radiomics_reference_derived.py --set internal_ct --workers 4
  python LP01_extract_CT_radiomics_reference_derived.py --set external_ct --workers 4
  python LP01_extract_CT_radiomics_reference_derived.py --set pet_spotcheck --workers 2
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import logging
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import os

import numpy as np
import pandas as pd

# -------------------------------------------------------------------------
# Locked paths
# -------------------------------------------------------------------------
PROJECT = Path(os.environ["PSMA_PROJECT_ROOT"])
V1 = PROJECT / "PSMA_AI_Project"
V2 = PROJECT / "PSMA_AI_Project_V2_PET_Derived"
V2_MODULE_FILE = Path(__file__).resolve().parents[1] / "02_V2_PET_only_candidate_pipeline" / "PSMA_objective9D9E_extract_radiomics_DIRECT_RESUMABLE_LONG_PATH_FIXED.py"
DEEP_PSMA = Path(os.environ["PSMA_EXTERNAL_DATA"])

V1_5C_DIR = V1 / "06_results" / "Objective5C FINAL-FINAL-FINALL routine plus_Radiomics"
V1_5C_PET_TABLE = V1_5C_DIR / "PSMA_objective5C_PET_radiomics_FULL.csv"
V1_5C_MASK_DIR = V1_5C_DIR / "candidate_masks"
V1_8B_TABLE = (V1 / "06_results" / "Objective8_Corrected_External_Validation"
               / "Objective8B_Exact_InternalStyle_Candidate_Dataset"
               / "PSMA_objective8B_exact_internal_style_candidates_20260714_125116.csv")
LOCK_TABLE = (V1 / "06_results" / "Objective7_DEEP_PSMA_External_Validation"
              / "Objective7A4_Locked_Eligibility"
              / "PSMA_objective7_external_validation_eligibility_LOCKED_95.csv")
V2_9B_INTERNAL = (V2 / "06_results" / "Objective9_PET_Derived_Candidate_Pipeline"
                  / "Objective9B_GT_Overlap_and_Label_Assignment" / "internal")
V2_CT_INTERNAL = (V2 / "06_results" / "Objective9_PET_Derived_Candidate_Pipeline"
                  / "Objective9D9E_Direct_PET_CT_Radiomics" / "internal" / "CT"
                  / "PSMA_objective9CT_internal_CT_radiomics_checkpoint.csv")
OUT_ROOT = V2 / "06_results" / "LeakagePaper_FeatureExtraction"

SPOTCHECK_N = 20
SPOTCHECK_SEED = 20261109
SPOTCHECK_REL_TOL = 0.005
EXPECTED = {"internal_ct": (107, 1970), "external_ct": (95, 4658), "pet_spotcheck": (None, 20)}

# -------------------------------------------------------------------------
# Import the V2 extraction module unchanged (also in spawned worker processes)
# -------------------------------------------------------------------------
_MOD_NAME = "psma_v2_9d9e"
_spec = importlib.util.spec_from_file_location(_MOD_NAME, str(V2_MODULE_FILE))
V2MOD = importlib.util.module_from_spec(_spec)
sys.modules[_MOD_NAME] = V2MOD
_spec.loader.exec_module(V2MOD)


# -------------------------------------------------------------------------
# Build candidate records
# -------------------------------------------------------------------------
def v2_internal_pet_paths() -> pd.Series:
    files = list(V2_9B_INTERNAL.glob("*/checkpoints/PSMA_objective9B_candidate_overlap_checkpoint.csv"))
    if len(files) != 3:
        raise RuntimeError(f"Expected 3 internal 9B checkpoints, found {len(files)}")
    t = pd.concat([pd.read_csv(f, usecols=["case_id", "PET_path"]) for f in files])
    return t.drop_duplicates("case_id").set_index("case_id")["PET_path"]


def internal_records(spotcheck: bool) -> pd.DataFrame:
    v1 = pd.read_csv(V1_5C_PET_TABLE, usecols=["candidate_global_id", "case_id",
                     "candidate_label_binary", "label_value_in_mask", "Mask"], low_memory=False)
    v1["candidate_labelmap_path"] = [str(V1_5C_MASK_DIR / Path(str(m).replace("\\", "/")).name) for m in v1.Mask]
    pet = v2_internal_pet_paths()
    v1["PET_path"] = v1.case_id.map(pet)
    if v1.PET_path.isna().any():
        raise RuntimeError("Some V1 internal studies have no V2 PET path")
    if spotcheck:
        pos = v1[v1.candidate_label_binary == 1]
        sel = pos.sample(n=SPOTCHECK_N, random_state=SPOTCHECK_SEED)
    else:
        ct = pd.read_csv(V2_CT_INTERNAL, usecols=["case_id", "CT_resampled_to_PET_grid"])
        resampled = set(ct.loc[ct.CT_resampled_to_PET_grid == 1, "case_id"])
        sel = v1[v1.case_id.isin(resampled)]
    out = pd.DataFrame({
        "candidate_global_id": sel.candidate_global_id.astype(str),
        "case_id": sel.case_id.astype(str),
        "cohort": "internal",
        "candidate_label_value": sel.label_value_in_mask.astype(int),
        "candidate_labelmap_path": sel.candidate_labelmap_path,
        "PET_path": sel.PET_path,
        "v1_candidate_label_binary": sel.candidate_label_binary.astype(int),
        "source_table": "V1_Objective5C",
    })
    return out


def external_records() -> pd.DataFrame:
    lock = pd.read_csv(LOCK_TABLE)
    inc = set(lock.loc[lock.include_external_validation.astype(str).str.lower().isin(["1", "true"]), "case_id"])
    b = pd.read_csv(V1_8B_TABLE, usecols=["candidate_global_id", "case_id", "candidate_label_value",
                    "candidate_label_binary", "candidate_label_map_path"], low_memory=False)
    if set(b.case_id) != inc:
        raise RuntimeError("V1 8B case set differs from the locked 95-case cohort")
    out = pd.DataFrame({
        "candidate_global_id": b.candidate_global_id.astype(str),
        "case_id": b.case_id.astype(str),
        "cohort": "external",
        "candidate_label_value": b.candidate_label_value.astype(int),
        "candidate_labelmap_path": b.candidate_label_map_path.astype(str),
        "PET_path": [str(DEEP_PSMA / c / "PSMA" / "PET.nii.gz") for c in b.case_id],
        "v1_candidate_label_binary": b.candidate_label_binary.astype(int),
        "source_table": "V1_Objective8B",
    })
    return out


def validate(rec: pd.DataFrame, which: str) -> None:
    n_studies, n_rows = EXPECTED[which]
    if len(rec) != n_rows:
        raise RuntimeError(f"{which}: expected {n_rows} candidates, found {len(rec)}")
    if n_studies is not None and rec.case_id.nunique() != n_studies:
        raise RuntimeError(f"{which}: expected {n_studies} studies, found {rec.case_id.nunique()}")
    if rec.candidate_global_id.duplicated().any():
        raise RuntimeError(f"{which}: duplicate candidate IDs")
    if rec.duplicated(["case_id", "candidate_label_value"]).any():
        raise RuntimeError(f"{which}: duplicate label values within a study")
    missing = [p for p in set(rec.candidate_labelmap_path) | set(rec.PET_path) if not Path(p).is_file()]
    if missing:
        raise RuntimeError(f"{which}: {len(missing)} input files missing, e.g. {missing[:3]}")


# -------------------------------------------------------------------------
# Output helpers
# -------------------------------------------------------------------------
def append_csv(path: Path, rows: list) -> None:
    if not rows:
        return
    df = pd.DataFrame(rows)
    if path.exists():
        header = pd.read_csv(path, nrows=0).columns.tolist()
        for c in header:
            if c not in df.columns:
                df[c] = ""
        extra = [c for c in df.columns if c not in header]
        if extra:
            raise RuntimeError(f"New columns not in existing checkpoint: {extra[:5]}")
        df = df[header]
        df.to_csv(path, mode="a", header=False, index=False, quoting=csv.QUOTE_MINIMAL)
    else:
        df.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def versions() -> dict:
    v = {"python": sys.version.split()[0], "platform": platform.platform()}
    for name in ["radiomics", "SimpleITK", "numpy", "pandas"]:
        try:
            v[name] = __import__(name).__version__
        except Exception as exc:
            v[name] = f"unavailable ({exc})"
    v["V2_module"] = str(V2_MODULE_FILE)
    v["PET_yaml"] = str(V2MOD.PET_YAML)
    v["CT_yaml"] = str(V2MOD.CT_YAML)
    v["CROP_PAD_VOX"] = V2MOD.CROP_PAD_VOX
    return v


# -------------------------------------------------------------------------
# PET spot-check comparison
# -------------------------------------------------------------------------
def compare_spotcheck(out_csv: Path, out_dir: Path) -> None:
    new = pd.read_csv(out_csv, low_memory=False)
    old = pd.read_csv(V1_5C_PET_TABLE, low_memory=False)
    feats = [c for c in new.columns if c.startswith("PET_original_") and c in old.columns]
    m = new[["candidate_global_id"] + feats].merge(old[["candidate_global_id"] + feats],
                                                   on="candidate_global_id", suffixes=("_v2", "_v1"))
    rows = []
    for f in feats:
        a = m[f + "_v1"].astype(float).to_numpy()
        b = m[f + "_v2"].astype(float).to_numpy()
        den = np.maximum(np.abs(a), np.abs(b))
        rel = np.divide(np.abs(a - b), den, out=np.zeros_like(den), where=den > 0)
        rows.append({"feature": f, "max_rel_diff": float(rel.max()), "median_rel_diff": float(np.median(rel))})
    r = pd.DataFrame(rows)
    r.to_csv(out_dir / "pet_spotcheck_comparison.csv", index=False)
    worst = r.max_rel_diff.max()
    verdict = "PASS" if (len(m) == SPOTCHECK_N and worst <= SPOTCHECK_REL_TOL) else "FAIL"
    logging.info("PET spot-check: %d candidates, %d features, max relative difference %.3g -> %s "
                 "(tolerance %.3g)", len(m), len(feats), worst, verdict, SPOTCHECK_REL_TOL)


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", required=True, choices=["internal_ct", "external_ct", "pet_spotcheck"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--retry-failures", action="store_true")
    args = ap.parse_args()

    which = args.set
    modality = "PET" if which == "pet_spotcheck" else "CT"
    out_dir = OUT_ROOT / which
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"LP01_{which}_{modality}_radiomics_checkpoint.csv"
    fail_csv = out_dir / f"LP01_{which}_{modality}_radiomics_failures.csv"

    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout),
                                  logging.FileHandler(out_dir / f"LP01_{which}.log", encoding="utf-8")])
    ver = versions()
    logging.info("LP01 %s | modality %s | versions %s", which, modality, json.dumps(ver))

    rec = external_records() if which == "external_ct" else internal_records(spotcheck=(which == "pet_spotcheck"))
    validate(rec, which)
    manifest = out_dir / f"LP01_{which}_input_manifest.csv"
    if not manifest.exists():
        rec.to_csv(manifest, index=False)
    logging.info("Candidates: %d in %d studies (positives %d, A-hist negatives %d)", len(rec),
                 rec.case_id.nunique(), int((rec.v1_candidate_label_binary == 1).sum()),
                 int((rec.v1_candidate_label_binary == 0).sum()))

    done = set()
    if out_csv.exists():
        done = set(pd.read_csv(out_csv, usecols=["candidate_global_id"]).candidate_global_id.astype(str))
    if fail_csv.exists() and not args.retry_failures:
        done |= set(pd.read_csv(fail_csv, usecols=["candidate_global_id"]).candidate_global_id.astype(str))
    todo = rec[~rec.candidate_global_id.isin(done)]
    logging.info("Already completed: %d | to process now: %d", len(rec) - len(todo), len(todo))

    payloads = [(cid, modality, g.to_dict("records"), str(V2MOD.PET_YAML), str(V2MOD.CT_YAML))
                for cid, g in todo.groupby("case_id", sort=True)]
    t0 = time.time(); n_ok = n_err = 0
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = [ex.submit(V2MOD.process_case, p) for p in payloads]
        for i, fu in enumerate(as_completed(futs), 1):
            case_id, ok, err, sec = fu.result()
            append_csv(out_csv, ok)
            append_csv(fail_csv, err)
            n_ok += len(ok); n_err += len(err)
            logging.info("[%d/%d] %s | success=%d failures=%d | %.1fs", i, len(payloads), case_id, len(ok), len(err), sec)

    total_ok = len(pd.read_csv(out_csv, usecols=["candidate_global_id"])) if out_csv.exists() else 0
    state = {"set": which, "modality": modality, "expected": len(rec), "success_total": total_ok,
             "new_success": n_ok, "new_failures": n_err, "minutes": round((time.time() - t0) / 60, 2),
             "versions": ver, "updated": time.strftime("%Y-%m-%d %H:%M:%S")}
    (out_dir / f"LP01_{which}_state.json").write_text(json.dumps(state, indent=2), encoding="utf-8")
    logging.info("DONE: %d / %d successful | new failures this run: %d", total_ok, len(rec), n_err)
    if which == "pet_spotcheck" and total_ok == len(rec):
        compare_spotcheck(out_csv, out_dir)
    return 0 if total_ok == len(rec) else 1


if __name__ == "__main__":
    sys.exit(main())
