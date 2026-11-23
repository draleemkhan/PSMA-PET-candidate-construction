#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PSMA Objective 9C — Final routine-feature extraction
=====================================================

Processes:
    internal/train
    internal/validation
    internal/test
    external/external

Key corrections in this final version:
1. Automatically creates/uses a short Q: drive alias for the Objective 9 root.
2. Supports internal nnU-Net naming:
       *_0001.nii.gz = PET
       *_0000.nii.gz = CT
3. Supports external naming:
       PET.nii.gz
       CT.nii.gz
4. Resolves stale candidate-label-map paths by exact filename/case-ID search.
5. Uses affine-aware SimpleITK CT-to-PET resampling.
6. Resumable: skips candidate IDs already present in a valid checkpoint.
7. Saves checkpoint after every completed case.
8. Handles empty/malformed old checkpoints safely.
9. Clear failure reporting; no misleading secondary KeyError.

Run examples:
    python PSMA_objective9C_extract_routine_features_FINAL.py --cohort validation
    python PSMA_objective9C_extract_routine_features_FINAL.py --cohort all
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import SimpleITK as sitk


# -------------------------------------------------------------------------
# PROJECT PATHS
# -------------------------------------------------------------------------

LONG_OBJ9_ROOT = (Path(os.environ["PSMA_PROJECT_ROOT"]) / "PSMA_AI_Project_V2_PET_Derived" / "06_results" / "Objective9_LockedExternal_Addendum_train0096_0100")

SHORT_DRIVE = "Q:"
SHORT_OBJ9_ROOT = Path("Q:/")

OBJ9A_FOLDER = "Objective9A_PET_Only_Candidate_Generation"
OBJ9C_FOLDER = "Objective9C_Routine_Feature_Extraction"


def ensure_short_drive() -> Path:
    """
    Return a short Objective 9 root.

    If O: is already mapped to the correct directory, use it.
    If O: does not exist, create the mapping automatically.
    If O: exists but points elsewhere, stop safely.
    """
    if not LONG_OBJ9_ROOT.is_dir():
        raise FileNotFoundError(f"Objective 9 root not found:\n{LONG_OBJ9_ROOT}")

    try:
        result = subprocess.run(
            ["subst"],
            capture_output=True,
            text=True,
            check=False,
            shell=True,
        )
        mappings = result.stdout or ""
    except Exception:
        mappings = ""

    normalized_target = str(LONG_OBJ9_ROOT).rstrip("\\").lower()

    for line in mappings.splitlines():
        lower = line.strip().lower()
        if lower.startswith("o:\\:"):
            if normalized_target in lower:
                return SHORT_OBJ9_ROOT
            raise RuntimeError(
                "Drive Q: is already mapped to a different location.\n"
                "Remove it with: subst Q: /D"
            )

    create = subprocess.run(
        ["subst", SHORT_DRIVE, str(LONG_OBJ9_ROOT)],
        capture_output=True,
        text=True,
        check=False,
        shell=True,
    )
    if create.returncode != 0:
        raise RuntimeError(
            "Could not create short Q: drive alias.\n"
            f"stdout: {create.stdout}\n"
            f"stderr: {create.stderr}"
        )

    if not SHORT_OBJ9_ROOT.is_dir():
        raise RuntimeError("Q: drive alias was created but cannot be accessed.")

    return SHORT_OBJ9_ROOT


OBJ9_ROOT = ensure_short_drive()
OBJ9A_ROOT = OBJ9_ROOT / OBJ9A_FOLDER
OUT_ROOT = OBJ9_ROOT / OBJ9C_FOLDER


# -------------------------------------------------------------------------
# TABLE DEFINITIONS
# -------------------------------------------------------------------------

REQUIRED_COLUMNS = [
    "candidate_global_id",
    "case_id",
    "candidate_label_value",
    "candidate_voxels",
    "candidate_volume_mL",
    "PET_path",
    "candidate_labelmap_path",
]

OUTPUT_COLUMNS = [
    "candidate_global_id",
    "dataset_name",
    "cohort",
    "split",
    "case_id",
    "candidate_label_value",
    "candidate_voxels",
    "candidate_volume_mL",
    "SUVmax",
    "SUVmean",
    "SUVmedian",
    "SUVstd",
    "CTmean",
    "CTmedian",
    "CTstd",
    "centroid_x_norm",
    "centroid_y_norm",
    "centroid_z_norm",
    "distance_from_center_mm",
    "bbox_x_mm",
    "bbox_y_mm",
    "bbox_z_mm",
    "lateral_zone",
    "AP_zone",
    "craniocaudal_zone",
    "scanner",
    "radionuclide",
    "contrast",
    "PET_path",
    "CT_path",
    "candidate_labelmap_path",
    "processing_status",
    "error_message",
]

NUMERIC_FEATURES = [
    "candidate_volume_mL",
    "SUVmax",
    "SUVmean",
    "SUVmedian",
    "SUVstd",
    "CTmean",
    "CTmedian",
    "CTstd",
    "centroid_x_norm",
    "centroid_y_norm",
    "centroid_z_norm",
    "distance_from_center_mm",
    "bbox_x_mm",
    "bbox_y_mm",
    "bbox_z_mm",
]


# -------------------------------------------------------------------------
# HELPERS
# -------------------------------------------------------------------------

def configure_logging(log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_path, encoding="utf-8"),
        ],
        force=True,
    )


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(temp, index=False)
    os.replace(temp, path)


def safe_text(row: pd.Series, name: str) -> str:
    if name not in row.index or pd.isna(row[name]):
        return ""
    return str(row[name])


def latest_candidate_csv(scope: str, split: str) -> Path:
    folder = OBJ9A_ROOT / scope / split
    if not folder.is_dir():
        raise FileNotFoundError(f"Objective 9A folder not found:\n{folder}")

    finals = [
        p for p in folder.glob("PSMA_objective9A_PET_only_candidate_index_*.csv")
        if "checkpoint" not in p.name.lower()
    ]
    if finals:
        return max(finals, key=lambda p: p.stat().st_mtime)

    checkpoint = folder / "checkpoints" / "PSMA_objective9A_candidate_index_checkpoint.csv"
    if checkpoint.is_file():
        return checkpoint

    raise FileNotFoundError(f"No Objective 9A candidate table found in:\n{folder}")


def read_candidate_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype={"candidate_global_id": str, "case_id": str})
    df.columns = [str(c).strip() for c in df.columns]

    if (
        "candidate_labelmap_path" not in df.columns
        and "candidate_label_map_path" in df.columns
    ):
        df = df.rename(
            columns={"candidate_label_map_path": "candidate_labelmap_path"}
        )

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise RuntimeError(f"Candidate table missing columns: {missing}")

    if df["candidate_global_id"].isna().any():
        raise RuntimeError("Candidate table contains blank candidate_global_id values.")

    duplicate_count = int(df["candidate_global_id"].duplicated().sum())
    if duplicate_count:
        raise RuntimeError(
            f"Candidate table contains {duplicate_count} duplicate candidate IDs."
        )

    return df


def convert_to_short_path(path: Path) -> Path:
    """
    Convert any path located under the long Objective 9 root to O:/...
    """
    text = str(path)
    long_text = str(LONG_OBJ9_ROOT)
    if text.lower().startswith(long_text.lower()):
        relative = text[len(long_text):].lstrip("\\/")
        return SHORT_OBJ9_ROOT / Path(relative)
    return path


def resolve_labelmap_path(
    csv_path: Path,
    case_id: str,
    scope: str,
    split: str,
) -> Path:
    """
    Resolve candidate label map robustly.

    Priority:
    1. CSV path converted to short O: path.
    2. Expected Objective 9A candidate_label_maps folder.
    3. Exact filename recursive search.
    4. Case-ID + labelmap recursive search.
    """
    short_csv_path = convert_to_short_path(csv_path)
    if short_csv_path.is_file():
        return short_csv_path

    root = OBJ9A_ROOT / scope / split
    expected_folder = root / "candidate_label_maps"

    exact_expected = expected_folder / csv_path.name
    if exact_expected.is_file():
        return exact_expected

    exact_matches = list(root.rglob(csv_path.name))
    exact_matches = [p for p in exact_matches if p.is_file()]
    if len(exact_matches) == 1:
        return exact_matches[0]
    if len(exact_matches) > 1:
        raise RuntimeError(
            f"{case_id}: multiple exact label-map matches found:\n"
            + "\n".join(str(p) for p in exact_matches)
        )

    case_matches: List[Path] = []
    for p in root.rglob("*.nii*"):
        name = p.name.lower()
        if (
            case_id.lower() in name
            and "candidate" in name
            and "labelmap" in name
        ):
            case_matches.append(p)

    if len(case_matches) == 1:
        return case_matches[0]
    if len(case_matches) > 1:
        raise RuntimeError(
            f"{case_id}: multiple candidate label-map matches found:\n"
            + "\n".join(str(p) for p in case_matches)
        )

    raise FileNotFoundError(
        f"{case_id}: candidate label map not found.\n"
        f"CSV path: {csv_path}\n"
        f"Searched root: {root}"
    )


def find_ct_path(pet_path: Path) -> Path:
    """
    Internal:
        *_0001.nii.gz PET -> *_0000.nii.gz CT

    External:
        PET.nii.gz -> CT.nii.gz
    """
    pet_path = convert_to_short_path(pet_path)

    name = pet_path.name
    lower_name = name.lower()

    if lower_name.endswith("_0001.nii.gz"):
        ct = pet_path.with_name(name[:-12] + "_0000.nii.gz")
        if ct.is_file():
            return ct

    if lower_name.endswith("_0001.nii"):
        ct = pet_path.with_name(name[:-9] + "_0000.nii")
        if ct.is_file():
            return ct

    external_options = [
        pet_path.parent / "CT.nii.gz",
        pet_path.parent / "ct.nii.gz",
        pet_path.parent / "CT.nii",
        pet_path.parent / "ct.nii",
    ]
    for p in external_options:
        if p.is_file():
            return p

    candidates: List[Path] = []
    for pattern in ("*.nii.gz", "*.nii"):
        for p in pet_path.parent.glob(pattern):
            n = p.name.lower()
            if (
                "ct" in n
                and "pet" not in n
                and "fdg" not in n
                and "ttb" not in n
                and "mask" not in n
                and "label" not in n
                and "seg" not in n
            ):
                candidates.append(p)

    candidates = list(dict.fromkeys(candidates))
    if len(candidates) == 1:
        return candidates[0]

    raise FileNotFoundError(
        f"CT not uniquely found beside PET:\n{pet_path}\n"
        f"Candidates: {candidates}"
    )


def same_geometry(a: sitk.Image, b: sitk.Image, tol: float = 1e-5) -> bool:
    return (
        a.GetSize() == b.GetSize()
        and np.allclose(a.GetSpacing(), b.GetSpacing(), atol=tol, rtol=0)
        and np.allclose(a.GetOrigin(), b.GetOrigin(), atol=tol, rtol=0)
        and np.allclose(a.GetDirection(), b.GetDirection(), atol=tol, rtol=0)
    )


def resample_ct_to_pet(ct: sitk.Image, pet: sitk.Image) -> sitk.Image:
    if same_geometry(ct, pet):
        return sitk.Cast(ct, sitk.sitkFloat32)

    return sitk.Resample(
        ct,
        pet,
        sitk.Transform(),
        sitk.sitkLinear,
        -1024.0,
        sitk.sitkFloat32,
    )


def lateral_zone(x: float) -> str:
    if x < 0.45:
        return "right"
    if x > 0.55:
        return "left"
    return "midline"


def ap_zone(y: float) -> str:
    if y < 0.45:
        return "posterior"
    if y > 0.55:
        return "anterior"
    return "central"


def cc_zone(z: float) -> str:
    if z < 1.0 / 3.0:
        return "superior"
    if z > 2.0 / 3.0:
        return "inferior"
    return "middle"


# -------------------------------------------------------------------------
# CASE EXTRACTION
# -------------------------------------------------------------------------

def extract_case(
    case_df: pd.DataFrame,
    scope: str,
    split: str,
) -> List[Dict[str, object]]:
    case_id = str(case_df["case_id"].iloc[0])

    pet_values = case_df["PET_path"].astype(str).unique()
    map_values = case_df["candidate_labelmap_path"].astype(str).unique()

    if len(pet_values) != 1:
        raise RuntimeError(f"{case_id}: multiple PET paths in candidate table.")
    if len(map_values) != 1:
        raise RuntimeError(f"{case_id}: multiple label-map paths in candidate table.")

    pet_path = convert_to_short_path(Path(pet_values[0]))
    label_path = resolve_labelmap_path(
        Path(map_values[0]), case_id, scope, split
    )
    ct_path = find_ct_path(pet_path)

    for p, label in [
        (pet_path, "PET"),
        (ct_path, "CT"),
        (label_path, "candidate label map"),
    ]:
        if not p.is_file():
            raise FileNotFoundError(f"{case_id}: {label} file not found:\n{p}")

    pet_img = sitk.ReadImage(str(pet_path))
    ct_img_native = sitk.ReadImage(str(ct_path))
    label_img = sitk.ReadImage(str(label_path))

    if not same_geometry(pet_img, label_img):
        raise RuntimeError(
            f"{case_id}: PET and candidate label-map geometry differ."
        )

    ct_img = resample_ct_to_pet(ct_img_native, pet_img)

    pet = sitk.GetArrayFromImage(pet_img).astype(np.float64, copy=False)
    ct = sitk.GetArrayFromImage(ct_img).astype(np.float64, copy=False)
    labels = sitk.GetArrayFromImage(label_img)

    expected_labels = set(
        case_df["candidate_label_value"].astype(int).tolist()
    )
    present_labels = set(np.unique(labels).astype(int).tolist())
    present_labels.discard(0)

    if expected_labels != present_labels:
        missing = sorted(expected_labels - present_labels)
        extra = sorted(present_labels - expected_labels)
        raise RuntimeError(
            f"{case_id}: candidate-label mismatch. "
            f"Missing={missing[:10]}, extra={extra[:10]}"
        )

    spacing_xyz = np.asarray(pet_img.GetSpacing(), dtype=float)
    spacing_zyx = spacing_xyz[::-1]
    size_zyx = np.asarray(labels.shape, dtype=float)
    image_center_zyx = (size_zyx - 1.0) / 2.0
    voxel_volume_ml = float(np.prod(spacing_xyz) / 1000.0)

    output: List[Dict[str, object]] = []

    for _, row in case_df.iterrows():
        label_value = int(row["candidate_label_value"])
        coords = np.argwhere(labels == label_value)

        if coords.size == 0:
            raise RuntimeError(
                f"{case_id}: label {label_value} has no voxels."
            )

        zyx = tuple(coords.T)
        pet_values_candidate = pet[zyx]
        ct_values_candidate = ct[zyx]

        pet_values_candidate = pet_values_candidate[
            np.isfinite(pet_values_candidate)
        ]
        ct_values_candidate = ct_values_candidate[
            np.isfinite(ct_values_candidate)
        ]

        if pet_values_candidate.size == 0:
            raise RuntimeError(
                f"{case_id}: label {label_value} has no finite PET values."
            )
        if ct_values_candidate.size == 0:
            raise RuntimeError(
                f"{case_id}: label {label_value} has no finite CT values."
            )

        candidate_voxels = int(coords.shape[0])
        volume_ml = candidate_voxels * voxel_volume_ml

        centroid_zyx = coords.mean(axis=0)
        centroid_norm_zyx = (centroid_zyx + 1.0) / size_zyx

        bbox_voxels_zyx = coords.max(axis=0) - coords.min(axis=0) + 1
        bbox_mm_zyx = bbox_voxels_zyx * spacing_zyx

        delta_mm_zyx = (
            centroid_zyx - image_center_zyx
        ) * spacing_zyx
        distance_mm = float(np.sqrt(np.sum(delta_mm_zyx ** 2)))

        cz, cy, cx = centroid_norm_zyx.tolist()
        bbox_z, bbox_y, bbox_x = bbox_mm_zyx.tolist()

        result = {
            "candidate_global_id": str(row["candidate_global_id"]),
            "dataset_name": safe_text(row, "dataset_name"),
            "cohort": safe_text(row, "cohort"),
            "split": safe_text(row, "split"),
            "case_id": case_id,
            "candidate_label_value": label_value,
            "candidate_voxels": candidate_voxels,
            "candidate_volume_mL": float(volume_ml),
            "SUVmax": float(np.max(pet_values_candidate)),
            "SUVmean": float(np.mean(pet_values_candidate)),
            "SUVmedian": float(np.median(pet_values_candidate)),
            "SUVstd": (
                float(np.std(pet_values_candidate, ddof=1))
                if pet_values_candidate.size > 1
                else 0.0
            ),
            "CTmean": float(np.mean(ct_values_candidate)),
            "CTmedian": float(np.median(ct_values_candidate)),
            "CTstd": (
                float(np.std(ct_values_candidate, ddof=1))
                if ct_values_candidate.size > 1
                else 0.0
            ),
            "centroid_x_norm": float(cx),
            "centroid_y_norm": float(cy),
            "centroid_z_norm": float(cz),
            "distance_from_center_mm": distance_mm,
            "bbox_x_mm": float(bbox_x),
            "bbox_y_mm": float(bbox_y),
            "bbox_z_mm": float(bbox_z),
            "lateral_zone": lateral_zone(float(cx)),
            "AP_zone": ap_zone(float(cy)),
            "craniocaudal_zone": cc_zone(float(cz)),
            "scanner": safe_text(row, "scanner"),
            "radionuclide": safe_text(row, "radionuclide"),
            "contrast": safe_text(row, "contrast"),
            "PET_path": str(pet_path),
            "CT_path": str(ct_path),
            "candidate_labelmap_path": str(label_path),
            "processing_status": "success",
            "error_message": "",
        }

        # Safety reconciliation against Objective 9A values.
        original_voxels = int(row["candidate_voxels"])
        if original_voxels != candidate_voxels:
            raise RuntimeError(
                f"{case_id}: voxel-count mismatch for "
                f"{result['candidate_global_id']}: "
                f"Objective9A={original_voxels}, recalculated={candidate_voxels}"
            )

        original_volume = float(row["candidate_volume_mL"])
        if abs(original_volume - volume_ml) > 1e-5:
            raise RuntimeError(
                f"{case_id}: volume mismatch for "
                f"{result['candidate_global_id']}: "
                f"Objective9A={original_volume}, recalculated={volume_ml}"
            )

        output.append(result)

    return output


# -------------------------------------------------------------------------
# COHORT RUNNER
# -------------------------------------------------------------------------

def load_valid_checkpoint(
    checkpoint_path: Path,
    expected_ids: set[str],
) -> Tuple[List[Dict[str, object]], set[str]]:
    if not checkpoint_path.is_file():
        return [], set()

    try:
        checkpoint = pd.read_csv(
            checkpoint_path,
            dtype={"candidate_global_id": str},
        )
    except Exception:
        return [], set()

    if checkpoint.empty:
        return [], set()

    if "candidate_global_id" not in checkpoint.columns:
        return [], set()

    if not set(OUTPUT_COLUMNS).issubset(checkpoint.columns):
        return [], set()

    checkpoint = checkpoint[OUTPUT_COLUMNS].copy()
    checkpoint = checkpoint.drop_duplicates(
        subset=["candidate_global_id"],
        keep="last",
    )

    checkpoint_ids = set(
        checkpoint["candidate_global_id"].astype(str)
    )

    unknown = checkpoint_ids - expected_ids
    if unknown:
        raise RuntimeError(
            f"Checkpoint contains {len(unknown)} unknown candidate IDs."
        )

    return checkpoint.to_dict("records"), checkpoint_ids


def run_one(scope: str, split: str) -> None:
    candidate_csv = latest_candidate_csv(scope, split)

    out_dir = OUT_ROOT / scope / split
    out_dir.mkdir(parents=True, exist_ok=True)

    checkpoint_path = (
        out_dir
        / f"PSMA_objective9C_{scope}_{split}_routine_features_checkpoint.csv"
    )
    failures_path = (
        out_dir
        / f"PSMA_objective9C_{scope}_{split}_routine_features_failures.csv"
    )
    summary_path = (
        out_dir
        / f"PSMA_objective9C_{scope}_{split}_summary.csv"
    )
    state_path = (
        out_dir
        / f"PSMA_objective9C_{scope}_{split}_state.json"
    )
    log_path = (
        out_dir
        / f"PSMA_objective9C_{scope}_{split}.log"
    )

    configure_logging(log_path)

    logging.info("=" * 80)
    logging.info("PSMA OBJECTIVE 9C — FINAL ROUTINE FEATURE EXTRACTION")
    logging.info("Scope/split: %s/%s", scope, split)
    logging.info("Objective 9 root: %s", OBJ9_ROOT)
    logging.info("Candidate table: %s", candidate_csv)
    logging.info("=" * 80)

    candidates = read_candidate_table(candidate_csv)
    expected_ids = set(
        candidates["candidate_global_id"].astype(str)
    )

    records, completed_ids = load_valid_checkpoint(
        checkpoint_path,
        expected_ids,
    )

    failures: List[Dict[str, object]] = []
    if failures_path.is_file():
        try:
            old_failures = pd.read_csv(failures_path)
            failures = old_failures.to_dict("records")
        except Exception:
            failures = []

    grouped = list(candidates.groupby("case_id", sort=False))
    total_cases = len(grouped)

    for case_number, (case_id, case_df) in enumerate(grouped, start=1):
        case_ids = set(
            case_df["candidate_global_id"].astype(str)
        )

        if case_ids.issubset(completed_ids):
            logging.info(
                "Case %d/%d already complete: %s",
                case_number,
                total_cases,
                case_id,
            )
            continue

        logging.info(
            "Case %d/%d: %s | candidates=%d",
            case_number,
            total_cases,
            case_id,
            len(case_df),
        )

        started = time.time()

        try:
            case_records = extract_case(
                case_df=case_df,
                scope=scope,
                split=split,
            )

            new_records = [
                row for row in case_records
                if row["candidate_global_id"] not in completed_ids
            ]

            records.extend(new_records)
            completed_ids.update(
                row["candidate_global_id"] for row in new_records
            )

            out = pd.DataFrame(records)
            out = out[OUTPUT_COLUMNS]

            original_order = {
                candidate_id: index
                for index, candidate_id in enumerate(
                    candidates["candidate_global_id"].astype(str)
                )
            }
            out["_sort_order"] = (
                out["candidate_global_id"]
                .map(original_order)
            )
            out = (
                out.sort_values("_sort_order")
                .drop(columns="_sort_order")
            )

            atomic_csv(out, checkpoint_path)

            state = {
                "scope": scope,
                "split": split,
                "candidate_input": str(candidate_csv),
                "expected_candidates": int(len(candidates)),
                "successful_candidates": int(len(completed_ids)),
                "completed_cases_seen": int(case_number),
                "total_cases": int(total_cases),
                "updated_at": pd.Timestamp.now().isoformat(),
            }
            state_path.write_text(
                json.dumps(state, indent=2),
                encoding="utf-8",
            )

            logging.info(
                "Progress saved: %d/%d candidates | %.1f sec",
                len(completed_ids),
                len(candidates),
                time.time() - started,
            )

        except Exception as exc:
            logging.error(
                "Case failed: %s | %s",
                case_id,
                exc,
            )
            logging.debug(traceback.format_exc())

            failures.append(
                {
                    "scope": scope,
                    "split": split,
                    "case_id": str(case_id),
                    "candidate_count": int(len(case_df)),
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "timestamp": pd.Timestamp.now().isoformat(),
                }
            )
            atomic_csv(pd.DataFrame(failures), failures_path)

    if checkpoint_path.is_file():
        final = pd.read_csv(
            checkpoint_path,
            dtype={"candidate_global_id": str},
        )
    else:
        final = pd.DataFrame(columns=OUTPUT_COLUMNS)

    success_ids = (
        set(final["candidate_global_id"].astype(str))
        if "candidate_global_id" in final.columns
        else set()
    )
    missing_ids = expected_ids - success_ids

    nonfinite_rows = 0
    if not final.empty:
        nonfinite_mask = np.zeros(len(final), dtype=bool)
        for column in NUMERIC_FEATURES:
            values = pd.to_numeric(
                final[column],
                errors="coerce",
            ).to_numpy(dtype=float)
            nonfinite_mask |= ~np.isfinite(values)
        nonfinite_rows = int(nonfinite_mask.sum())

    summary = pd.DataFrame(
        {
            "metric": [
                "Expected candidate rows",
                "Successful candidate rows",
                "Missing candidate rows",
                "Unique candidate IDs",
                "Rows with non-finite numeric features",
                "Cases",
                "Objective 9 root used",
                "Candidate input",
            ],
            "value": [
                len(candidates),
                len(final),
                len(missing_ids),
                (
                    final["candidate_global_id"].nunique()
                    if "candidate_global_id" in final.columns
                    else 0
                ),
                nonfinite_rows,
                candidates["case_id"].nunique(),
                str(OBJ9_ROOT),
                str(candidate_csv),
            ],
        }
    )
    atomic_csv(summary, summary_path)

    logging.info("=" * 80)
    logging.info(
        "Success total: %d/%d",
        len(final),
        len(candidates),
    )
    logging.info("Missing: %d", len(missing_ids))
    logging.info("Non-finite rows: %d", nonfinite_rows)
    logging.info("Checkpoint: %s", checkpoint_path)
    logging.info("Failures: %s", failures_path)
    logging.info("Summary: %s", summary_path)
    logging.info("=" * 80)

    if missing_ids or nonfinite_rows:
        raise RuntimeError(
            f"Safety stop {scope}/{split}: "
            f"missing={len(missing_ids)}, "
            f"nonfinite={nonfinite_rows}"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cohort",
        choices=[
            "train",
            "validation",
            "test",
            "external",
            "all",
        ],
        default="all",
    )
    args = parser.parse_args()

    jobs = {
        "train": ("internal", "train"),
        "validation": ("internal", "validation"),
        "test": ("internal", "test"),
        "external": ("external", "external"),
    }

    selected = (
        list(jobs.keys())
        if args.cohort == "all"
        else [args.cohort]
    )

    for name in selected:
        scope, split = jobs[name]
        run_one(scope, split)


if __name__ == "__main__":
    main()
