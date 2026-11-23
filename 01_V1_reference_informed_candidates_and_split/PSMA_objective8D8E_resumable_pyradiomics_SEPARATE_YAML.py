#!/usr/bin/env python
"""
PSMA Objectives 8D + 8E
Resumable PET and CT PyRadiomics extraction for the corrected external cohort.

Designed for Windows PowerShell.

Key behavior
------------
- Reads the latest Objective 8C merged candidate CSV unless paths are supplied.
- Uses Objective 8B candidate label maps.
- Extracts PET and CT radiomics candidate-by-candidate.
- Resumes from existing checkpoint CSVs.
- Saves progress every N candidates and every M minutes.
- Writes failures separately and continues.
- Uses one exact PyRadiomics parameter YAML for reproducibility.
- CT is resampled to PET/candidate-label-map geometry using SimpleITK.
- PET is used in its native PET geometry.
- Candidate labels are read from candidate_label_value.
- Diagnostic PyRadiomics fields are excluded from final feature tables.
- Existing successful candidate IDs are skipped on resume.

Requirements
------------
pip install pyradiomics SimpleITK pandas numpy openpyxl

Important
---------
Use the SAME YAML parameter file that was used for the internal radiomics
pipeline. Do not invent a new parameterization for definitive inference.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor


DEFAULT_PROJECT_ROOT = (Path(os.environ["PSMA_PROJECT_ROOT"]) / "PSMA_AI_Project")
DEFAULT_INPUT_ROOT = DEFAULT_PROJECT_ROOT / "06_results" / \
    "Objective8_Corrected_External_Validation" / \
    "Objective8C_InternalStyle_Routine_Features"
DEFAULT_OUTPUT_ROOT = DEFAULT_PROJECT_ROOT / "06_results" / \
    "Objective8_Corrected_External_Validation" / \
    "Objective8D8E_Resumable_PET_CT_Radiomics"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Resumable Objective 8D PET + 8E CT PyRadiomics extraction."
    )
    p.add_argument(
        "--input-csv",
        type=Path,
        default=None,
        help="Objective 8C merged candidate CSV. If omitted, latest matching CSV is used.",
    )
    p.add_argument(
        "--pet-params-yaml",
        type=Path,
        default=None,
        help="Exact internal PET PyRadiomics YAML parameter file.",
    )
    p.add_argument(
        "--ct-params-yaml",
        type=Path,
        default=None,
        help="Exact internal CT PyRadiomics YAML parameter file.",
    )
    p.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Output directory.",
    )
    p.add_argument(
        "--mode",
        choices=["pet", "ct", "both"],
        default="both",
        help="Run PET only, CT only, or both.",
    )
    p.add_argument(
        "--save-every",
        type=int,
        default=25,
        help="Checkpoint after this many newly completed candidates.",
    )
    p.add_argument(
        "--save-minutes",
        type=float,
        default=5.0,
        help="Checkpoint at least this often in minutes.",
    )
    p.add_argument(
        "--retry-failures",
        action="store_true",
        help="Retry candidates already present in the failure CSV.",
    )
    p.add_argument(
        "--max-candidates",
        type=int,
        default=None,
        help="Optional test limit for newly processed candidates.",
    )
    p.add_argument(
        "--case-id",
        type=str,
        default=None,
        help="Optional single case ID for testing.",
    )
    p.add_argument(
        "--verbosity",
        type=int,
        choices=[10, 20, 30, 40],
        default=20,
        help="Logging level: 10 DEBUG, 20 INFO, 30 WARNING, 40 ERROR.",
    )
    return p.parse_args()


def latest_file(folder: Path, pattern: str) -> Path:
    matches = sorted(folder.glob(pattern), key=lambda x: x.stat().st_mtime)
    if not matches:
        raise FileNotFoundError(f"No files matching {pattern!r} under {folder}")
    return matches[-1]


def atomic_write_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_write_json(data: Dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
    os.replace(tmp, path)


def read_existing(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, low_memory=False)


def normalize_path(value: Any) -> Path:
    text = str(value).strip().strip('"')
    if not text:
        raise ValueError("Blank path.")
    return Path(text)


def get_ct_path(row: pd.Series) -> Path:
    for name in ("CT_path", "ct_path"):
        if name in row.index and pd.notna(row[name]) and str(row[name]).strip():
            return normalize_path(row[name])
    raise KeyError("CT_path is missing from Objective 8C table.")


def get_pet_path(row: pd.Series) -> Path:
    for name in ("PET_path", "pet_path"):
        if name in row.index and pd.notna(row[name]) and str(row[name]).strip():
            return normalize_path(row[name])
    raise KeyError("PET_path is missing from Objective 8C table.")


def get_map_path(row: pd.Series) -> Path:
    for name in ("candidate_label_map_path", "label_map_path"):
        if name in row.index and pd.notna(row[name]) and str(row[name]).strip():
            return normalize_path(row[name])
    raise KeyError("candidate_label_map_path is missing.")


def same_geometry(a: sitk.Image, b: sitk.Image, tol: float = 1e-5) -> bool:
    return (
        a.GetSize() == b.GetSize()
        and np.allclose(a.GetSpacing(), b.GetSpacing(), atol=tol, rtol=0)
        and np.allclose(a.GetOrigin(), b.GetOrigin(), atol=tol, rtol=0)
        and np.allclose(a.GetDirection(), b.GetDirection(), atol=tol, rtol=0)
    )


def resample_to_reference(
    image: sitk.Image,
    reference: sitk.Image,
    interpolator: int = sitk.sitkLinear,
    default_value: float = 0.0,
) -> sitk.Image:
    rf = sitk.ResampleImageFilter()
    rf.SetReferenceImage(reference)
    rf.SetTransform(sitk.Transform())
    rf.SetInterpolator(interpolator)
    rf.SetDefaultPixelValue(float(default_value))
    rf.SetOutputPixelType(sitk.sitkFloat32)
    return rf.Execute(image)


def make_binary_mask(label_map: sitk.Image, label_value: int) -> sitk.Image:
    arr = sitk.GetArrayFromImage(label_map)
    mask_arr = (arr == int(label_value)).astype(np.uint8)
    if int(mask_arr.sum()) == 0:
        raise ValueError(f"Label {label_value} is absent from candidate label map.")
    mask = sitk.GetImageFromArray(mask_arr)
    mask.CopyInformation(label_map)
    return mask


def clean_feature_value(v: Any) -> Any:
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, (list, tuple, np.ndarray)):
        a = np.asarray(v)
        if a.size == 1:
            return a.ravel()[0].item()
        return json.dumps(a.tolist())
    return v


def extract_features(
    extractor: featureextractor.RadiomicsFeatureExtractor,
    image: sitk.Image,
    mask: sitk.Image,
) -> Dict[str, Any]:
    result = extractor.execute(image, mask, label=1)
    cleaned: Dict[str, Any] = {}
    for k, v in result.items():
        if str(k).startswith("diagnostics_"):
            continue
        cleaned[str(k)] = clean_feature_value(v)
    if not cleaned:
        raise RuntimeError("PyRadiomics returned no non-diagnostic features.")
    return cleaned


class CheckpointStore:
    def __init__(self, output_root: Path, modality: str) -> None:
        self.modality = modality.upper()
        self.success_path = output_root / f"PSMA_objective8{modality.upper()}_{modality}_radiomics_checkpoint.csv"
        self.failure_path = output_root / f"PSMA_objective8{modality.upper()}_{modality}_radiomics_failures.csv"
        self.state_path = output_root / f"PSMA_objective8{modality.upper()}_{modality}_state.json"
        self.success = read_existing(self.success_path)
        self.failures = read_existing(self.failure_path)

    @property
    def completed_ids(self) -> set[int]:
        if self.success.empty or "candidate_global_id" not in self.success.columns:
            return set()
        return set(pd.to_numeric(self.success["candidate_global_id"], errors="coerce").dropna().astype(int))

    @property
    def failed_ids(self) -> set[int]:
        if self.failures.empty or "candidate_global_id" not in self.failures.columns:
            return set()
        return set(pd.to_numeric(self.failures["candidate_global_id"], errors="coerce").dropna().astype(int))

    def append_success(self, row: Dict[str, Any]) -> None:
        new = pd.DataFrame([row])
        self.success = pd.concat([self.success, new], ignore_index=True, sort=False)
        self.success = self.success.drop_duplicates("candidate_global_id", keep="last")

    def append_failure(self, row: Dict[str, Any]) -> None:
        new = pd.DataFrame([row])
        self.failures = pd.concat([self.failures, new], ignore_index=True, sort=False)
        self.failures = self.failures.drop_duplicates("candidate_global_id", keep="last")

    def remove_failure(self, candidate_id: int) -> None:
        if not self.failures.empty and "candidate_global_id" in self.failures.columns:
            ids = pd.to_numeric(self.failures["candidate_global_id"], errors="coerce")
            self.failures = self.failures.loc[ids != candidate_id].copy()

    def save(self, total_rows: int, input_csv: Path, params_yaml: Path) -> None:
        atomic_write_csv(self.success, self.success_path)
        atomic_write_csv(self.failures, self.failure_path)
        atomic_write_json(
            {
                "modality": self.modality,
                "updated": datetime.now().isoformat(timespec="seconds"),
                "input_csv": str(input_csv),
                "params_yaml": str(params_yaml),
                "total_candidate_rows": int(total_rows),
                "successful_candidates": int(len(self.completed_ids)),
                "failed_candidates": int(len(self.failed_ids)),
            },
            self.state_path,
        )


def process_modality(
    df: pd.DataFrame,
    modality: str,
    args: argparse.Namespace,
    extractor: featureextractor.RadiomicsFeatureExtractor,
    params_yaml: Path,
) -> None:
    modality_lower = modality.lower()
    store = CheckpointStore(args.output_root, modality_lower)
    completed = store.completed_ids
    failed = store.failed_ids

    pending = df.loc[~df["candidate_global_id"].astype(int).isin(completed)].copy()
    if not args.retry_failures:
        pending = pending.loc[~pending["candidate_global_id"].astype(int).isin(failed)].copy()

    logging.info(
        "%s: total=%d completed=%d previous_failures=%d pending=%d",
        modality,
        len(df),
        len(completed),
        len(failed),
        len(pending),
    )

    current_case: Optional[str] = None
    pet_img: Optional[sitk.Image] = None
    ct_img: Optional[sitk.Image] = None
    map_img: Optional[sitk.Image] = None
    case_ct_resampled = False

    newly_processed = 0
    since_save = 0
    last_save = time.monotonic()
    started = time.monotonic()

    for _, row in pending.iterrows():
        if args.max_candidates is not None and newly_processed >= args.max_candidates:
            logging.info("Reached --max-candidates=%d", args.max_candidates)
            break

        cid = int(row["candidate_global_id"])
        case_id = str(row["case_id"])
        label_value = int(row["candidate_label_value"])

        try:
            if case_id != current_case:
                pet_path = get_pet_path(row)
                ct_path = get_ct_path(row)
                map_path = get_map_path(row)

                for p in (pet_path, ct_path, map_path):
                    if not p.is_file():
                        raise FileNotFoundError(str(p))

                pet_img = sitk.ReadImage(str(pet_path), sitk.sitkFloat32)
                ct_img = sitk.ReadImage(str(ct_path), sitk.sitkFloat32)
                map_img = sitk.ReadImage(str(map_path))

                if not same_geometry(pet_img, map_img):
                    raise ValueError(
                        f"PET and label map geometry differ for case {case_id}."
                    )

                case_ct_resampled = not same_geometry(ct_img, map_img)
                if case_ct_resampled:
                    # Use physical-coordinate resampling, not matrix-only resizing.
                    ct_img = resample_to_reference(
                        ct_img,
                        map_img,
                        interpolator=sitk.sitkLinear,
                        default_value=-1024.0,
                    )

                current_case = case_id
                logging.info(
                    "%s case loaded: %s | CT resampled=%s",
                    modality,
                    case_id,
                    case_ct_resampled,
                )

            assert pet_img is not None
            assert ct_img is not None
            assert map_img is not None

            mask = make_binary_mask(map_img, label_value)
            image = pet_img if modality == "PET" else ct_img
            feats = extract_features(extractor, image, mask)

            out: Dict[str, Any] = {
                "candidate_global_id": cid,
                "case_id": case_id,
                "candidate_label_value": label_value,
                "candidate_label_binary": int(row["candidate_label_binary"]),
                "candidate_volume_mL": float(row["candidate_volume_mL"]),
                "modality": modality,
                "CT_resampled_to_PET": bool(case_ct_resampled) if modality == "CT" else False,
            }
            out.update(feats)
            store.append_success(out)
            store.remove_failure(cid)

        except Exception as exc:
            logging.error(
                "%s failed: case=%s candidate=%s label=%s | %s",
                modality,
                case_id,
                cid,
                label_value,
                exc,
            )
            store.append_failure(
                {
                    "candidate_global_id": cid,
                    "case_id": case_id,
                    "candidate_label_value": label_value,
                    "modality": modality,
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "traceback": traceback.format_exc(limit=5),
                    "last_attempt": datetime.now().isoformat(timespec="seconds"),
                }
            )

        newly_processed += 1
        since_save += 1

        elapsed_since_save = time.monotonic() - last_save
        if since_save >= args.save_every or elapsed_since_save >= args.save_minutes * 60:
            store.save(len(df), args.input_csv, params_yaml)
            done = len(store.completed_ids)
            rate = newly_processed / max((time.monotonic() - started) / 60, 1e-9)
            logging.info(
                "%s checkpoint: success=%d/%d failures=%d new_rate=%.2f candidates/min",
                modality,
                done,
                len(df),
                len(store.failed_ids),
                rate,
            )
            since_save = 0
            last_save = time.monotonic()

    store.save(len(df), args.input_csv, params_yaml)

    logging.info(
        "%s complete/current state: success=%d/%d failures=%d",
        modality,
        len(store.completed_ids),
        len(df),
        len(store.failed_ids),
    )


def main() -> int:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)

    log_path = args.output_root / "PSMA_objective8D8E_radiomics.log"
    logging.basicConfig(
        level=args.verbosity,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_path, encoding="utf-8"),
        ],
    )

    if args.input_csv is None:
        args.input_csv = latest_file(
            DEFAULT_INPUT_ROOT,
            "PSMA_objective8C_candidate_table_plus_routine18_*.csv",
        )

    args.input_csv = args.input_csv.resolve()
    args.output_root = args.output_root.resolve()

    if args.pet_params_yaml is not None:
        args.pet_params_yaml = args.pet_params_yaml.resolve()
    if args.ct_params_yaml is not None:
        args.ct_params_yaml = args.ct_params_yaml.resolve()

    if not args.input_csv.is_file():
        raise FileNotFoundError(args.input_csv)

    if args.mode in ("pet", "both"):
        if args.pet_params_yaml is None or not args.pet_params_yaml.is_file():
            raise FileNotFoundError(
                f"Valid --pet-params-yaml is required for mode {args.mode!r}."
            )

    if args.mode in ("ct", "both"):
        if args.ct_params_yaml is None or not args.ct_params_yaml.is_file():
            raise FileNotFoundError(
                f"Valid --ct-params-yaml is required for mode {args.mode!r}."
            )

    logging.info("Input CSV: %s", args.input_csv)
    logging.info("PET parameters YAML: %s", args.pet_params_yaml)
    logging.info("CT parameters YAML: %s", args.ct_params_yaml)
    logging.info("Output root: %s", args.output_root)

    df = pd.read_csv(args.input_csv, low_memory=False)

    required = {
        "candidate_global_id",
        "case_id",
        "candidate_label_value",
        "candidate_label_binary",
        "candidate_volume_mL",
        "PET_path",
        "CT_path",
        "candidate_label_map_path",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    df["candidate_global_id"] = pd.to_numeric(
        df["candidate_global_id"], errors="raise"
    ).astype(int)
    df["candidate_label_value"] = pd.to_numeric(
        df["candidate_label_value"], errors="raise"
    ).astype(int)

    if args.case_id:
        df = df.loc[df["case_id"].astype(str) == args.case_id].copy()
        if df.empty:
            raise ValueError(f"No rows found for --case-id {args.case_id!r}")

    if df["candidate_global_id"].duplicated().any():
        dup = df.loc[df["candidate_global_id"].duplicated(), "candidate_global_id"]
        raise ValueError(f"Duplicate candidate_global_id values: {dup.head().tolist()}")

    if args.mode in ("pet", "both"):
        pet_extractor = featureextractor.RadiomicsFeatureExtractor(
            str(args.pet_params_yaml)
        )
        logging.info("PET PyRadiomics settings loaded successfully.")
        process_modality(
            df, "PET", args, pet_extractor, args.pet_params_yaml
        )

    if args.mode in ("ct", "both"):
        ct_extractor = featureextractor.RadiomicsFeatureExtractor(
            str(args.ct_params_yaml)
        )
        logging.info("CT PyRadiomics settings loaded successfully.")
        process_modality(
            df, "CT", args, ct_extractor, args.ct_params_yaml
        )

    logging.info("Objective 8D/8E script finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
