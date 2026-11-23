#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PSMA Objective 9D/9E
Direct, resumable PET or CT PyRadiomics extraction from Objective 9A
PET-derived multi-label candidate maps.

Key design
----------
- Candidate geometry is NEVER reconstructed in Python.
- Each ROI is read directly from:
      candidate_labelmap_path + candidate_label_value
- PET and CT can be run simultaneously in separate PowerShell windows.
- Subject-level parallel processing avoids repeated full-volume loading.
- Cropped extraction reduces runtime and memory.
- Successful candidate IDs are skipped on rerun.
- Results and errors are appended periodically.
- CT is resampled to the PET/candidate-map grid only when geometry differs.
- Locked PSMA YAML parameter files are used.

Examples
--------
python PSMA_objective9D9E_extract_radiomics_DIRECT_RESUMABLE.py --modality PET --scope internal --workers 4
python PSMA_objective9D9E_extract_radiomics_DIRECT_RESUMABLE.py --modality CT  --scope internal --workers 4
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import logging
import os
import sys
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor

warnings.filterwarnings("ignore")
logging.getLogger("radiomics").setLevel(logging.ERROR)

# =========================================================================
# LOCKED PROJECT PATHS
# =========================================================================

SOURCE_PROJECT_ROOT = (Path(os.environ["PSMA_PROJECT_ROOT"]) / "PSMA_AI_Project")
V2_PROJECT_ROOT = (Path(os.environ["PSMA_PROJECT_ROOT"]) / "PSMA_AI_Project_V2_PET_Derived")

OBJECTIVE9A_ROOT = (
    V2_PROJECT_ROOT
    / "06_results"
    / "Objective9_PET_Derived_Candidate_Pipeline"
    / "Objective9A_PET_Only_Candidate_Generation"
)

OUTPUT_ROOT = (
    V2_PROJECT_ROOT
    / "06_results"
    / "Objective9_PET_Derived_Candidate_Pipeline"
    / "Objective9D9E_Direct_PET_CT_Radiomics"
)

PET_YAML = (Path(__file__).resolve().parents[1] / "config" / "PSMA_objective5C_PET_radiomics_params.yaml")

CT_YAML = (Path(__file__).resolve().parents[1] / "config" / "PSMA_objective5C_CT_radiomics_params.yaml")

# =========================================================================
# EXTRACTION SETTINGS
# =========================================================================

SAVE_EVERY_ROWS = 50
CROP_PAD_VOX = 3

REQUIRED_COLUMNS = [
    "candidate_global_id",
    "dataset_name",
    "cohort",
    "split",
    "case_id",
    "candidate_label_value",
    "candidate_volume_mL",
    "PET_path",
    "candidate_labelmap_path",
]

TEXT_COLUMNS = {
    "candidate_global_id",
    "dataset_name",
    "cohort",
    "split",
    "case_id",
    "PET_path",
    "candidate_labelmap_path",
}


# =========================================================================
# CSV AND TABLE UTILITIES
# =========================================================================

def normalize_column_name(value: str) -> str:
    return str(value).strip()


def read_objective9_csv(path: Path) -> pd.DataFrame:
    """Read comma-delimited Objective 9 CSV without delimiter guessing."""
    df = pd.read_csv(path, sep=",", dtype={c: "string" for c in TEXT_COLUMNS})
    df.columns = [normalize_column_name(c) for c in df.columns]
    return df


def find_latest_candidate_index(split_dir: Path) -> Path:
    finals = sorted(
        split_dir.glob("PSMA_objective9A_PET_only_candidate_index_*.csv"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    checkpoint = (
        split_dir
        / "checkpoints"
        / "PSMA_objective9A_candidate_index_checkpoint.csv"
    )

    candidates = finals + ([checkpoint] if checkpoint.exists() else [])
    diagnostics: List[str] = []

    for path in candidates:
        try:
            df = read_objective9_csv(path)
            missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
            if not missing and len(df) > 0:
                return path
            diagnostics.append(
                f"{path} | rows={len(df)} | missing={','.join(missing)}"
            )
        except Exception as exc:
            diagnostics.append(f"{path} | error={exc}")

    raise RuntimeError(
        "No valid Objective 9A candidate index found in "
        f"{split_dir}\n" + "\n".join(diagnostics)
    )


def append_rows_csv(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    if not rows:
        return

    path.parent.mkdir(parents=True, exist_ok=True)

    if path.exists() and path.stat().st_size > 0:
        # Preserve the already-locked output schema.
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            fieldnames = next(reader)
    else:
        fieldnames: List[str] = []
        for row in rows:
            for key in row.keys():
                if key not in fieldnames:
                    fieldnames.append(key)

    file_exists = path.exists() and path.stat().st_size > 0

    with path.open("a", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        if not file_exists:
            writer.writeheader()
        writer.writerows(rows)


def read_completed_ids(output_csv: Path) -> set[str]:
    if not output_csv.exists() or output_csv.stat().st_size == 0:
        return set()

    try:
        done = pd.read_csv(
            output_csv,
            sep=",",
            usecols=["candidate_global_id"],
            dtype={"candidate_global_id": "string"},
        )
        return set(done["candidate_global_id"].dropna().astype(str))
    except Exception:
        return set()



def normalize_filesystem_path(value: Any) -> Path:
    """
    Convert a table value to a clean filesystem Path.

    On Windows, absolute local and UNC paths are converted to extended-length
    paths (\\\\?\\...) so files beyond the traditional MAX_PATH limit can be
    accessed consistently by pathlib and SimpleITK.
    """
    raw = str(value).strip().strip('"').strip("'")
    if not raw:
        raise ValueError("Empty filesystem path")

    if os.name != "nt":
        return Path(raw).expanduser()

    absolute = os.path.abspath(os.path.expanduser(raw))

    if absolute.startswith("\\\\?\\"):
        return Path(absolute)

    if absolute.startswith("\\\\"):
        # UNC path: \\server\share\folder -> \\?\UNC\server\share\folder
        return Path("\\\\?\\UNC\\" + absolute[2:])

    # Local drive path: D:\folder -> \\?\D:\folder
    return Path("\\\\?\\" + absolute)


# =========================================================================
# IMAGE AND GEOMETRY UTILITIES
# =========================================================================

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
    interpolator: int,
    default_value: float,
) -> sitk.Image:
    return sitk.Resample(
        image,
        reference,
        sitk.Transform(),
        interpolator,
        default_value,
        image.GetPixelID(),
    )


def resolve_ct_path(row: pd.Series) -> Path:
    pet_path = normalize_filesystem_path(row["PET_path"])
    cohort = str(row["cohort"]).lower()

    if cohort == "internal":
        name = pet_path.name
        if "_0001.nii.gz" in name:
            ct_name = name.replace("_0001.nii.gz", "_0000.nii.gz")
        elif "_0001.nii" in name:
            ct_name = name.replace("_0001.nii", "_0000.nii")
        else:
            raise ValueError(
                f"Cannot derive internal CT path from PET filename: {pet_path}"
            )
        ct_path = pet_path.with_name(ct_name)
    elif cohort == "external":
        # Corrected DEEP-PSMA branch must explicitly use the PSMA subfolder.
        if pet_path.parent.name.upper() != "PSMA":
            raise ValueError(
                "External PET path is not inside the PSMA subfolder: "
                f"{pet_path}"
            )
        gz = pet_path.parent / "CT.nii.gz"
        nii = pet_path.parent / "CT.nii"
        ct_path = gz if gz.exists() else nii
    else:
        raise ValueError(f"Unknown cohort: {cohort}")

    if not ct_path.exists():
        raise FileNotFoundError(f"CT image not found: {ct_path}")

    return ct_path


def crop_to_label(
    image: sitk.Image,
    label_map: sitk.Image,
    label_value: int,
    pad: int,
) -> Tuple[sitk.Image, sitk.Image]:
    label_arr = sitk.GetArrayFromImage(label_map)
    mask = label_arr == int(label_value)

    coords = np.argwhere(mask)
    if coords.size == 0:
        raise ValueError(f"Candidate label {label_value} has zero voxels")

    z0, y0, x0 = coords.min(axis=0)
    z1, y1, x1 = coords.max(axis=0)

    size_xyz = label_map.GetSize()
    x0 = max(0, int(x0) - pad)
    y0 = max(0, int(y0) - pad)
    z0 = max(0, int(z0) - pad)
    x1 = min(size_xyz[0] - 1, int(x1) + pad)
    y1 = min(size_xyz[1] - 1, int(y1) + pad)
    z1 = min(size_xyz[2] - 1, int(z1) + pad)

    index = [x0, y0, z0]
    size = [x1 - x0 + 1, y1 - y0 + 1, z1 - z0 + 1]

    cropped_image = sitk.RegionOfInterest(image, size=size, index=index)
    cropped_labels = sitk.RegionOfInterest(label_map, size=size, index=index)

    binary_mask = sitk.Cast(cropped_labels == int(label_value), sitk.sitkUInt8)
    binary_mask.CopyInformation(cropped_labels)

    return cropped_image, binary_mask


def clean_feature_dictionary(
    values: Dict[str, Any],
    modality: str,
) -> Dict[str, Any]:
    output: Dict[str, Any] = {}

    for key, value in values.items():
        key = str(key)
        if key.startswith("diagnostics_"):
            continue

        output_key = f"{modality}_{key}"

        if isinstance(value, np.ndarray):
            if value.size != 1:
                continue
            value = value.ravel()[0]

        if isinstance(value, np.generic):
            value = value.item()

        if isinstance(value, (float, int, np.floating, np.integer)):
            value = float(value)
            if not np.isfinite(value):
                value = ""

        output[output_key] = value

    return output


# =========================================================================
# WORKER
# =========================================================================

def make_extractor(yaml_path: Path) -> featureextractor.RadiomicsFeatureExtractor:
    if not yaml_path.exists():
        raise FileNotFoundError(f"Locked PyRadiomics YAML not found: {yaml_path}")
    extractor = featureextractor.RadiomicsFeatureExtractor(str(yaml_path))
    return extractor


def process_case(
    payload: Tuple[str, str, List[Dict[str, Any]], str, str]
) -> Tuple[str, List[Dict[str, Any]], List[Dict[str, Any]], float]:
    case_id, modality, records, pet_yaml_str, ct_yaml_str = payload
    started = time.time()
    successes: List[Dict[str, Any]] = []
    errors: List[Dict[str, Any]] = []

    try:
        case_df = pd.DataFrame(records)
        first = case_df.iloc[0]

        pet_path = normalize_filesystem_path(first["PET_path"])
        labelmap_path = normalize_filesystem_path(
            first["candidate_labelmap_path"]
        )

        if not pet_path.exists():
            raise FileNotFoundError(
                f"PET image not found: {pet_path} "
                f"(characters={len(str(pet_path))})"
            )

        if not labelmap_path.exists():
            raise FileNotFoundError(
                f"Candidate label map not found: {labelmap_path} "
                f"(characters={len(str(labelmap_path))})"
            )

        pet_image = sitk.ReadImage(str(pet_path))
        label_map = sitk.ReadImage(str(labelmap_path))

        if not same_geometry(pet_image, label_map):
            raise ValueError(
                "PET and Objective 9A candidate label-map geometries differ"
            )

        modality_upper = modality.upper()

        if modality_upper == "PET":
            source_image = sitk.Cast(pet_image, sitk.sitkFloat32)
            extractor = make_extractor(Path(pet_yaml_str))
            image_path = pet_path
            geometry_resampled = False

        elif modality_upper == "CT":
            ct_path = resolve_ct_path(first)
            ct_image = sitk.ReadImage(str(ct_path))
            geometry_resampled = not same_geometry(ct_image, pet_image)

            if geometry_resampled:
                ct_image = resample_to_reference(
                    image=ct_image,
                    reference=pet_image,
                    interpolator=sitk.sitkBSpline,
                    default_value=-1024.0,
                )

            source_image = sitk.Cast(ct_image, sitk.sitkFloat32)
            extractor = make_extractor(Path(ct_yaml_str))
            image_path = ct_path

        else:
            raise ValueError(f"Unsupported modality: {modality}")

        for _, row in case_df.iterrows():
            candidate_id = str(row["candidate_global_id"])
            label_value = int(float(row["candidate_label_value"]))

            try:
                cropped_image, cropped_mask = crop_to_label(
                    source_image,
                    label_map,
                    label_value,
                    CROP_PAD_VOX,
                )

                features = extractor.execute(
                    cropped_image,
                    cropped_mask,
                    label=1,
                )

                result = row.to_dict()
                for key, value in list(result.items()):
                    if pd.isna(value):
                        result[key] = ""

                result.update(
                    {
                        "radiomics_status": "success",
                        "radiomics_modality": modality_upper,
                        "radiomics_image_path": str(image_path),
                        "candidate_label_value_verified": label_value,
                        "CT_resampled_to_PET_grid": (
                            int(geometry_resampled)
                            if modality_upper == "CT"
                            else ""
                        ),
                    }
                )
                result.update(clean_feature_dictionary(features, modality_upper))
                successes.append(result)

            except Exception as exc:
                errors.append(
                    {
                        "case_id": case_id,
                        "candidate_global_id": candidate_id,
                        "candidate_label_value": label_value,
                        "modality": modality_upper,
                        "error_message": str(exc),
                        "traceback": traceback.format_exc(limit=3),
                    }
                )

        del source_image, label_map, pet_image
        gc.collect()

    except Exception as exc:
        for record in records:
            errors.append(
                {
                    "case_id": case_id,
                    "candidate_global_id": str(
                        record.get("candidate_global_id", "")
                    ),
                    "candidate_label_value": record.get(
                        "candidate_label_value", ""
                    ),
                    "modality": modality.upper(),
                    "error_message": f"CASE_FAILURE: {exc}",
                    "traceback": traceback.format_exc(limit=5),
                }
            )

    return case_id, successes, errors, time.time() - started


# =========================================================================
# SCOPE AND OUTPUT LOGIC
# =========================================================================

def selected_split_directories(scope: str) -> List[Path]:
    internal = [
        OBJECTIVE9A_ROOT / "internal" / "train",
        OBJECTIVE9A_ROOT / "internal" / "validation",
        OBJECTIVE9A_ROOT / "internal" / "test",
    ]
    external = [OBJECTIVE9A_ROOT / "external" / "external"]

    if scope == "internal":
        return internal
    if scope == "external":
        return external
    if scope == "all":
        return internal + external
    raise ValueError(f"Unknown scope: {scope}")


def load_scope_candidates(scope: str) -> Tuple[pd.DataFrame, List[str]]:
    frames: List[pd.DataFrame] = []
    source_files: List[str] = []

    for split_dir in selected_split_directories(scope):
        if not split_dir.exists():
            raise FileNotFoundError(
                f"Objective 9A split folder not found: {split_dir}"
            )

        index_path = find_latest_candidate_index(split_dir)
        frame = read_objective9_csv(index_path)

        missing = [c for c in REQUIRED_COLUMNS if c not in frame.columns]
        if missing:
            raise ValueError(
                f"Missing columns in {index_path}: {', '.join(missing)}"
            )

        frames.append(frame)
        source_files.append(str(index_path))

    combined = pd.concat(frames, ignore_index=True)

    if combined["candidate_global_id"].duplicated().any():
        duplicates = combined.loc[
            combined["candidate_global_id"].duplicated(),
            "candidate_global_id",
        ].head(10)
        raise ValueError(
            "Duplicate candidate IDs detected: "
            + ", ".join(duplicates.astype(str))
        )

    return combined, source_files


def output_paths(
    modality: str,
    scope: str,
) -> Tuple[Path, Path, Path, Path]:
    run_dir = OUTPUT_ROOT / scope / modality.upper()
    run_dir.mkdir(parents=True, exist_ok=True)

    output_csv = (
        run_dir
        / f"PSMA_objective9{modality.upper()}_"
        f"{scope}_{modality.upper()}_radiomics_checkpoint.csv"
    )
    error_csv = (
        run_dir
        / f"PSMA_objective9{modality.upper()}_"
        f"{scope}_{modality.upper()}_radiomics_failures.csv"
    )
    state_json = (
        run_dir
        / f"PSMA_objective9{modality.upper()}_"
        f"{scope}_{modality.upper()}_state.json"
    )
    log_path = (
        run_dir
        / f"PSMA_objective9{modality.upper()}_"
        f"{scope}_{modality.upper()}_radiomics.log"
    )

    return output_csv, error_csv, state_json, log_path


def write_state(
    state_json: Path,
    modality: str,
    scope: str,
    input_rows: int,
    already_done: int,
    new_success: int,
    new_failures: int,
    source_files: Sequence[str],
    started: float,
) -> None:
    payload = {
        "objective": "9D/9E",
        "modality": modality.upper(),
        "scope": scope,
        "input_candidates": input_rows,
        "already_successful_at_start": already_done,
        "new_successful_this_run": new_success,
        "new_failures_this_run": new_failures,
        "current_success_total": already_done + new_success,
        "elapsed_minutes": round((time.time() - started) / 60.0, 3),
        "source_candidate_tables": list(source_files),
        "PET_yaml": str(PET_YAML),
        "CT_yaml": str(CT_YAML),
        "candidate_geometry_source": "Objective 9A PET-derived label maps",
        "updated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    state_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")


# =========================================================================
# MAIN
# =========================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Direct, resumable PSMA Objective 9 PET/CT PyRadiomics extraction"
        )
    )
    parser.add_argument(
        "--modality",
        required=True,
        choices=["PET", "CT", "pet", "ct"],
    )
    parser.add_argument(
        "--scope",
        default="internal",
        choices=["internal", "external", "all"],
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, min(4, os.cpu_count() or 1)),
    )
    parser.add_argument(
        "--retry-failures",
        action="store_true",
        help=(
            "Successful candidates remain skipped. Existing failure rows "
            "do not block retry."
        ),
    )
    args = parser.parse_args()

    modality = args.modality.upper()
    scope = args.scope.lower()
    workers = max(1, int(args.workers))

    if modality == "PET" and not PET_YAML.exists():
        raise FileNotFoundError(f"PET YAML not found: {PET_YAML}")
    if modality == "CT" and not CT_YAML.exists():
        raise FileNotFoundError(f"CT YAML not found: {CT_YAML}")

    candidates, source_files = load_scope_candidates(scope)

    output_csv, error_csv, state_json, log_path = output_paths(
        modality,
        scope,
    )

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_path, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )
    logger = logging.getLogger("objective9_radiomics")

    started = time.time()

    logger.info("=" * 72)
    logger.info("PSMA OBJECTIVE 9D/9E DIRECT RADIOMICS EXTRACTION")
    logger.info("Modality: %s | Scope: %s | Workers: %d", modality, scope, workers)
    logger.info("Candidate geometry: Objective 9A label maps")
    logger.info("Candidate reconstruction in Python: NO")
    logger.info("Input candidates: %d", len(candidates))
    logger.info("PET YAML: %s", PET_YAML)
    logger.info("CT YAML: %s", CT_YAML)
    logger.info("=" * 72)

    done_ids = read_completed_ids(output_csv)
    already_done = len(done_ids)

    todo = candidates.loc[
        ~candidates["candidate_global_id"].astype(str).isin(done_ids)
    ].copy()

    logger.info("Already successful: %d", already_done)
    logger.info("Candidates remaining: %d", len(todo))

    if todo.empty:
        logger.info("Nothing to do.")
        write_state(
            state_json,
            modality,
            scope,
            len(candidates),
            already_done,
            0,
            0,
            source_files,
            started,
        )
        return 0

    groups: List[Tuple[str, str, List[Dict[str, Any]], str, str]] = []
    for case_id, group in todo.groupby("case_id", sort=False):
        groups.append(
            (
                str(case_id),
                modality,
                group.to_dict(orient="records"),
                str(PET_YAML),
                str(CT_YAML),
            )
        )

    logger.info("Cases remaining: %d", len(groups))

    result_buffer: List[Dict[str, Any]] = []
    error_buffer: List[Dict[str, Any]] = []
    new_success = 0
    new_failures = 0
    completed_cases = 0

    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(process_case, payload): payload[0]
            for payload in groups
        }

        for future in as_completed(future_map):
            case_id = future_map[future]
            completed_cases += 1

            try:
                _, successes, errors, elapsed = future.result()
            except Exception as exc:
                successes = []
                errors = [
                    {
                        "case_id": case_id,
                        "candidate_global_id": "",
                        "candidate_label_value": "",
                        "modality": modality,
                        "error_message": f"FUTURE_FAILURE: {exc}",
                        "traceback": traceback.format_exc(limit=5),
                    }
                ]
                elapsed = 0.0

            result_buffer.extend(successes)
            error_buffer.extend(errors)
            new_success += len(successes)
            new_failures += len(errors)

            logger.info(
                "[%d/%d] %s | success=%d failures=%d time=%.1fs "
                "| cumulative success=%d failures=%d",
                completed_cases,
                len(groups),
                case_id,
                len(successes),
                len(errors),
                elapsed,
                new_success,
                new_failures,
            )

            if len(result_buffer) + len(error_buffer) >= SAVE_EVERY_ROWS:
                if result_buffer:
                    append_rows_csv(output_csv, result_buffer)
                    result_buffer.clear()
                if error_buffer:
                    append_rows_csv(error_csv, error_buffer)
                    error_buffer.clear()

                write_state(
                    state_json,
                    modality,
                    scope,
                    len(candidates),
                    already_done,
                    new_success,
                    new_failures,
                    source_files,
                    started,
                )
                logger.info("Progress saved.")

    if result_buffer:
        append_rows_csv(output_csv, result_buffer)
    if error_buffer:
        append_rows_csv(error_csv, error_buffer)

    write_state(
        state_json,
        modality,
        scope,
        len(candidates),
        already_done,
        new_success,
        new_failures,
        source_files,
        started,
    )

    logger.info("=" * 72)
    logger.info("OBJECTIVE 9%s COMPLETED/CURRENT STATE", modality)
    logger.info(
        "Success total: %d/%d",
        already_done + new_success,
        len(candidates),
    )
    logger.info("New failures this run: %d", new_failures)
    logger.info("Output: %s", output_csv)
    logger.info("Failures: %s", error_csv)
    logger.info("State: %s", state_json)
    logger.info("=" * 72)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
