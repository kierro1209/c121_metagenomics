#!/usr/bin/env python3
"""Build the sequencing-record-level Accuracy manifest and metadata crosswalks.

This script intentionally stops at metadata assembly. It does not download or
audit sequence reads and does not adjudicate mNGS/RVP discrepancies.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests


ARTICLE_URL = "https://www.nature.com/articles/s41467-024-51470-y"
ZENODO_RECORD_URL = "https://zenodo.org/records/10553379"
ZENODO_API_URL = "https://zenodo.org/api/records/10553379"
SRA_RUNINFO_URL = (
    "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/"
    "runinfo?acc=PRJNA1084017"
)
ZENODO_MASTER_URL = (
    "https://zenodo.org/api/records/10553379/files/"
    "respmNGS-sequencing-mastersheet.xlsx/content"
)
SUPPLEMENTARY_INFO_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2F"
    "s41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM1_ESM.pdf"
)
SUPPLEMENTARY_DATA_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2F"
    "s41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM4_ESM.xlsx"
)
SOURCE_DATA_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2F"
    "s41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM6_ESM.xlsx"
)

ZENODO_MASTER_MD5 = "6c315aadb1c80d9433e8cdcd625e711a"
EXPECTED_ACCURACY_RECORDS = 281
EXPECTED_CLINICAL_SPECIMENS = 191

SOURCES = {
    "zenodo_master": {
        "filename": "respmNGS-sequencing-mastersheet.xlsx",
        "url": ZENODO_MASTER_URL,
        "landing_page": ZENODO_RECORD_URL,
        "role": "Sequencing record master sheet; Purpose filter and run identifiers",
        "expected_md5": ZENODO_MASTER_MD5,
    },
    "supplementary_data_1": {
        "filename": "41467_2024_51470_MOESM4_ESM.xlsx",
        "url": SUPPLEMENTARY_DATA_URL,
        "landing_page": ARTICLE_URL,
        "role": "Clinical and laboratory metadata (Supplementary Data 1)",
    },
    "supplementary_information": {
        "filename": "41467_2024_51470_MOESM1_ESM.pdf",
        "url": SUPPLEMENTARY_INFO_URL,
        "landing_page": ARTICLE_URL,
        "role": "Supplementary Tables 1-8, including discrepancy Tables 6-7",
    },
    "figure_source_data": {
        "filename": "41467_2024_51470_MOESM6_ESM.xlsx",
        "url": SOURCE_DATA_URL,
        "landing_page": ARTICLE_URL,
        "role": "Figure source data; contextual, not complete clinical metadata",
    },
    "sra_runinfo": {
        "filename": "PRJNA1084017_runinfo.csv",
        "url": SRA_RUNINFO_URL,
        "landing_page": "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1084017",
        "role": "SRA/BioSample accessions and read metadata for sequence retrieval",
    },
}


def digest(path: Path, algorithm: str) -> str:
    h = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    with requests.get(
        url,
        stream=True,
        timeout=(30, 180),
        headers={"User-Agent": "c121-metagenomics-metadata/1.0"},
    ) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                if chunk:
                    handle.write(chunk)
    temporary.replace(destination)


def acquire_sources(raw_dir: Path, force: bool) -> dict[str, dict[str, Any]]:
    provenance: dict[str, dict[str, Any]] = {}
    retrieved_at = datetime.now(timezone.utc).isoformat()
    for source_name, source in SOURCES.items():
        destination = raw_dir / str(source["filename"])
        if force or not destination.exists():
            download(str(source["url"]), destination)
        if source.get("expected_md5"):
            observed_md5 = digest(destination, "md5")
            if observed_md5 != source["expected_md5"]:
                raise ValueError(
                    f"{source_name} MD5 mismatch: {observed_md5} "
                    f"!= {source['expected_md5']}"
                )
        provenance[source_name] = {
            **source,
            "local_path": destination.as_posix(),
            "bytes": destination.stat().st_size,
            "sha256": digest(destination, "sha256"),
            "retrieved_or_verified_at_utc": retrieved_at,
        }
    return provenance


def clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def canonical_specimen_id(original_id: str) -> tuple[str, str]:
    """Return canonical key and mapping method without speculative joins."""
    original_id = clean_text(original_id)
    np_match = re.fullmatch(r"NP0*(\d+)", original_id, flags=re.IGNORECASE)
    if np_match:
        # NP is a laboratory identifier family, not a numeric alias for the
        # article's P-NP identifiers (for example, NP419 is not P-NP-419).
        return f"NP{int(np_match.group(1)):03d}", "NP_lab_id_article_crosswalk_unresolved"
    pnp_match = re.fullmatch(r"P-NP-0*(\d+)", original_id, flags=re.IGNORECASE)
    if pnp_match:
        return f"P-NP-{int(pnp_match.group(1))}", "P-NP_canonical"
    nnp_match = re.fullmatch(r"N-NP-0*(\d+)", original_id, flags=re.IGNORECASE)
    if nnp_match:
        return f"N-NP-{int(nnp_match.group(1))}", "N-NP_canonical"
    bal_match = re.fullmatch(r"BAL0*(\d+)", original_id, flags=re.IGNORECASE)
    if bal_match:
        return f"BAL{int(bal_match.group(1)):03d}", "BAL_canonical"
    nuc_match = re.fullmatch(r"NUC0*(\d+)", original_id, flags=re.IGNORECASE)
    if nuc_match:
        # NUC identifiers are valid specimen identifiers. The article workbook
        # uses P-NP aliases but provides no explicit NUC/P-NP crosswalk.
        return f"NUC{int(nuc_match.group(1)):05d}", "NUC_canonical_crosswalk_unresolved"
    if re.fullmatch(r"[PN]C-\d+", original_id, flags=re.IGNORECASE):
        return original_id.upper(), "technical_control_label"
    raise ValueError(f"Unrecognized Accuracy Original_ID format: {original_id!r}")


def article_id(value: Any) -> str:
    value = clean_text(value)
    if not value:
        return ""
    canonical, _ = canonical_specimen_id(value)
    return canonical


def record_role(original_id: str) -> str:
    if re.fullmatch(r"PC-\d+", original_id, flags=re.IGNORECASE):
        return "positive_control"
    if re.fullmatch(r"NC-\d+", original_id, flags=re.IGNORECASE):
        return "negative_control"
    return "clinical_sample"


def normalize_clinical_metadata(
    workbook: Path, output_dir: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_excel(
        workbook, sheet_name="Supplementary Data 1", dtype=str
    ).rename(columns=lambda column: clean_text(column))
    source.insert(0, "source_row", range(2, len(source) + 2))
    source["article_sample_id"] = source["Sample ID"].ffill().map(article_id)
    if (source["article_sample_id"] == "").any():
        rows = source.loc[source["article_sample_id"] == "", "source_row"].tolist()
        raise ValueError(f"Clinical continuation rows lack a preceding ID: {rows}")

    renamed = {
        "Sample Type": "sample_type",
        "PCR Results": "pcr_results",
        "mNGS Results": "mngs_results",
        "Viral Load (copies/ml)": "viral_load_copies_ml",
        "Clinical Diagnosis": "clinical_diagnosis",
        "Symptoms": "symptoms",
        "Discrepancy Testing (Ct value)": "discrepancy_testing_ct",
        "Comments": "comments",
    }
    severity_column = next(
        column for column in source.columns if column.startswith("Disease Severity")
    )
    source = source.rename(columns={**renamed, severity_column: "disease_severity"})
    line_columns = [
        "source_row",
        "article_sample_id",
        "sample_type",
        "pcr_results",
        "mngs_results",
        "viral_load_copies_ml",
        "clinical_diagnosis",
        "disease_severity",
        "symptoms",
        "discrepancy_testing_ct",
        "comments",
    ]
    lines = source[line_columns].copy()
    for column in line_columns[2:]:
        lines[column] = lines[column].map(clean_text)
    lines.to_csv(output_dir / "clinical_metadata_lines.csv", index=False)

    def combine(series: pd.Series) -> str:
        values = [clean_text(value) for value in series]
        return " || ".join(dict.fromkeys(value for value in values if value))

    samples = (
        lines.groupby("article_sample_id", sort=False)
        .agg(
            source_row_start=("source_row", "min"),
            source_row_end=("source_row", "max"),
            metadata_line_count=("source_row", "size"),
            sample_type=("sample_type", combine),
            pcr_results=("pcr_results", combine),
            mngs_results=("mngs_results", combine),
            viral_load_copies_ml=("viral_load_copies_ml", combine),
            clinical_diagnosis=("clinical_diagnosis", combine),
            disease_severity=("disease_severity", combine),
            symptoms=("symptoms", combine),
            discrepancy_testing_ct=("discrepancy_testing_ct", combine),
            comments=("comments", combine),
        )
        .reset_index()
    )
    samples.to_csv(output_dir / "clinical_metadata_samples.csv", index=False)
    return lines, samples


def build_manifest(
    master_workbook: Path,
    sra_runinfo_csv: Path,
    clinical_samples: pd.DataFrame,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    master = pd.read_excel(master_workbook, sheet_name="Sheet1", dtype=str)
    master = master.rename(columns=lambda column: clean_text(column))
    required = {
        "File_Name",
        "Sample_Name",
        "Instrument",
        "mNGS_run",
        "Sample_ID",
        "Purpose",
        "Original_ID",
        "Original_Run",
    }
    missing = required - set(master.columns)
    if missing:
        raise ValueError(f"Zenodo master sheet is missing columns: {sorted(missing)}")

    accuracy = master.loc[master["Purpose"].map(clean_text).eq("Accuracy")].copy()
    accuracy.insert(0, "master_row", accuracy.index + 2)
    accuracy = accuracy.reset_index(drop=True)
    if len(accuracy) != EXPECTED_ACCURACY_RECORDS:
        raise ValueError(
            f"Expected {EXPECTED_ACCURACY_RECORDS} Accuracy rows; got {len(accuracy)}"
        )
    if accuracy[list(required)].isna().any().any():
        raise ValueError("Accuracy rows contain missing required master-sheet values")
    if not accuracy["File_Name"].is_unique:
        raise ValueError("Accuracy File_Name values are not unique")
    if not accuracy["Sample_Name"].is_unique:
        raise ValueError("Accuracy Sample_Name values are not unique")

    canonical = accuracy["Original_ID"].map(canonical_specimen_id)
    accuracy["specimen_key"] = canonical.map(lambda pair: pair[0])
    accuracy["mapping_method"] = canonical.map(lambda pair: pair[1])
    accuracy["record_role"] = accuracy["Original_ID"].map(record_role)
    accuracy["is_control"] = accuracy["record_role"].ne("clinical_sample")
    accuracy["article_sample_id"] = accuracy["specimen_key"].where(
        accuracy["specimen_key"].str.match(r"^(?:P-NP-|N-NP-|BAL)"), ""
    )
    clinical_ids = set(clinical_samples["article_sample_id"])
    accuracy["clinical_metadata_join"] = "not_in_supplement"
    accuracy.loc[accuracy["is_control"], "clinical_metadata_join"] = "technical_control"
    accuracy.loc[
        (~accuracy["is_control"]) & accuracy["article_sample_id"].eq(""),
        "clinical_metadata_join",
    ] = "unresolved_NP_NUC_alias"
    accuracy.loc[
        accuracy["article_sample_id"].isin(clinical_ids), "clinical_metadata_join"
    ] = "matched"

    clinical_mask = ~accuracy["is_control"]
    accuracy["specimen_record_number"] = 0
    accuracy.loc[clinical_mask, "specimen_record_number"] = (
        accuracy.loc[clinical_mask].groupby("specimen_key").cumcount() + 1
    )
    accuracy["is_rerun_or_replicate"] = (
        clinical_mask & accuracy["specimen_record_number"].gt(1)
    )
    accuracy["is_explicit_repeat_run"] = accuracy["Original_Run"].str.contains(
        r"rpt|repeat", case=False, regex=True
    )

    initial_accuracy = master.loc[
        master["Purpose"].eq("Initial_Accuracy"),
        ["Original_ID", "File_Name", "Sample_Name", "Sample_ID", "mNGS_run"],
    ].drop_duplicates("Original_ID")
    initial_lookup = initial_accuracy.set_index("Original_ID")
    accuracy["linked_initial_accuracy_file"] = accuracy["Original_ID"].map(
        initial_lookup["File_Name"]
    ).fillna("")
    accuracy["linked_initial_accuracy_sample"] = accuracy["Original_ID"].map(
        initial_lookup["Sample_Name"]
    ).fillna("")

    runinfo = pd.read_csv(sra_runinfo_csv, dtype=str)
    sra_columns = [
        "SampleName",
        "Run",
        "Experiment",
        "BioSample",
        "spots",
        "bases",
        "avgLength",
        "size_MB",
        "download_path",
    ]
    missing_sra_columns = set(sra_columns) - set(runinfo.columns)
    if missing_sra_columns:
        raise ValueError(
            f"SRA runinfo is missing columns: {sorted(missing_sra_columns)}"
        )
    sra = runinfo[sra_columns].drop_duplicates("SampleName")
    accuracy = accuracy.merge(
        sra, how="left", left_on="Sample_Name", right_on="SampleName", validate="1:1"
    ).drop(columns="SampleName")
    if accuracy["Run"].isna().any():
        missing_names = accuracy.loc[accuracy["Run"].isna(), "Sample_Name"].tolist()
        raise ValueError(f"Accuracy records missing SRA joins: {missing_names}")
    primary_file = (
        accuracy.loc[clinical_mask]
        .drop_duplicates("specimen_key")
        .set_index("specimen_key")["File_Name"]
    )
    accuracy["primary_file_name"] = accuracy["specimen_key"].map(primary_file).fillna("")
    accuracy["fastq_stem"] = accuracy["File_Name"].str.replace(
        r"\.fastq\.gz$", "", regex=True
    )

    main_mask = clinical_mask & ~accuracy["is_explicit_repeat_run"]
    main_specimen_keys = accuracy.loc[main_mask, "specimen_key"].nunique()
    # The master structurally reconciles the article's 191 as 186 unambiguous
    # non-run-039 keys plus five members of the 20-record repeat run. The
    # published metadata do not identify which five, so no speculative collapse
    # is performed.
    if main_specimen_keys != 186:
        raise ValueError(
            f"Expected 186 unambiguous non-repeat-run specimen keys; "
            f"got {main_specimen_keys}"
        )

    output_columns = [
        "master_row",
        "File_Name",
        "fastq_stem",
        "Sample_Name",
        "Sample_ID",
        "Instrument",
        "mNGS_run",
        "Original_Run",
        "Purpose",
        "Original_ID",
        "specimen_key",
        "article_sample_id",
        "mapping_method",
        "record_role",
        "is_control",
        "specimen_record_number",
        "is_rerun_or_replicate",
        "primary_file_name",
        "clinical_metadata_join",
        "is_explicit_repeat_run",
        "linked_initial_accuracy_file",
        "linked_initial_accuracy_sample",
        "Run",
        "Experiment",
        "BioSample",
        "spots",
        "bases",
        "avgLength",
        "size_MB",
        "download_path",
    ]
    accuracy[output_columns].to_csv(output_dir / "accuracy_manifest.csv", index=False)

    specimens = (
        accuracy.loc[clinical_mask]
        .groupby("specimen_key", sort=False)
        .agg(
            article_sample_id=("article_sample_id", "first"),
            source_original_ids=("Original_ID", lambda x: " || ".join(dict.fromkeys(x))),
            sequencing_record_count=("File_Name", "size"),
            file_names=("File_Name", lambda x: " || ".join(x)),
            sample_names=("Sample_Name", lambda x: " || ".join(x)),
            sample_ids=("Sample_ID", lambda x: " || ".join(x)),
            instruments=("Instrument", lambda x: " || ".join(dict.fromkeys(x))),
            mngs_runs=("mNGS_run", lambda x: " || ".join(dict.fromkeys(x))),
            original_runs=("Original_Run", lambda x: " || ".join(dict.fromkeys(x))),
            clinical_metadata_join=("clinical_metadata_join", "first"),
            explicit_repeat_run=("is_explicit_repeat_run", "max"),
            sra_runs=("Run", lambda x: " || ".join(x)),
            biosamples=("BioSample", lambda x: " || ".join(x)),
        )
        .reset_index()
    )
    specimens.to_csv(output_dir / "accuracy_specimens.csv", index=False)

    identifier_map = (
        accuracy.groupby("Original_ID", sort=False)
        .agg(
            specimen_key=("specimen_key", "first"),
            article_sample_id=("article_sample_id", "first"),
            mapping_method=("mapping_method", "first"),
            record_role=("record_role", "first"),
            sequencing_record_count=("File_Name", "size"),
            clinical_metadata_join=("clinical_metadata_join", "first"),
        )
        .reset_index()
        .rename(columns={"Original_ID": "source_identifier"})
    )
    identifier_map.to_csv(output_dir / "accuracy_identifier_map.csv", index=False)
    return accuracy, specimens, identifier_map


def write_unresolved(
    identifier_map: pd.DataFrame,
    clinical_samples: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    unresolved = identifier_map.loc[
        identifier_map["mapping_method"].isin(
            [
                "NUC_canonical_crosswalk_unresolved",
                "NP_lab_id_article_crosswalk_unresolved",
            ]
        )
    ].copy()
    mapped_article_ids = set(
        identifier_map.loc[
            identifier_map["article_sample_id"].ne(""), "article_sample_id"
        ]
    )
    unmatched_article_ids = sorted(
        set(clinical_samples["article_sample_id"]) - mapped_article_ids
    )
    unresolved["candidate_article_id_pool"] = " || ".join(unmatched_article_ids)
    unresolved["issue"] = (
        "Published Zenodo and article supplements provide no explicit "
        "NP/NUC-to-P-NP crosswalk; no pairing was inferred."
    )
    unresolved.to_csv(output_dir / "unresolved_identifier_mappings.csv", index=False)
    return unresolved


def normalize_figure_context(workbook: Path, output_dir: Path) -> int:
    frame = pd.read_excel(
        workbook, sheet_name="Figure_2C_viral_loads_severity", dtype=str
    )
    frame.to_csv(output_dir / "figure2c_viral_load_context.csv", index=False)
    return len(frame)


def write_provenance(
    provenance: dict[str, dict[str, Any]],
    accuracy: pd.DataFrame,
    specimens: pd.DataFrame,
    identifier_map: pd.DataFrame,
    clinical_lines: pd.DataFrame,
    clinical_samples: pd.DataFrame,
    unresolved: pd.DataFrame,
    figure_rows: int,
    output_dir: Path,
) -> None:
    role_counts = Counter(accuracy["record_role"])
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "assemble-manifest only; no FASTQ audit and no final discrepancy report"
        ),
        "selection_rule": "Zenodo master sheet Purpose exactly equals 'Accuracy'",
        "counts": {
            "accuracy_sequencing_records": len(accuracy),
            "unique_fastq_file_names": accuracy["File_Name"].nunique(),
            "clinical_sequencing_records": int((~accuracy["is_control"]).sum()),
            "control_sequencing_records": int(accuracy["is_control"].sum()),
            "positive_control_records": role_counts["positive_control"],
            "negative_control_records": role_counts["negative_control"],
            "observed_clinical_specimen_keys_before_alias_resolution": len(specimens),
            "rerun_or_replicate_records": int(
                accuracy["is_rerun_or_replicate"].sum()
            ),
            "unique_master_original_ids": accuracy["Original_ID"].nunique(),
            "unique_normalized_identifiers": len(identifier_map),
            "supplementary_clinical_lines": len(clinical_lines),
            "supplementary_clinical_sample_ids": len(clinical_samples),
            "expected_article_clinical_specimens": EXPECTED_CLINICAL_SPECIMENS,
            "unambiguous_non_repeat_run_specimen_keys": int(
                accuracy.loc[
                    (~accuracy["is_control"]) & (~accuracy["is_explicit_repeat_run"]),
                    "specimen_key",
                ].nunique()
            ),
            "explicit_repeat_run_clinical_records": int(
                ((~accuracy["is_control"]) & accuracy["is_explicit_repeat_run"]).sum()
            ),
            "repeat_run_records_linked_to_initial_accuracy": int(
                (
                    (~accuracy["is_control"])
                    & accuracy["is_explicit_repeat_run"]
                    & accuracy["linked_initial_accuracy_file"].ne("")
                ).sum()
            ),
            "unresolved_NP_NUC_identifiers": len(unresolved),
            "figure2c_context_rows": figure_rows,
        },
        "known_mapping_limitations": [
            (
                "Supplementary Data 1 contains 109 explicit sample IDs across "
                "134 result lines, while the paper describes 110 virus-positive "
                "samples and the normalized sequencing master contains 191 "
                "clinical specimens."
            ),
            (
                f"{len(unresolved)} NP/NUC identifiers cannot be paired to the "
                "article's P-NP identifiers from the published files alone. "
                "They remain stable specimen keys and are not orphaned FASTQ rows."
            ),
            (
                "The 191 article-cohort count structurally reconciles as 186 "
                "unambiguous non-repeat-run specimen keys plus 5 members of the "
                "20 clinical records in explicit repeat run 039. The published "
                "metadata do not identify which 5; all 20 are retained and flagged."
            ),
        ],
        "sources": provenance,
    }
    (output_dir / "accuracy_manifest_provenance.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    markdown = [
        "# Accuracy manifest provenance",
        "",
        "Scope: metadata assembly only. No sequence audit or final adjudication was run.",
        "",
        "Selection: `Purpose == \"Accuracy\"` in the Zenodo sequencing master sheet.",
        "",
        "## Source URLs",
        "",
        f"- Article: {ARTICLE_URL}",
        f"- Zenodo record: {ZENODO_RECORD_URL}",
        f"- Zenodo API record: {ZENODO_API_URL}",
    ]
    for source in provenance.values():
        markdown.append(f"- {source['role']}: {source['url']}")
    markdown += [
        "",
        "## Verified counts",
        "",
        f"- Accuracy sequencing records: {len(accuracy)}",
        f"- Clinical sequencing records: {(~accuracy['is_control']).sum()}",
        f"- Control sequencing records: {accuracy['is_control'].sum()}",
        f"- Observed clinical specimen keys (before unresolved alias collapse): "
        f"{len(specimens)}",
        f"- Article clinical cohort: {EXPECTED_CLINICAL_SPECIMENS} "
        "(186 unambiguous keys + 5 unresolved members of repeat run 039)",
        f"- Rerun/replicate records beyond each specimen's first record: "
        f"{accuracy['is_rerun_or_replicate'].sum()}",
        f"- Unresolved NP/NUC-to-article identifier crosswalks: {len(unresolved)}",
        "",
        "See `accuracy_manifest_provenance.json` for hashes and machine-readable details.",
    ]
    (output_dir / "README_accuracy_manifest.md").write_text(
        "\n".join(markdown) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Project root (default: parent of analysis/)",
    )
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Redownload raw source files even if already present",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.project_root.resolve()
    raw_dir = root / "data" / "raw"
    output_dir = root / "data" / "normalized"
    raw_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    provenance = acquire_sources(raw_dir, args.force_download)
    clinical_lines, clinical_samples = normalize_clinical_metadata(
        raw_dir / str(SOURCES["supplementary_data_1"]["filename"]), output_dir
    )
    accuracy, specimens, identifier_map = build_manifest(
        raw_dir / str(SOURCES["zenodo_master"]["filename"]),
        raw_dir / str(SOURCES["sra_runinfo"]["filename"]),
        clinical_samples,
        output_dir,
    )
    unresolved = write_unresolved(identifier_map, clinical_samples, output_dir)
    figure_rows = normalize_figure_context(
        raw_dir / str(SOURCES["figure_source_data"]["filename"]), output_dir
    )
    write_provenance(
        provenance,
        accuracy,
        specimens,
        identifier_map,
        clinical_lines,
        clinical_samples,
        unresolved,
        figure_rows,
        output_dir,
    )
    print(
        json.dumps(
            {
                "accuracy_records": len(accuracy),
                "clinical_records": int((~accuracy["is_control"]).sum()),
                "control_records": int(accuracy["is_control"].sum()),
                "observed_clinical_specimen_keys": len(specimens),
                "expected_article_clinical_specimens": EXPECTED_CLINICAL_SPECIMENS,
                "rerun_or_replicate_records": int(
                    accuracy["is_rerun_or_replicate"].sum()
                ),
                "unresolved_NP_NUC_identifiers": len(unresolved),
                "output_directory": str(output_dir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
