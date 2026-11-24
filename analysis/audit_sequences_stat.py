#!/usr/bin/env python3
"""Independently audit Accuracy runs via NCBI STAT public taxonomy API.

Uses https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy for each SRR.
Applies a >=3 total_count threshold as a conservative proxy for the paper's
>=3 non-overlapping read/contig rule. STAT reports k-mer/spot taxonomy counts,
not SURPI+ non-overlapping alignments; evidence is labeled accordingly.
"""

from __future__ import annotations

import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/normalized/accuracy_manifest.csv"
CLINICAL = ROOT / "data/normalized/clinical_metadata_samples.csv"
TABLE6 = ROOT / "data/normalized/supplementary_table6_discrepancies.csv"
TABLE7 = ROOT / "data/normalized/supplementary_table7_discrepancies.csv"
OUT_DIR = ROOT / "data/normalized/stat"
CACHE_DIR = OUT_DIR / "run_taxonomy_cache"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)

API = "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy"
THRESHOLD = 3  # proxy for >=3 non-overlapping reads/contigs
WORKERS = 12

# Target families used for clinical RVP comparison and control checks.
TARGET_PATTERNS: dict[str, list[str]] = {
    "RSV": [
        r"\borthopneumovirus\b",
        r"respiratory syncytial",
        r"\bhuman orthopneumovirus\b",
        r"\brsv\b",
    ],
    "Influenza_A": [r"influenza a\b", r"\bh1n1\b", r"\bh3n2\b", r"alphainfluenzavirus"],
    "Influenza_B": [r"influenza b\b", r"betainfluenzavirus"],
    "SARS_CoV_2": [
        r"severe acute respiratory syndrome coronavirus 2",
        r"\bsars-cov-2\b",
        r"\bsarbecovirus\b",
    ],
    "Rhinovirus_Enterovirus": [
        r"\brhinovirus\b",
        r"\benterovirus\b",
        r"\bev-d68\b",
        r"\bevd68\b",
    ],
    "Adenovirus": [r"\badenovirus\b", r"\bmastadenovirus\b"],
    "Metapneumovirus": [r"\bmetapneumovirus\b"],
    "Parainfluenza": [
        r"\bparainfluenza\b",
        r"\brespirovirus\b",
        r"human respirovirus",
        r"\brubulavirus\b",
    ],
    "Coronavirus_OC43": [r"coronavirus oc43", r"betacoronavirus 1", r"\boc43\b"],
    "Coronavirus_HKU1": [r"coronavirus hku1", r"\bhku1\b"],
    "Coronavirus_NL63": [r"coronavirus nl63", r"\bnl63\b"],
    "Bocavirus": [r"\bbocavirus\b", r"\bbocaparvovirus\b"],
    "MS2_control": [r"emesvirus zinderi", r"\bescherichia virus ms2\b", r"\bms2\b"],
}


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat()


def fetch_taxonomy(acc: str, session: requests.Session) -> dict[str, Any]:
    cache = CACHE_DIR / f"{acc}.json"
    if cache.exists() and cache.stat().st_size > 10:
        return json.loads(cache.read_text())
    url = f"{API}?acc={acc}&cluster_name=public"
    for attempt in range(5):
        try:
            r = session.get(url, timeout=180)
            if r.status_code in {429, 500, 502, 503, 504}:
                time.sleep(2**attempt)
                continue
            r.raise_for_status()
            data = r.json()
            obj = data[0] if isinstance(data, list) else data
            cache.write_text(json.dumps(obj))
            return obj
        except Exception:
            if attempt == 4:
                raise
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to fetch {acc}")


def match_targets(org: str) -> list[str]:
    name = (org or "").lower()
    hits = []
    for target, patterns in TARGET_PATTERNS.items():
        if any(re.search(p, name) for p in patterns):
            hits.append(target)
    return hits


def summarize_run(acc: str, obj: dict[str, Any]) -> dict[str, Any]:
    table = obj.get("tax_table") or []
    totals = obj.get("tax_totals") or {}
    virus_node = next((r for r in table if int(r.get("tax_id", -1)) == 10239), None)
    target_hits: dict[str, dict[str, Any]] = {}
    viral_rows = []
    for row in table:
        org = row.get("org") or ""
        self_c = int(row.get("self_count") or 0)
        total_c = int(row.get("total_count") or 0)
        if total_c < 1 and self_c < 1:
            continue
        targets = match_targets(org)
        if not targets and "virus" not in org.lower() and "viridae" not in org.lower():
            continue
        viral_rows.append(
            {
                "tax_id": row.get("tax_id"),
                "org": org,
                "self_count": self_c,
                "total_count": total_c,
                "targets": targets,
            }
        )
        for t in targets:
            prev = target_hits.get(t)
            if prev is None or total_c > prev["total_count"]:
                target_hits[t] = {
                    "org": org,
                    "tax_id": row.get("tax_id"),
                    "self_count": self_c,
                    "total_count": total_c,
                    "pass_threshold": total_c >= THRESHOLD,
                }

    passing = sorted(
        [t for t, v in target_hits.items() if v["pass_threshold"] and t != "MS2_control"]
    )
    top_virus = sorted(viral_rows, key=lambda x: -x["total_count"])[:15]
    return {
        "Run": acc,
        "stat_available": True,
        "stat_version": totals.get("version"),
        "stat_analysed_spots": totals.get("analysed"),
        "stat_identified_spots": totals.get("identified"),
        "stat_total_spots": totals.get("total"),
        "virus_total_count": int((virus_node or {}).get("total_count") or 0),
        "virus_self_count": int((virus_node or {}).get("self_count") or 0),
        "passing_targets": ";".join(passing),
        "n_passing_targets": len(passing),
        "ms2_total_count": int((target_hits.get("MS2_control") or {}).get("total_count") or 0),
        "ms2_pass": bool((target_hits.get("MS2_control") or {}).get("pass_threshold")),
        "target_hit_details_json": json.dumps(target_hits, sort_keys=True),
        "top_virus_taxa_json": json.dumps(top_virus, sort_keys=True),
        "evidence_type": "sequence_STAT_kmer_counts",
        "threshold_rule": f">= {THRESHOLD} STAT total_count (proxy for >=3 non-overlapping reads/contigs)",
    }


def expected_targets_from_text(text: str) -> set[str]:
    if not isinstance(text, str) or not text.strip() or text.strip() in {"-", "NA", "nan"}:
        return set()
    t = text.lower()
    out: set[str] = set()
    if re.search(r"\brsv\b|orthopneumovirus|respiratory syncytial", t):
        out.add("RSV")
    if re.search(r"flu\s*a|influenza a|h3n2|h1n1", t):
        out.add("Influenza_A")
    if re.search(r"flu\s*b|influenza b", t):
        out.add("Influenza_B")
    if re.search(r"sars-cov-2|sars.?cov.?2|covid", t):
        out.add("SARS_CoV_2")
    if re.search(r"rhino|enterovirus|evd68|ev-d68", t):
        out.add("Rhinovirus_Enterovirus")
    if re.search(r"adeno|mastadenovirus", t):
        out.add("Adenovirus")
    if re.search(r"metapneumo|\bmpv\b", t):
        out.add("Metapneumovirus")
    if re.search(r"parainfluenza|respirovirus", t):
        out.add("Parainfluenza")
    if re.search(r"oc43", t):
        out.add("Coronavirus_OC43")
    if re.search(r"hku1", t):
        out.add("Coronavirus_HKU1")
    if re.search(r"nl63", t):
        out.add("Coronavirus_NL63")
    if re.search(r"bocavirus|bocaparvo", t):
        out.add("Bocavirus")
    # Explicit not-detected phrases remove that target from expected positives.
    if re.search(r"respiratory virus not detected", t):
        out.clear()
    return out


def classify_row(row: pd.Series) -> dict[str, Any]:
    role = row["record_role"]
    passing = set([x for x in str(row.get("passing_targets") or "").split(";") if x])
    notes: list[str] = []
    sequence_vs_metadata = "not_assessed"
    inconsistency_flags: list[str] = []
    confidence = "n/a"
    audit_call = "unknown"

    if not row.get("stat_available"):
        return {
            "audit_call": "stat_unavailable",
            "sequence_vs_metadata": "unavailable",
            "inconsistency_flags": "stat_unavailable",
            "audit_notes": "No NCBI STAT taxonomy payload for this SRR.",
            "confidence": "low",
        }

    if role == "positive_control":
        expected = {"SARS_CoV_2", "Influenza_A", "Influenza_B", "RSV"}
        missing = sorted(expected - passing)
        unexpected = sorted(passing - expected - {"MS2_control"})
        if not missing:
            sequence_vs_metadata = "consistent"
            audit_call = "pc_targets_detected"
            confidence = "high"
            notes.append("Positive-control Accuplex targets detected at threshold.")
        elif len(missing) <= 1 and (row.get("virus_total_count") or 0) >= THRESHOLD:
            sequence_vs_metadata = "partial"
            audit_call = "pc_partial"
            confidence = "medium"
            notes.append(f"Positive control missing threshold calls for: {', '.join(missing)}.")
            inconsistency_flags.append("pc_missing_target")
        else:
            sequence_vs_metadata = "conflict"
            audit_call = "pc_unexpected"
            confidence = "high"
            notes.append(f"Positive control missing: {', '.join(missing) or 'n/a'}.")
            inconsistency_flags.append("pc_missing_target")
        if unexpected:
            notes.append(f"Additional non-panel targets: {', '.join(unexpected)}.")
    elif role == "negative_control":
        pathogen_hits = sorted(passing)
        if not pathogen_hits:
            sequence_vs_metadata = "consistent"
            audit_call = "nc_clean"
            confidence = "high"
            notes.append("Negative control has no threshold respiratory-pathogen STAT hits.")
        else:
            sequence_vs_metadata = "conflict"
            audit_call = "nc_contamination_or_mislabeled"
            confidence = "medium"
            notes.append(
                "Negative control has threshold STAT hits for: " + ", ".join(pathogen_hits)
            )
            inconsistency_flags.append("nc_pathogen_signal")
    else:
        expected_mngs = expected_targets_from_text(str(row.get("mngs_results") or ""))
        expected_pcr = expected_targets_from_text(str(row.get("pcr_results") or ""))
        declared = expected_mngs or expected_pcr
        not_detected_declared = (
            isinstance(row.get("mngs_results"), str)
            and "not detected" in row["mngs_results"].lower()
            and not expected_mngs
        )
        if pd.isna(row.get("article_sample_id")) or (
            not declared and not not_detected_declared
        ):
            if passing:
                sequence_vs_metadata = "sequence_only_no_metadata_target"
                audit_call = "clinical_sequence_signal_without_joined_metadata"
                confidence = "low"
                notes.append(
                    "STAT pathogen signal present, but no joined Supp Data 1 target list "
                    "for this specimen key (often unresolved NP/NUC alias)."
                )
            else:
                sequence_vs_metadata = "sequence_only_no_metadata_target"
                audit_call = "clinical_no_threshold_pathogen"
                confidence = "low"
                notes.append(
                    "No threshold pathogen STAT hits; clinical metadata join incomplete."
                )
        else:
            missing = sorted(declared - passing)
            extra = sorted(passing - declared)
            if not_detected_declared and not passing:
                sequence_vs_metadata = "consistent"
                audit_call = "clinical_nd_supported"
                confidence = "medium"
                notes.append(
                    "Declared mNGS not-detected and no threshold STAT pathogen hits."
                )
            elif not_detected_declared and passing:
                sequence_vs_metadata = "conflict"
                audit_call = "clinical_nd_but_stat_positive"
                confidence = "medium"
                notes.append(
                    "Supp Data 1 reports respiratory virus not detected, but STAT has: "
                    + ", ".join(passing)
                )
                inconsistency_flags.append("metadata_nd_vs_stat_positive")
            elif not missing:
                sequence_vs_metadata = "consistent"
                audit_call = "clinical_targets_supported"
                confidence = "high"
                notes.append(
                    "Declared clinical respiratory targets supported by STAT threshold hits."
                )
                if extra:
                    notes.append(
                        "Additional STAT targets beyond declared list (may be co-detections "
                        f"or non-RVP organisms): {', '.join(extra)}."
                    )
            elif missing and (set(passing) & declared):
                sequence_vs_metadata = "partial"
                audit_call = "clinical_partial_support"
                confidence = "medium"
                notes.append(
                    f"Partial support; missing STAT threshold for: {', '.join(missing)}."
                )
                if extra:
                    notes.append(f"Extra STAT targets: {', '.join(extra)}.")
                inconsistency_flags.append("partial_target_support")
            else:
                sequence_vs_metadata = "conflict"
                audit_call = "clinical_target_mismatch"
                confidence = "medium"
                notes.append(
                    f"Declared targets {sorted(declared)} not supported by STAT; "
                    f"STAT passing={sorted(passing)}."
                )
                inconsistency_flags.append("declared_vs_stat_mismatch")

    disc_class = row.get("article_discrepancy_class")
    if isinstance(disc_class, str) and disc_class.strip():
        inconsistency_flags.append("article_reported_mngs_rvp_discrepancy")
        notes.append(f"Article discrepancy class: {disc_class}.")

    return {
        "audit_call": audit_call,
        "sequence_vs_metadata": sequence_vs_metadata,
        "inconsistency_flags": ";".join(sorted(set(inconsistency_flags)))
        if inconsistency_flags
        else "",
        "audit_notes": " ".join(notes),
        "confidence": confidence,
    }


def build_discrepancy_index() -> pd.DataFrame:
    frames = []
    for path, label in [(TABLE6, "mNGS_positive_RVP_negative"), (TABLE7, "mNGS_negative_RVP_positive")]:
        if path.exists():
            df = pd.read_csv(path)
            df["article_discrepancy_class"] = label
            frames.append(df)
    if not frames:
        return pd.DataFrame(columns=["article_sample_id", "article_discrepancy_class"])
    d = pd.concat(frames, ignore_index=True)
    # Keep one row per sample id with concatenated explanations.
    d = d[d["article_sample_id"].notna() & (d["article_sample_id"].astype(str).str.len() > 0)]
    grouped = (
        d.groupby("article_sample_id", as_index=False)
        .agg(
            article_discrepancy_class=("article_discrepancy_class", lambda s: ";".join(sorted(set(s)))),
            article_final_mngs_by_dtca=("final_mngs_by_dtca", lambda s: " || ".join(map(str, s))),
            article_id_mapping_status=("id_mapping_status", lambda s: ";".join(sorted(set(map(str, s))))),
            article_explanation=("explanation", lambda s: " || ".join(map(str, s))),
        )
    )
    return grouped


def main() -> None:
    manifest = pd.read_csv(MANIFEST)
    clinical = pd.read_csv(CLINICAL)
    disc = build_discrepancy_index()

    session_local = requests.Session()
    session_local.headers.update({"User-Agent": "metagenomics-audit/1.0"})

    runs = manifest["Run"].dropna().unique().tolist()
    print(f"fetching STAT for {len(runs)} runs...", flush=True)

    summaries: list[dict[str, Any]] = []
    errors: list[dict[str, str]] = []

    def worker(acc: str) -> tuple[str, dict[str, Any] | None, str | None]:
        sess = requests.Session()
        sess.headers.update({"User-Agent": "metagenomics-audit/1.0"})
        try:
            obj = fetch_taxonomy(acc, sess)
            return acc, summarize_run(acc, obj), None
        except Exception as e:  # noqa: BLE001
            return acc, None, str(e)

    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = {ex.submit(worker, acc): acc for acc in runs}
        for fut in as_completed(futs):
            acc, summary, err = fut.result()
            done += 1
            if err:
                errors.append({"Run": acc, "error": err})
                summaries.append(
                    {
                        "Run": acc,
                        "stat_available": False,
                        "passing_targets": "",
                        "n_passing_targets": 0,
                        "evidence_type": "unavailable",
                        "threshold_rule": f">= {THRESHOLD} STAT total_count",
                        "error": err,
                    }
                )
            else:
                assert summary is not None
                summaries.append(summary)
            if done % 25 == 0 or done == len(runs):
                print(f"  {done}/{len(runs)}", flush=True)

    stat_df = pd.DataFrame(summaries)
    merged = manifest.merge(stat_df, on="Run", how="left")
    merged = merged.merge(clinical, on="article_sample_id", how="left", suffixes=("", "_clin"))
    merged = merged.merge(disc, on="article_sample_id", how="left")

    audit_cols = []
    for _, row in merged.iterrows():
        audit_cols.append(classify_row(row))
    audit_df = pd.concat([merged.reset_index(drop=True), pd.DataFrame(audit_cols)], axis=1)

    # Persist detailed and compact outputs
    detail_path = OUT_DIR / "accuracy_stat_audit_detail.csv"
    compact_path = OUT_DIR / "accuracy_stat_audit_summary.csv"
    errors_path = OUT_DIR / "accuracy_stat_fetch_errors.csv"
    provenance_path = OUT_DIR / "accuracy_stat_audit_provenance.json"

    detail_cols = [
        "File_Name",
        "Run",
        "BioSample",
        "Original_ID",
        "specimen_key",
        "article_sample_id",
        "record_role",
        "is_control",
        "is_rerun_or_replicate",
        "is_explicit_repeat_run",
        "mapping_method",
        "clinical_metadata_join",
        "pcr_results",
        "mngs_results",
        "spots",
        "size_MB",
        "stat_available",
        "stat_version",
        "stat_analysed_spots",
        "stat_identified_spots",
        "virus_total_count",
        "ms2_total_count",
        "ms2_pass",
        "passing_targets",
        "n_passing_targets",
        "evidence_type",
        "threshold_rule",
        "audit_call",
        "sequence_vs_metadata",
        "inconsistency_flags",
        "confidence",
        "audit_notes",
        "article_discrepancy_class",
        "article_final_mngs_by_dtca",
        "article_id_mapping_status",
        "target_hit_details_json",
        "top_virus_taxa_json",
    ]
    for c in detail_cols:
        if c not in audit_df.columns:
            audit_df[c] = pd.NA
    audit_df[detail_cols].to_csv(detail_path, index=False)
    audit_df[
        [
            "File_Name",
            "Run",
            "article_sample_id",
            "specimen_key",
            "record_role",
            "passing_targets",
            "audit_call",
            "sequence_vs_metadata",
            "inconsistency_flags",
            "confidence",
        ]
    ].to_csv(compact_path, index=False)
    pd.DataFrame(errors).to_csv(errors_path, index=False)

    provenance = {
        "generated_at_utc": utcnow(),
        "source_api": API,
        "threshold": THRESHOLD,
        "n_runs_requested": len(runs),
        "n_stat_available": int(audit_df["stat_available"].fillna(False).astype(bool).sum()),
        "n_fetch_errors": len(errors),
        "counts_by_record_role": audit_df["record_role"].value_counts().to_dict(),
        "counts_by_sequence_vs_metadata": audit_df["sequence_vs_metadata"]
        .value_counts(dropna=False)
        .to_dict(),
        "counts_by_audit_call": audit_df["audit_call"].value_counts(dropna=False).to_dict(),
        "limitations": [
            "STAT total_count is a k-mer/spot taxonomy count, not SURPI+ non-overlapping genomic coverage.",
            "Threshold >=3 total_count is used as a proxy for the paper's >=3 non-overlapping reads/contigs rule.",
            "Deposited FASTQs are human-depleted (nohuman); STAT reflects that composition.",
            "Do not invent NP/NUC to P-NP crosswalks; unresolved aliases remain metadata-incomplete.",
        ],
        "outputs": {
            "detail": str(detail_path),
            "summary": str(compact_path),
            "errors": str(errors_path),
            "cache_dir": str(CACHE_DIR),
        },
    }
    provenance_path.write_text(json.dumps(provenance, indent=2))
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()
