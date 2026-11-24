#!/usr/bin/env python3
"""Independent biological-content audit of Accuracy FASTQ / SRA records.

Resource-aware strategy
-----------------------
Full local re-analysis of all 281 human-subtracted FASTQs would require
~11 GB of downloads plus a viral classifier database and tools that are not
installed on this host (no kraken2/centrifuge/minimap2/sra-toolkit; WSL bash
unavailable). NCBI S3 STAT dumps are ~209 GB. Instead this script:

1. Uses SRA runinfo already in the Accuracy manifest for QC/read summaries.
2. Fetches accession-provided NCBI STAT taxonomy via the public Run Browser
   backend (`run_taxonomy`) for every SRR — sequence-derived evidence, not
   clinical metadata inference.
3. Applies a paper-consistent >=3 assigned-spot threshold at clinically
   relevant ranks (proxy for >=3 non-overlapping reads/contigs; see
   limitations).
4. Compares observed respiratory-virus calls to expected clinical/control
   labels, with batch/control context.
5. Optionally stages ENA FASTQ downloads for candidate mismatches when
   `--confirm-candidates` is set (alignment tools still required for true
   coverage confirmation).

This phase does NOT write the final narrative report.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import platform
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "data" / "normalized" / "accuracy_manifest.csv"
DEFAULT_CLINICAL = ROOT / "data" / "normalized" / "clinical_metadata_samples.csv"
DEFAULT_CLINICAL_LINES = ROOT / "data" / "normalized" / "clinical_metadata_lines.csv"
OUT_DIR = ROOT / "results" / "intermediate" / "sequence_audit"
STAT_CACHE = OUT_DIR / "stat_cache"
STAT_URL = (
    "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy"
    "?acc={acc}&cluster_name=public"
)
ENA_FILEREPORT = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    "?accession={acc}&result=read_run"
    "&fields=run_accession,fastq_ftp,submitted_ftp&format=tsv"
)

VIRUSES_TAX_ID = 10239
MIN_EVIDENCE = 3  # paper: >=3 non-overlapping reads/contigs; STAT proxy

# Accuplex Verification Panel (paper Methods): SARS-CoV-2, Flu A, Flu B, RSV.
# MS2 phage is the internal qualitative control spiked into all samples.
PC_EXPECTED_KEYS = {
    "sars_cov_2",
    "influenza_a",
    "influenza_b",
    "rsv",
}
INTERNAL_CONTROL_KEYS = {"ms2"}

# Clinically actionable / panel-relevant respiratory virus families for audit.
CLINICAL_VIRUS_DEFS: list[dict[str, Any]] = [
    {
        "key": "sars_cov_2",
        "label": "SARS-CoV-2",
        "call_patterns": [
            r"sars[- ]?cov[- ]?2",
            r"covid",
        ],
        "org_patterns": [
            r"severe acute respiratory syndrome coronavirus 2",
            r"sars[- ]?cov[- ]?2",
            r"\bbetacoronavirus\b",
        ],
        "tax_ids": {2697049, 694009, 2509511},
    },
    {
        "key": "influenza_a",
        "label": "Influenza A",
        "call_patterns": [r"flu\s*a", r"influenza\s*a", r"h3n2", r"h1n1"],
        "org_patterns": [
            r"influenza a virus",
            r"alphainfluenzavirus",
            r"influenza a\b",
        ],
        "tax_ids": {11320, 197911, 197912},
    },
    {
        "key": "influenza_b",
        "label": "Influenza B",
        "call_patterns": [r"flu\s*b", r"influenza\s*b"],
        "org_patterns": [r"influenza b virus", r"betainfluenzavirus", r"influenza b\b"],
        "tax_ids": {11520, 197913},
    },
    {
        "key": "rsv",
        "label": "RSV",
        "call_patterns": [r"\brsv\b", r"respiratory syncytial"],
        "org_patterns": [
            r"respiratory syncytial virus",
            r"orthopneumovirus",
            r"human orthopneumovirus",
            r"\brsv\b",
        ],
        "tax_ids": {11250, 208893, 208895, 1868221},
    },
    {
        "key": "hmpv",
        "label": "Metapneumovirus",
        "call_patterns": [r"metapneumovirus", r"\bhmpv\b", r"\bmpv\b"],
        "org_patterns": [r"metapneumovirus", r"human metapneumovirus"],
        "tax_ids": {162145, 162387},
    },
    {
        "key": "adenovirus",
        "label": "Adenovirus",
        "call_patterns": [r"adenovirus", r"mastadenovirus"],
        "org_patterns": [r"adenovirus", r"mastadenovirus", r"adenoviridae"],
        "tax_ids": {10508, 10509, 130308},
    },
    {
        "key": "rhinovirus",
        "label": "Rhinovirus",
        "call_patterns": [r"rhinovirus"],
        "org_patterns": [r"rhinovirus", r"enterovirus.*(rhino|b|c)\b"],
        "tax_ids": {147711, 147712, 138951, 138949},
    },
    {
        "key": "enterovirus",
        "label": "Enterovirus",
        "call_patterns": [r"enterovirus", r"\bevd?68\b", r"rhinovirus/enterovirus"],
        "org_patterns": [r"enterovirus", r"ev[- ]?d68", r"enterovirus d"],
        "tax_ids": {12059, 138950},
    },
    {
        "key": "parainfluenza",
        "label": "Parainfluenza",
        "call_patterns": [r"parainfluenza", r"parainFluenza"],
        "org_patterns": [
            r"parainfluenza",
            r"respirovirus",
            r"rubulavirus",
            r"human respirovirus",
        ],
        "tax_ids": {11216, 12730, 11203},
    },
    {
        "key": "ms2",
        "label": "MS2 bacteriophage",
        "call_patterns": [r"\bms2\b"],
        "org_patterns": [
            r"\bms2\b",
            r"emesvirus zinderi",
            r"enterobacteria phage ms2",
            r"escherichia virus ms2",
        ],
        "tax_ids": {329852, 2842664},
    },
]

# Background / non-panel taxa commonly seen; flag but do not treat as
# clinical mismatches by themselves.
BACKGROUND_ORG_PATTERNS = [
    r"anellovir",
    r"torque teno",
    r"ttv",
    r"caudoviricetes",
    r"bacteriophage",
    r"phage",
    r"papillomavirus",
    r"polyomavirus",
    r"herpesvir",
    r"cytomegalovirus",
    r"epstein",
    r"sindbis",
    r"alphavirus",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_versions() -> dict[str, Any]:
    versions: dict[str, Any] = {
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "requests": requests.__version__,
    }
    for tool in (
        "fasterq-dump",
        "prefetch",
        "kraken2",
        "centrifuge",
        "minimap2",
        "bowtie2",
        "blastn",
        "samtools",
        "aws",
        "curl",
    ):
        path = shutil.which(tool)
        versions[f"which_{tool}"] = path
        if path and tool in {"aws"}:
            try:
                out = subprocess.check_output(
                    [path, "--version"], text=True, stderr=subprocess.STDOUT, timeout=15
                )
                versions[f"{tool}_version"] = out.strip().splitlines()[0]
            except Exception as exc:  # noqa: BLE001
                versions[f"{tool}_version_error"] = str(exc)
    return versions


def compile_defs() -> list[dict[str, Any]]:
    compiled = []
    for item in CLINICAL_VIRUS_DEFS:
        compiled.append(
            {
                **item,
                "call_re": [re.compile(p, re.I) for p in item["call_patterns"]],
                "org_re": [re.compile(p, re.I) for p in item["org_patterns"]],
            }
        )
    return compiled


DEFS = compile_defs()
BACKGROUND_RE = [re.compile(p, re.I) for p in BACKGROUND_ORG_PATTERNS]


def parse_expected_keys(text: str | float | None) -> set[str]:
    if text is None or (isinstance(text, float) and pd.isna(text)):
        return set()
    s = str(text).strip()
    if not s or s in {"-", "nan", "None"}:
        return set()
    # Treat "not detected" phrases as absence for that token only when whole-call.
    keys: set[str] = set()
    # Split multi-call strings.
    parts = re.split(r"\|\||\band\b|;", s, flags=re.I)
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if re.search(r"not detected|undetected|negative", part, re.I):
            continue
        for d in DEFS:
            if any(rx.search(part) for rx in d["call_re"]):
                keys.add(d["key"])
    # Rhinovirus/Enterovirus panel call implies both families possible.
    if re.search(r"rhinovirus\s*/\s*enterovirus", s, re.I):
        keys.add("rhinovirus")
        keys.add("enterovirus")
    return keys


def is_background_org(org: str | None) -> bool:
    if not org:
        return False
    return any(rx.search(org) for rx in BACKGROUND_RE)


def match_virus_keys(org: str | None, tax_id: int | None) -> set[str]:
    keys: set[str] = set()
    for d in DEFS:
        if tax_id is not None and tax_id in d["tax_ids"]:
            keys.add(d["key"])
            continue
        if org and any(rx.search(org) for rx in d["org_re"]):
            keys.add(d["key"])
    return keys


def virus_subtree(tax_table: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_id: dict[int, dict[str, Any]] = {}
    for row in tax_table:
        try:
            tid = int(row["tax_id"])
        except (KeyError, TypeError, ValueError):
            continue
        by_id[tid] = row
    if VIRUSES_TAX_ID not in by_id:
        # Fallback: rows whose org mentions virus-like names already covered.
        return [r for r in tax_table if match_virus_keys(r.get("org"), None)]

    virus_ids = {VIRUSES_TAX_ID}
    changed = True
    while changed:
        changed = False
        for tid, row in by_id.items():
            parent = row.get("parent")
            if parent is None:
                continue
            try:
                parent_i = int(parent)
            except (TypeError, ValueError):
                continue
            if parent_i in virus_ids and tid not in virus_ids:
                virus_ids.add(tid)
                changed = True
    out = []
    for tid in virus_ids:
        row = dict(by_id[tid])
        row["tax_id"] = tid
        out.append(row)
    return out


def supported_clinical_calls(
    tax_table: list[dict[str, Any]], min_evidence: int = MIN_EVIDENCE
) -> tuple[set[str], list[dict[str, Any]], list[dict[str, Any]]]:
    """Return supported clinical keys, evidence rows, and low-confidence flags."""
    viruses = virus_subtree(tax_table)
    evidence_rows: list[dict[str, Any]] = []
    low_conf: list[dict[str, Any]] = []
    supported: set[str] = set()

    for row in viruses:
        org = row.get("org")
        try:
            tid = int(row.get("tax_id"))
        except (TypeError, ValueError):
            tid = None
        self_count = int(row.get("self_count") or 0)
        total_count = int(row.get("total_count") or 0)
        rank = (row.get("rank") or "").lower()
        keys = match_virus_keys(org, tid)
        rec = {
            "tax_id": tid,
            "org": org,
            "rank": rank,
            "self_count": self_count,
            "total_count": total_count,
            "matched_keys": sorted(keys),
            "background": is_background_org(org),
        }
        # Prefer self_count at species/no-rank leaves; else total_count at genus+.
        evidence = self_count if self_count >= min_evidence else (
            total_count if rank in {"species", "genus", "no rank", ""} else 0
        )
        # For higher ranks mapped to clinical keys, allow total_count.
        if keys and evidence < min_evidence and total_count >= min_evidence:
            evidence = total_count
        rec["evidence_count"] = evidence
        if keys and evidence >= min_evidence:
            supported |= keys
            evidence_rows.append(rec)
        elif keys and (self_count > 0 or total_count > 0):
            low_conf.append(rec)
        elif (
            not rec["background"]
            and rank in {"species", "genus"}
            and total_count >= min_evidence
            and org
            and "virus" in org.lower()
        ):
            # Unexpected non-panel virus with evidence.
            low_conf.append({**rec, "unexpected_nonpanel": True})
    return supported, evidence_rows, low_conf


def fetch_stat(acc: str, session: requests.Session, cache_dir: Path, retries: int = 4) -> dict[str, Any]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{acc}.json"
    if cache_path.exists() and cache_path.stat().st_size > 10:
        return json.loads(cache_path.read_text(encoding="utf-8"))

    url = STAT_URL.format(acc=acc)
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            resp = session.get(url, timeout=180)
            if resp.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            resp.raise_for_status()
            data = resp.json()
            obj = data[0] if isinstance(data, list) else data
            cache_path.write_text(json.dumps(obj), encoding="utf-8")
            return obj
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"STAT fetch failed for {acc}: {last_err}")


def load_clinical_expectations(
    clinical_samples: Path, clinical_lines: Path
) -> dict[str, dict[str, Any]]:
    samples = pd.read_csv(clinical_samples)
    lines = pd.read_csv(clinical_lines)
    out: dict[str, dict[str, Any]] = {}
    for _, row in samples.iterrows():
        sid = row["article_sample_id"]
        pcr_keys = parse_expected_keys(row.get("pcr_results"))
        mngs_keys = parse_expected_keys(row.get("mngs_results"))
        # Prefer union of PCR + reported mNGS as "metadata-declared content".
        # For consistency vs clinical label, PCR is primary expected clinical call;
        # mNGS is published sequence call (not independent).
        out[sid] = {
            "pcr_results": row.get("pcr_results"),
            "mngs_results": row.get("mngs_results"),
            "pcr_keys": sorted(pcr_keys),
            "mngs_keys": sorted(mngs_keys),
            "expected_clinical_keys": sorted(pcr_keys),
            "published_mngs_keys": sorted(mngs_keys),
            "comments": row.get("comments"),
        }
    # Enrich from lines if sample-level aggregation missed tokens.
    for sid, grp in lines.groupby("article_sample_id"):
        pcr_keys: set[str] = set()
        mngs_keys: set[str] = set()
        for _, row in grp.iterrows():
            pcr_keys |= parse_expected_keys(row.get("pcr_results"))
            mngs_keys |= parse_expected_keys(row.get("mngs_results"))
        if sid not in out:
            out[sid] = {
                "pcr_results": None,
                "mngs_results": None,
                "pcr_keys": sorted(pcr_keys),
                "mngs_keys": sorted(mngs_keys),
                "expected_clinical_keys": sorted(pcr_keys),
                "published_mngs_keys": sorted(mngs_keys),
                "comments": None,
            }
        else:
            out[sid]["pcr_keys"] = sorted(set(out[sid]["pcr_keys"]) | pcr_keys)
            out[sid]["mngs_keys"] = sorted(set(out[sid]["mngs_keys"]) | mngs_keys)
            out[sid]["expected_clinical_keys"] = sorted(
                set(out[sid]["expected_clinical_keys"]) | pcr_keys
            )
            out[sid]["published_mngs_keys"] = sorted(
                set(out[sid]["published_mngs_keys"]) | mngs_keys
            )
    return out


def classify_record(
    row: pd.Series,
    clinical: dict[str, dict[str, Any]],
    supported: set[str],
    evidence_rows: list[dict[str, Any]],
    low_conf: list[dict[str, Any]],
    tax_totals: dict[str, Any] | None,
    batch_pc_keys: set[str],
) -> dict[str, Any]:
    role = row["record_role"]
    article_id = row.get("article_sample_id")
    if isinstance(article_id, float) and pd.isna(article_id):
        article_id = None
    clin = clinical.get(str(article_id), {}) if article_id else {}

    expected: set[str] = set()
    expectation_source = "none"
    if role == "positive_control":
        expected = set(PC_EXPECTED_KEYS)
        expectation_source = "paper_accuplex_PC"
    elif role == "negative_control":
        expected = set()
        expectation_source = "paper_pooled_negative_NC"
    elif clin.get("expected_clinical_keys") is not None and article_id:
        expected = set(clin.get("expected_clinical_keys") or [])
        expectation_source = "supplementary_PCR"
    elif row.get("clinical_metadata_join") in {
        "unresolved_NP_NUC_alias",
        "not_in_supplement",
    }:
        expectation_source = f"metadata_unavailable:{row.get('clinical_metadata_join')}"
    else:
        expectation_source = "metadata_unavailable"

    # Observed clinical-relevant keys exclude MS2 for clinical mismatch logic,
    # but keep MS2 separately for internal-control checks.
    observed_clinical = {k for k in supported if k not in INTERNAL_CONTROL_KEYS}
    observed_ms2 = "ms2" in supported

    missing = sorted(expected - observed_clinical) if expectation_source.startswith(
        ("paper_", "supplementary_")
    ) else []
    unexpected = (
        sorted(observed_clinical - expected)
        if expectation_source.startswith(("paper_", "supplementary_"))
        else []
    )

    # Soften unexpected when same taxon is in batch PC (possible contamination).
    contamination_context = []
    for key in unexpected:
        if key in batch_pc_keys:
            contamination_context.append(key)

    status = "not_assessable"
    confidence = "low"
    explanation_parts: list[str] = []

    if tax_totals is None:
        status = "stat_fetch_failed"
        confidence = "none"
        explanation_parts.append("NCBI STAT taxonomy could not be retrieved.")
    elif expectation_source.startswith("metadata_unavailable"):
        status = "sequence_profiled_metadata_unlinked"
        confidence = "medium"
        explanation_parts.append(
            "STAT viral profile obtained, but no linked PCR expectation for this record."
        )
        if observed_clinical:
            explanation_parts.append(
                "Supported clinical taxa: " + ", ".join(sorted(observed_clinical))
            )
    elif role == "positive_control":
        if not missing and not unexpected:
            status = "consistent_control"
            confidence = "high"
            explanation_parts.append(
                "PC Accuplex targets detected at >=3 STAT spots; no unexpected panel viruses."
            )
        elif not missing and unexpected:
            status = "control_extra_panel_virus"
            confidence = "medium"
            explanation_parts.append(
                "All Accuplex targets present; additional panel viruses observed: "
                + ", ".join(unexpected)
            )
        elif missing:
            status = "candidate_inconsistency"
            confidence = "medium"
            explanation_parts.append(
                "PC missing expected Accuplex targets: " + ", ".join(missing)
            )
        if not observed_ms2:
            explanation_parts.append("MS2 internal-control signal not supported at threshold.")
    elif role == "negative_control":
        if not observed_clinical:
            status = "consistent_control"
            confidence = "high"
            explanation_parts.append(
                "No clinical panel virus supported at >=3 STAT spots."
            )
        else:
            status = "candidate_inconsistency"
            confidence = "medium"
            explanation_parts.append(
                "NC has clinical panel virus signal: " + ", ".join(sorted(observed_clinical))
            )
            if contamination_context:
                explanation_parts.append(
                    "Same taxon present in batch PC: " + ", ".join(contamination_context)
                )
    else:
        # Clinical sample with PCR expectation.
        if not missing and not unexpected:
            status = "consistent_with_PCR"
            confidence = "high"
            explanation_parts.append("Supported STAT taxa match PCR-expected panel viruses.")
        elif missing and not unexpected:
            status = "candidate_inconsistency"
            confidence = "medium"
            explanation_parts.append(
                "PCR-expected virus not supported by STAT >=3 spots: "
                + ", ".join(missing)
            )
            # If published mNGS also missed it, may be true low titer / FN.
            pub = set(clin.get("published_mngs_keys") or [])
            if missing and set(missing).isdisjoint(pub):
                explanation_parts.append(
                    "Published mNGS also lacked these taxa (possible mNGS FN / low titer)."
                )
        elif unexpected and not missing:
            status = "candidate_inconsistency"
            confidence = "medium"
            explanation_parts.append(
                "STAT supports panel viruses not in PCR call: " + ", ".join(unexpected)
            )
            pub = set(clin.get("published_mngs_keys") or [])
            if set(unexpected) <= pub:
                explanation_parts.append(
                    "Matches published mNGS extras (often RVP-missed or non-targeted)."
                )
                status = "sequence_supports_published_mNGS_extra"
                confidence = "high"
            if contamination_context:
                explanation_parts.append(
                    "Batch PC also contains: " + ", ".join(contamination_context)
                )
        else:
            status = "candidate_inconsistency"
            confidence = "medium"
            explanation_parts.append(
                f"Missing={missing}; unexpected={unexpected}"
            )

    # QC summary from SRA metadata (not local FastQC).
    spots = int(row.get("spots") or 0)
    qc_flag = "pass_sra_spots" if spots >= 5_000_000 else "below_paper_5M_preprocessed_spots_proxy"
    # Paper QC uses preprocessed reads >=5M; deposited nohuman spots are post human subtraction,
    # so this is a conservative proxy only.

    return {
        "Run": row["Run"],
        "File_Name": row["File_Name"],
        "Sample_Name": row["Sample_Name"],
        "Original_ID": row["Original_ID"],
        "specimen_key": row.get("specimen_key"),
        "article_sample_id": article_id,
        "record_role": role,
        "is_control": bool(row.get("is_control")),
        "is_rerun_or_replicate": bool(row.get("is_rerun_or_replicate")),
        "mNGS_run": row.get("mNGS_run"),
        "clinical_metadata_join": row.get("clinical_metadata_join"),
        "spots": spots,
        "bases": int(row.get("bases") or 0),
        "avgLength": row.get("avgLength"),
        "size_MB": row.get("size_MB"),
        "qc_flag_sra_metadata": qc_flag,
        "stat_analysed_spots": (tax_totals or {}).get("analysed"),
        "stat_identified_spots": (tax_totals or {}).get("identified"),
        "stat_total_spots": (tax_totals or {}).get("total"),
        "stat_version": (tax_totals or {}).get("version"),
        "stat_last_update": (tax_totals or {}).get("last_update"),
        "expectation_source": expectation_source,
        "expected_keys": sorted(expected),
        "observed_clinical_keys": sorted(observed_clinical),
        "observed_ms2": observed_ms2,
        "missing_keys": missing,
        "unexpected_keys": unexpected,
        "batch_pc_overlap_unexpected": contamination_context,
        "pcr_results": clin.get("pcr_results"),
        "mngs_results": clin.get("mngs_results"),
        "published_mngs_keys": clin.get("published_mngs_keys"),
        "evidence_json": json.dumps(evidence_rows, sort_keys=True),
        "low_confidence_json": json.dumps(low_conf, sort_keys=True),
        "consistency_status": status,
        "confidence": confidence,
        "evidence_type": "ncbi_stat_accession_taxonomy",
        "explanation": " ".join(explanation_parts),
    }


def build_batch_pc_index(
    manifest: pd.DataFrame, per_run_supported: dict[str, set[str]]
) -> dict[str, set[str]]:
    """Map mNGS_run -> clinical keys supported in that run's positive controls."""
    out: dict[str, set[str]] = defaultdict(set)
    pcs = manifest[manifest["record_role"] == "positive_control"]
    for _, row in pcs.iterrows():
        keys = per_run_supported.get(row["Run"], set())
        out[str(row["mNGS_run"])] |= {k for k in keys if k not in INTERNAL_CONTROL_KEYS}
    return out


def download_ena_fastq(acc: str, dest_dir: Path, session: requests.Session) -> list[Path]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    url = ENA_FILEREPORT.format(acc=acc)
    resp = session.get(url, timeout=120)
    resp.raise_for_status()
    lines = resp.text.strip().splitlines()
    if len(lines) < 2:
        return []
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    ftp = row.get("fastq_ftp") or row.get("submitted_ftp") or ""
    paths: list[Path] = []
    for part in ftp.split(";"):
        part = part.strip()
        if not part:
            continue
        if part.startswith("ftp://"):
            http = "https://" + part[len("ftp://") :]
        elif part.startswith("ftp."):
            http = "https://" + part
        else:
            http = "https://" + part
        name = Path(part).name
        out = dest_dir / name
        if out.exists() and out.stat().st_size > 0:
            paths.append(out)
            continue
        with session.get(http, stream=True, timeout=600) as r:
            r.raise_for_status()
            tmp = out.with_suffix(out.suffix + ".part")
            with tmp.open("wb") as handle:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
            tmp.replace(out)
        paths.append(out)
    return paths


def summarize_fastq_reads(path: Path, max_reads: int = 20000) -> dict[str, Any]:
    """Lightweight QC on a staged FASTQ.gz (read count sample + lengths)."""
    n = 0
    lengths: list[int] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        while n < max_reads:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            handle.readline()
            handle.readline()
            if not seq:
                break
            lengths.append(len(seq.strip()))
            n += 1
    return {
        "file": str(path),
        "sampled_reads": n,
        "mean_length": (sum(lengths) / n) if n else None,
        "min_length": min(lengths) if lengths else None,
        "max_length": max(lengths) if lengths else None,
        "note": f"sampled first {max_reads} reads only; not full-file count",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--clinical-samples", type=Path, default=DEFAULT_CLINICAL)
    parser.add_argument("--clinical-lines", type=Path, default=DEFAULT_CLINICAL_LINES)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--limit", type=int, default=0, help="Optional debug limit")
    parser.add_argument(
        "--confirm-candidates",
        action="store_true",
        help="Stage ENA FASTQs for candidate inconsistencies (QC sample only)",
    )
    parser.add_argument(
        "--max-confirm",
        type=int,
        default=10,
        help="Max candidate FASTQs to stage when confirming",
    )
    args = parser.parse_args()

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "stat_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    resource_note = {
        "assessed_at_utc": utc_now(),
        "accuracy_records": None,
        "manifest_size_MB_sum": None,
        "manifest_size_GiB": None,
        "disk_free_bytes_approx": shutil.disk_usage(ROOT).free,
        "ncbi_stat_s3_tax_analysis_GiB": 208.8,
        "ncbi_stat_s3_tax_analysis_info_GiB": 2.1,
        "local_classifier_tools": "missing",
        "chosen_strategy": (
            "accession-provided NCBI STAT taxonomy for all 281 runs; "
            "SRA runinfo QC; optional staged ENA FASTQ for candidates"
        ),
        "full_fastq_download_deferred": True,
        "reason": (
            "~11 GB FASTQs are within free disk, but no local viral classifier/"
            "aligner is installed and WSL bash is unavailable; NCBI STAT provides "
            "sequence-derived taxonomy without transferring 209 GB STAT dumps or "
            "building a DB."
        ),
    }

    manifest = pd.read_csv(args.manifest)
    assert len(manifest) == 281, f"Expected 281 Accuracy rows, found {len(manifest)}"
    if args.limit:
        manifest = manifest.head(args.limit).copy()

    resource_note["accuracy_records"] = int(len(pd.read_csv(args.manifest)))
    size_sum = float(pd.read_csv(args.manifest)["size_MB"].fillna(0).sum())
    resource_note["manifest_size_MB_sum"] = size_sum
    resource_note["manifest_size_GiB"] = round(size_sum / 1024, 2)

    clinical = load_clinical_expectations(args.clinical_samples, args.clinical_lines)
    versions = tool_versions()

    session = requests.Session()
    session.headers.update(
        {"User-Agent": "c121-metagenomics-accuracy-audit/1.0 (research; local)"}
    )

    runs = manifest["Run"].astype(str).tolist()
    fetch_errors: dict[str, str] = {}
    stat_objects: dict[str, dict[str, Any]] = {}

    def _job(acc: str) -> tuple[str, dict[str, Any] | None, str | None]:
        local = requests.Session()
        local.headers.update(session.headers)
        try:
            return acc, fetch_stat(acc, local, cache_dir), None
        except Exception as exc:  # noqa: BLE001
            return acc, None, str(exc)

    print(f"Fetching STAT taxonomy for {len(runs)} runs with {args.workers} workers...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futs = [pool.submit(_job, acc) for acc in runs]
        for i, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            acc, obj, err = fut.result()
            if err:
                fetch_errors[acc] = err
            else:
                assert obj is not None
                stat_objects[acc] = obj
            if i % 25 == 0 or i == len(runs):
                print(f"  STAT progress {i}/{len(runs)} (errors={len(fetch_errors)})")

    per_run_supported: dict[str, set[str]] = {}
    per_run_evidence: dict[str, list[dict[str, Any]]] = {}
    per_run_low: dict[str, list[dict[str, Any]]] = {}
    virus_rows_out: list[dict[str, Any]] = []

    for acc, obj in stat_objects.items():
        supported, evidence, low = supported_clinical_calls(obj.get("tax_table") or [])
        per_run_supported[acc] = supported
        per_run_evidence[acc] = evidence
        per_run_low[acc] = low
        for row in virus_subtree(obj.get("tax_table") or []):
            self_count = int(row.get("self_count") or 0)
            total_count = int(row.get("total_count") or 0)
            if self_count < 1 and total_count < MIN_EVIDENCE:
                continue
            try:
                tid = int(row.get("tax_id"))
            except (TypeError, ValueError):
                tid = None
            virus_rows_out.append(
                {
                    "Run": acc,
                    "tax_id": tid,
                    "org": row.get("org"),
                    "rank": row.get("rank"),
                    "self_count": self_count,
                    "total_count": total_count,
                    "matched_keys": "|".join(sorted(match_virus_keys(row.get("org"), tid))),
                    "background": is_background_org(row.get("org")),
                }
            )

    batch_pc = build_batch_pc_index(manifest, per_run_supported)

    findings: list[dict[str, Any]] = []
    for _, row in manifest.iterrows():
        acc = str(row["Run"])
        obj = stat_objects.get(acc)
        supported = per_run_supported.get(acc, set())
        evidence = per_run_evidence.get(acc, [])
        low = per_run_low.get(acc, [])
        tax_totals = (obj or {}).get("tax_totals")
        finding = classify_record(
            row,
            clinical,
            supported,
            evidence,
            low,
            tax_totals,
            batch_pc.get(str(row.get("mNGS_run")), set()),
        )
        if acc in fetch_errors:
            finding["stat_error"] = fetch_errors[acc]
        findings.append(finding)

    findings_df = pd.DataFrame(findings)
    findings_path = out_dir / "accuracy_sequence_audit_findings.csv"
    findings_df.to_csv(findings_path, index=False)

    virus_df = pd.DataFrame(virus_rows_out)
    virus_path = out_dir / "accuracy_stat_virus_taxa.csv"
    virus_df.to_csv(virus_path, index=False)

    candidates = findings_df[
        findings_df["consistency_status"].isin(
            [
                "candidate_inconsistency",
                "control_extra_panel_virus",
            ]
        )
    ].copy()
    candidates_path = out_dir / "accuracy_sequence_candidate_mismatches.csv"
    candidates.to_csv(candidates_path, index=False)

    confirm_results: list[dict[str, Any]] = []
    if args.confirm_candidates and len(candidates):
        fastq_dir = out_dir / "staged_fastq"
        subset = candidates.head(args.max_confirm)
        print(f"Staging up to {len(subset)} candidate FASTQs from ENA...")
        for _, crow in subset.iterrows():
            acc = crow["Run"]
            try:
                paths = download_ena_fastq(acc, fastq_dir / acc, session)
                for p in paths:
                    qc = summarize_fastq_reads(p)
                    confirm_results.append(
                        {
                            "Run": acc,
                            "consistency_status": crow["consistency_status"],
                            "fastq_path": str(p),
                            "bytes": p.stat().st_size,
                            "sha256": sha256_file(p),
                            **qc,
                            "alignment_confirmation": "blocked_no_local_aligner",
                            "note": (
                                "FASTQ staged for candidate; targeted genome alignment/"
                                "coverage not executed (minimap2/blastn absent)."
                            ),
                        }
                    )
            except Exception as exc:  # noqa: BLE001
                confirm_results.append(
                    {
                        "Run": acc,
                        "consistency_status": crow["consistency_status"],
                        "error": str(exc),
                        "alignment_confirmation": "download_failed",
                    }
                )
        pd.DataFrame(confirm_results).to_csv(
            out_dir / "accuracy_candidate_fastq_staging.csv", index=False
        )

    # Coverage / status tallies
    status_counts = findings_df["consistency_status"].value_counts().to_dict()
    coverage = {
        "accuracy_records_targeted": int(len(manifest)),
        "stat_taxonomy_retrieved": int(len(stat_objects)),
        "stat_fetch_failures": int(len(fetch_errors)),
        "coverage_fraction": round(len(stat_objects) / max(len(manifest), 1), 4),
        "candidate_inconsistencies": int(
            (findings_df["consistency_status"] == "candidate_inconsistency").sum()
        ),
        "sequence_supports_published_mNGS_extra": int(
            (
                findings_df["consistency_status"]
                == "sequence_supports_published_mNGS_extra"
            ).sum()
        ),
        "consistent_with_PCR": int(
            (findings_df["consistency_status"] == "consistent_with_PCR").sum()
        ),
        "consistent_control": int(
            (findings_df["consistency_status"] == "consistent_control").sum()
        ),
        "sequence_profiled_metadata_unlinked": int(
            (
                findings_df["consistency_status"]
                == "sequence_profiled_metadata_unlinked"
            ).sum()
        ),
        "status_counts": status_counts,
    }

    commands = {
        "resource_assessment": [
            "python -c \"import pandas as pd; d=pd.read_csv('data/normalized/accuracy_manifest.csv'); print(len(d), d.size_MB.sum())\"",
            "aws s3 ls s3://sra-pub-metadata-us-east-1/sra_tax_analysis_tool/tax_analysis/ --no-sign-request --summarize --human-readable",
            "Get-Command fasterq-dump,kraken2,minimap2,blastn -ErrorAction SilentlyContinue",
        ],
        "stat_endpoint": STAT_URL,
        "audit_command": "python analysis/audit_accuracy_fastq.py --confirm-candidates",
    }

    limitations = [
        "Evidence is NCBI STAT k-mer spot assignments on deposited human-subtracted SRA runs, not a local SURPI+/blast reanalysis of FASTQs.",
        "The paper threshold is >=3 non-overlapping reads/contigs to a target genome; STAT self_count/total_count >=3 is a reproducible proxy but does not prove non-overlapping genomic coverage.",
        "No local aligner/classifier was available; candidate targeted alignment/coverage confirmation is blocked unless tools are installed.",
        "Full ENA/SRA FASTQ download (~11 GB) was deferred after resource assessment; optional staging covers candidates only when --confirm-candidates is used.",
        "SRA spot counts are post-deposit metadata and are not identical to the paper's preprocessed-read QC metric (minimum 5M).",
        "79 NP/NUC Accuracy records lack Supplementary PCR crosswalks; they are profiled but not judged against clinical PCR expectations.",
        "Positive-control expectation uses Accuplex panel composition from the paper Methods (SARS-CoV-2, Flu A, Flu B, RSV) plus MS2 internal control.",
        "Unexpected non-panel viruses (anellovirus, phage, herpesvirus, papillomavirus, Sindbis) are recorded in taxon tables but not auto-labeled as clinical mislabels.",
    ]

    blockers = []
    if fetch_errors:
        blockers.append(
            {
                "blocker": "stat_fetch_failures",
                "n": len(fetch_errors),
                "examples": dict(list(fetch_errors.items())[:10]),
            }
        )
    blockers.append(
        {
            "blocker": "no_local_aligner_for_targeted_coverage_confirmation",
            "detail": (
                "minimap2/blastn/bowtie2/kraken2/sra-toolkit missing; WSL bash unavailable. "
                "Candidate mismatches retain STAT evidence only unless FASTQs are staged "
                "and an aligner is provided."
            ),
        }
    )
    if not args.confirm_candidates:
        blockers.append(
            {
                "blocker": "full_fastq_corpus_not_downloaded",
                "detail": (
                    "By design after size assessment (~11 GB FASTQ; ~209 GB STAT S3). "
                    "Re-run with --confirm-candidates to stage candidate FASTQs."
                ),
            }
        )

    provenance = {
        "generated_at_utc": utc_now(),
        "scope": "audit-sequences only; no final narrative report",
        "inputs": {
            "manifest": str(args.manifest),
            "clinical_samples": str(args.clinical_samples),
            "clinical_lines": str(args.clinical_lines),
        },
        "resource_assessment": resource_note,
        "tool_versions": versions,
        "commands": commands,
        "coverage": coverage,
        "limitations": limitations,
        "blockers": blockers,
        "outputs": {
            "findings": str(findings_path),
            "virus_taxa": str(virus_path),
            "candidates": str(candidates_path),
            "stat_cache_dir": str(cache_dir),
        },
        "min_evidence_threshold": MIN_EVIDENCE,
        "evidence_rule_note": (
            "Supported call requires matched clinical virus key with STAT "
            "self_count or qualifying total_count >= 3."
        ),
    }
    prov_path = out_dir / "accuracy_sequence_audit_provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    summary_md = out_dir / "README_sequence_audit.md"
    summary_md.write_text(
        "\n".join(
            [
                "# Accuracy sequence audit (intermediate)",
                "",
                f"Generated: `{provenance['generated_at_utc']}`",
                "",
                "## Coverage",
                f"- STAT profiles: **{coverage['stat_taxonomy_retrieved']}/{coverage['accuracy_records_targeted']}**",
                f"- Candidate inconsistencies: **{coverage['candidate_inconsistencies']}**",
                f"- Consistent with PCR: **{coverage['consistent_with_PCR']}**",
                f"- Consistent controls: **{coverage['consistent_control']}**",
                f"- Profiled but PCR-unlinked: **{coverage['sequence_profiled_metadata_unlinked']}**",
                "",
                "## Evidence type",
                "NCBI STAT accession taxonomy (sequence-derived). Not fabricated. "
                "Not equivalent to local SURPI+ non-overlapping contig evidence.",
                "",
                "## Key outputs",
                f"- `{findings_path.name}`",
                f"- `{virus_path.name}`",
                f"- `{candidates_path.name}`",
                f"- `{prov_path.name}`",
                "",
                "## Limitations",
                *[f"- {x}" for x in limitations],
                "",
                "## Blockers",
                *[f"- {b['blocker']}: {b.get('detail', b)}" for b in blockers],
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(json.dumps(coverage, indent=2))
    print(f"Wrote {findings_path}")
    print(f"Wrote {candidates_path}")
    print(f"Wrote {prov_path}")
    return 0 if not fetch_errors else 2


if __name__ == "__main__":
    raise SystemExit(main())
