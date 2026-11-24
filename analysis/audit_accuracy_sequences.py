#!/usr/bin/env python3
"""Audit Accuracy-cohort sequence content via NCBI STAT taxonomy.

Fetches public STAT run_taxonomy for every Purpose=Accuracy SRR, applies a
>=3-count evidence proxy for non-overlapping viral reads, and compares sequence
evidence against declared clinical/control metadata and article discrepancies.

Does not invent NP/NUC <-> P-NP crosswalks. Sequence evidence and metadata-only
inferences are kept distinct in the outputs.
"""

from __future__ import annotations

import argparse
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
DISC_SAMPLE = ROOT / "data/normalized/article_discrepancies_sample_level.csv"
OUT_DIR = ROOT / "results/intermediate"
CACHE_DIR = OUT_DIR / "stat_cache"
FINDINGS_CSV = OUT_DIR / "accuracy_sequence_audit.csv"
FINDINGS_JSON = OUT_DIR / "accuracy_sequence_audit.json"
VIRUS_HITS_CSV = OUT_DIR / "accuracy_stat_virus_hits.csv"
SUMMARY_JSON = OUT_DIR / "accuracy_sequence_audit_summary.json"

USER_AGENT = "metagenomics-audit/1.0 (+https://github.com/kierro1209/c121_metagenomics)"
MIN_EVIDENCE = 3
VIRUS_ROOT = 10239

# Background / non-respiratory taxa commonly present and not treated as
# respiratory-panel pathogen evidence unless expected by metadata.
BACKGROUND_NAME_RE = re.compile(
    r"(anello|torque\s*teno|\bttv\b|phage|caudovir|fiersvir|emesvirus|levivir|"
    r"lenarvir|unclassified bacterial viruses|gokushovir|\bmicrovir|"
    r"circovir|smacovir|genomovir|\binovir|siphovir|myovir|podovir)",
    re.I,
)

# Map normalized clinical tokens -> STAT organism name match patterns.
ORGANISM_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("rsv", re.compile(r"(orthopneumovirus|respiratory syncytial|human orthopneumovirus|rsv)", re.I)),
    ("influenza_a", re.compile(r"(influenza a|alphainfluenzavirus)", re.I)),
    ("influenza_b", re.compile(r"(influenza b|betainfluenzavirus)", re.I)),
    (
        "sars_cov_2",
        re.compile(
            r"(severe acute respiratory syndrome coronavirus 2|"
            r"severe acute respiratory syndrome-related coronavirus|"
            r"sars.?cov.?2|sars coronavirus 2)",
            re.I,
        ),
    ),
    ("coronavirus_oc43", re.compile(r"(oc43|betacoronavirus 1)", re.I)),
    ("coronavirus_hku1", re.compile(r"(hku1)", re.I)),
    ("coronavirus_nl63", re.compile(r"(nl63)", re.I)),
    ("coronavirus_229e", re.compile(r"(229e)", re.I)),
    # Generic coronavirus family/genus support used as secondary evidence for SARS expectations.
    ("coronaviridae", re.compile(r"^(coronaviridae|orthocoronavirinae|betacoronavirus)$", re.I)),
    ("rhinovirus", re.compile(r"(rhinovirus)", re.I)),
    ("enterovirus", re.compile(r"(enterovirus|ev.?d.?68|coxsackie|echovirus)", re.I)),
    ("parainfluenza", re.compile(r"(parainfluenza|respirovirus|orthorubulavirus|rubulavirus|human respirovirus)", re.I)),
    ("metapneumovirus", re.compile(r"(metapneumovirus|hmpv)", re.I)),
    ("adenovirus", re.compile(r"(adenovirus|mastadenovirus)", re.I)),
    ("bocavirus", re.compile(r"(bocavirus|bocaparvovirus)", re.I)),
    ("cmv", re.compile(r"(cytomegalovirus|human betaherpesvirus 5)", re.I)),
    ("herpesvirus", re.compile(r"(herpesvirus|herpesviridae|roseolovirus|betaherpesvirus|simplexvirus)", re.I)),
    ("papillomavirus", re.compile(r"(papillomavirus|papillomaviridae)", re.I)),
]

# Tokens extracted from pcr/mngs free text that are RVP-panel relevant.
CLINICAL_TOKEN_RE = re.compile(
    r"(rsv\s*[ab]?|flu(?:enza)?\s*[ab]?(?:\s*\(?h3n?2?\)?)?|sars[-\s]?cov[-\s]?2|"
    r"covid|metapneumovirus|mpv|hmpv|adenovirus|mastadenovirus|"
    r"rhinovirus(?:/enterovirus)?|enterovirus|evd?68|parain(?:f|fl)uenza\s*[1-4]?|"
    r"piv\s*[1-4]?|bocavirus|oc43|hku1|nl63|229e|"
    r"cmv|cytomegalovirus|papillomavirus|betaherpesvirus|"
    r"respiratory virus not detected|not detected|nd\b)",
    re.I,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def session() -> requests.Session:
    s = requests.Session()
    s.headers.update({"User-Agent": USER_AGENT})
    return s


def fetch_taxonomy(sess: requests.Session, acc: str, retries: int = 4) -> dict[str, Any]:
    url = (
        "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy"
        f"?acc={acc}&cluster_name=public"
    )
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            r = sess.get(url, timeout=180)
            if r.status_code == 429:
                time.sleep(2 ** attempt)
                continue
            r.raise_for_status()
            data = r.json()
            return data[0] if isinstance(data, list) else data
        except Exception as exc:  # noqa: BLE001 - collect and retry network errors
            last_err = exc
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"STAT fetch failed for {acc}: {last_err}")


def virus_rows(obj: dict[str, Any]) -> list[dict[str, Any]]:
    table = obj.get("tax_table") or []
    by_id = {int(x["tax_id"]): x for x in table if "tax_id" in x}
    if VIRUS_ROOT not in by_id:
        return []
    virus_ids = {VIRUS_ROOT}
    changed = True
    while changed:
        changed = False
        for tid, row in by_id.items():
            parent = row.get("parent")
            if parent is None:
                continue
            parent = int(parent)
            if parent in virus_ids and tid not in virus_ids:
                virus_ids.add(tid)
                changed = True
    rows: list[dict[str, Any]] = []
    for tid in virus_ids:
        row = by_id[tid]
        rows.append(
            {
                "tax_id": tid,
                "org": row.get("org") or "",
                "rank": row.get("rank"),
                "self_count": int(row.get("self_count") or 0),
                "total_count": int(row.get("total_count") or 0),
            }
        )
    rows.sort(key=lambda x: (-x["total_count"], x["org"]))
    return rows


def is_background(org: str) -> bool:
    return bool(BACKGROUND_NAME_RE.search(org or ""))


def evidence_ok(row: dict[str, Any]) -> bool:
    return max(row["self_count"], row["total_count"]) >= MIN_EVIDENCE


def tokenize_clinical(text: str | float | None) -> list[str]:
    if text is None or (isinstance(text, float) and pd.isna(text)):
        return []
    s = str(text).strip()
    if not s or s.lower() in {"-", "nan", "none"}:
        return []
    # Split co-detections; keep only positive-looking clauses.
    parts = re.split(r"\|\||;|,| and ", s)
    tokens: list[str] = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        low = part.lower()
        if re.search(r"not detected|\bnd\b", low) and not re.search(r"\+", low):
            continue
        for m in CLINICAL_TOKEN_RE.finditer(part):
            tok = m.group(0).lower()
            tok = re.sub(r"\s+", " ", tok).strip()
            if tok in {"not detected", "nd", "respiratory virus not detected"}:
                continue
            # Normalize aliases.
            if tok.startswith("rsv"):
                tokens.append("rsv")
            elif tok.startswith("flu") or "influenza" in tok:
                tokens.append("influenza_b" if re.search(r"\bb\b", tok) else "influenza_a")
            elif "sars" in tok or "covid" in tok:
                tokens.append("sars_cov_2")
            elif tok in {"mpv", "hmpv"} or "metapneumovirus" in tok:
                tokens.append("metapneumovirus")
            elif "adeno" in tok or "mastadenovirus" in tok:
                tokens.append("adenovirus")
            elif "rhinovirus" in tok:
                tokens.append("rhinovirus")
            elif "enterovirus" in tok or "evd68" in tok or "ev-d68" in tok or "ev d68" in tok:
                tokens.append("enterovirus")
            elif "parain" in tok or tok.startswith("piv"):
                tokens.append("parainfluenza")
            elif "bocavirus" in tok:
                tokens.append("bocavirus")
            elif "oc43" in tok:
                tokens.append("coronavirus_oc43")
            elif "hku1" in tok:
                tokens.append("coronavirus_hku1")
            elif "nl63" in tok:
                tokens.append("coronavirus_nl63")
            elif "229e" in tok:
                tokens.append("coronavirus_229e")
            elif "cmv" in tok or "cytomegalovirus" in tok:
                tokens.append("cmv")
            elif "papillomavirus" in tok:
                tokens.append("papillomavirus")
            elif "betaherpesvirus" in tok or "herpes" in tok:
                tokens.append("herpesvirus")
    # Preserve order, unique.
    seen: set[str] = set()
    out: list[str] = []
    for t in tokens:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def matched_organisms(virus_hits: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Return best STAT hit per organism key among evidence-passing non-background taxa."""
    found: dict[str, dict[str, Any]] = {}
    for hit in virus_hits:
        if not evidence_ok(hit):
            continue
        if is_background(hit["org"]):
            continue
        # Skip pure hierarchy nodes that are too broad unless species/genus/family useful.
        for key, pat in ORGANISM_PATTERNS:
            if pat.search(hit["org"]):
                prev = found.get(key)
                if prev is None or hit["total_count"] > prev["total_count"]:
                    found[key] = hit
    return found


def respiratory_hits(virus_hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = {
        "rsv",
        "influenza_a",
        "influenza_b",
        "sars_cov_2",
        "coronavirus_oc43",
        "coronavirus_hku1",
        "coronavirus_nl63",
        "coronavirus_229e",
        "rhinovirus",
        "enterovirus",
        "parainfluenza",
        "metapneumovirus",
        "adenovirus",
        "bocavirus",
    }
    matched = matched_organisms(virus_hits)
    return [
        {
            "organism_key": k,
            "org": v["org"],
            "tax_id": v["tax_id"],
            "self_count": v["self_count"],
            "total_count": v["total_count"],
            "rank": v["rank"],
        }
        for k, v in matched.items()
        if k in keys
    ]


def expected_supported(expected_key: str, found_map: dict[str, dict[str, Any]]) -> bool:
    """Return True if STAT evidence supports an expected clinical organism key."""
    if expected_key in found_map:
        return True
    # Allow SARS clinical expectations to be supported by SARS-related / coronaviridae nodes.
    if expected_key == "sars_cov_2":
        return "coronaviridae" in found_map or "sars_cov_2" in found_map
    # Rhinovirus/Enterovirus panel channel is often co-reported.
    if expected_key == "rhinovirus" and "enterovirus" in found_map:
        return True
    if expected_key == "enterovirus" and "rhinovirus" in found_map:
        return True
    return False


def classify_record(
    row: pd.Series,
    clinical: dict[str, Any] | None,
    virus_hits: list[dict[str, Any]],
    disc_by_sample: dict[str, list[dict[str, Any]]],
    fetch_error: str | None,
) -> dict[str, Any]:
    expected_pcr = tokenize_clinical(clinical.get("pcr_results") if clinical else None)
    expected_mngs = tokenize_clinical(clinical.get("mngs_results") if clinical else None)
    # Prefer mNGS declared detections for sequence comparison; fall back to PCR.
    expected = expected_mngs or expected_pcr
    found_map = matched_organisms(virus_hits)
    resp = respiratory_hits(virus_hits)
    found_keys = [h["organism_key"] for h in resp]

    article_id = row.get("article_sample_id")
    if isinstance(article_id, float) and pd.isna(article_id):
        article_id = ""
    article_id = str(article_id or "").strip()
    disc_rows = disc_by_sample.get(article_id, []) if article_id else []

    notes: list[str] = []
    flags: list[str] = []

    if fetch_error:
        flags.append("stat_fetch_failed")
        notes.append(fetch_error)
        consistency = "audit_incomplete"
        evidence_basis = "none"
    elif row["record_role"] == "negative_control":
        # Negative controls: respiratory pathogens at >=3 are suspicious.
        if resp:
            flags.append("sequence_metadata_conflict")
            notes.append(
                "Negative control has respiratory pathogen STAT evidence >=3: "
                + ", ".join(f"{h['organism_key']}({h['org']}:{h['total_count']})" for h in resp)
            )
            consistency = "inconsistent_sequence_vs_control_label"
        else:
            consistency = "consistent_negative_control"
            notes.append("No respiratory pathogen STAT evidence >=3.")
        evidence_basis = "sequence"
    elif row["record_role"] == "positive_control":
        if resp:
            consistency = "consistent_positive_control_virus_detected"
            notes.append(
                "Positive control shows respiratory/viral STAT evidence: "
                + ", ".join(f"{h['organism_key']}({h['total_count']})" for h in resp[:8])
            )
        else:
            # Still may have non-respiratory virus; check any non-background virus evidence.
            non_bg = [
                h
                for h in virus_hits
                if evidence_ok(h) and not is_background(h["org"]) and h["tax_id"] != VIRUS_ROOT
            ]
            if non_bg:
                consistency = "consistent_positive_control_virus_detected"
                notes.append(
                    "Positive control virus evidence present (non-panel taxa): "
                    + ", ".join(f"{h['org']}:{h['total_count']}" for h in non_bg[:5])
                )
            else:
                flags.append("sequence_metadata_conflict")
                consistency = "inconsistent_sequence_vs_control_label"
                notes.append("Positive control lacks STAT virus evidence >=3.")
        evidence_basis = "sequence"
    elif not expected:
        # No joined clinical organism expectation.
        if row.get("clinical_metadata_join") == "unresolved_NP_NUC_alias":
            consistency = "sequence_profiled_metadata_unresolved"
            notes.append(
                "NP/NUC identifier lacks published P-NP crosswalk; sequence profile retained without clinical comparison."
            )
        elif row.get("clinical_metadata_join") == "not_in_supplement":
            consistency = "sequence_profiled_no_supp_data1_join"
            notes.append("No Supplementary Data 1 organism labels joined for this record.")
        else:
            consistency = "sequence_profiled_no_expected_organisms"
            notes.append("No parseable expected organisms in joined clinical metadata.")
        if resp:
            notes.append(
                "STAT respiratory evidence: "
                + ", ".join(f"{h['organism_key']}({h['total_count']})" for h in resp)
            )
        evidence_basis = "sequence"
    else:
        missing = [k for k in expected if not expected_supported(k, found_map)]
        # Extra respiratory finds beyond expected (and beyond common off-panel co-detections
        # already declared in mNGS text) are informative but not automatic mislabels.
        unexpected = [k for k in found_keys if k not in expected]
        # Off-panel extras are not treated as RVP false positives / mislabels by themselves.
        off_panel = {"cmv", "herpesvirus", "papillomavirus", "bocavirus"}
        unexpected_panel = [k for k in unexpected if k not in off_panel]

        def _support_note(k: str) -> str:
            if k in found_map:
                h = found_map[k]
                return f"{k}({h['org']}:{h['total_count']})"
            if k == "sars_cov_2" and "coronaviridae" in found_map:
                h = found_map["coronaviridae"]
                return f"{k} via {h['org']}:{h['total_count']}"
            if k == "rhinovirus" and "enterovirus" in found_map:
                h = found_map["enterovirus"]
                return f"{k} via {h['org']}:{h['total_count']}"
            if k == "enterovirus" and "rhinovirus" in found_map:
                h = found_map["rhinovirus"]
                return f"{k} via {h['org']}:{h['total_count']}"
            return k

        if not missing:
            consistency = "consistent_sequence_supports_declared_organisms"
            notes.append(
                "STAT evidence supports declared organisms: "
                + ", ".join(_support_note(k) for k in expected)
            )
            if unexpected_panel:
                notes.append(
                    "Additional on-panel STAT detections not listed in joined metadata: "
                    + ", ".join(unexpected_panel)
                )
                flags.append("extra_on_panel_stat_detection")
        else:
            # If sample is a known article discrepancy, missing STAT support may still
            # align with mNGS-negative adjudication — flag carefully.
            if disc_rows:
                flags.append("article_reported_discrepancy")
                notes.append(
                    "Sample appears in reconstructed Supplementary Tables 6/7; "
                    "missing STAT support for: " + ", ".join(missing)
                )
                consistency = "article_discrepancy_sequence_checked"
            else:
                flags.append("possible_sequence_metadata_conflict")
                notes.append(
                    "Declared organisms lacking STAT evidence >=3: "
                    + ", ".join(missing)
                    + f"; found={found_keys or 'none'}"
                )
                near = []
                for h in virus_hits:
                    if is_background(h["org"]):
                        continue
                    if 0 < max(h["self_count"], h["total_count"]) < MIN_EVIDENCE:
                        for key, pat in ORGANISM_PATTERNS:
                            if key in missing and pat.search(h["org"]):
                                near.append(
                                    f"{h['org']}:{h['total_count']}/{h['self_count']}"
                                )
                                break
                if near:
                    notes.append("Below-threshold near-miss taxa: " + ", ".join(near[:8]))
                    flags.append("below_threshold_near_miss")
                consistency = "possible_inconsistent_sequence_vs_metadata"
        evidence_basis = "sequence"

    if disc_rows and "article_reported_discrepancy" not in flags:
        flags.append("article_reported_discrepancy")

    top_viruses = [
        {
            "org": h["org"],
            "tax_id": h["tax_id"],
            "self_count": h["self_count"],
            "total_count": h["total_count"],
            "rank": h["rank"],
            "background": is_background(h["org"]),
        }
        for h in virus_hits
        if evidence_ok(h) and h["tax_id"] != VIRUS_ROOT
    ][:15]

    return {
        "File_Name": row["File_Name"],
        "Run": row["Run"],
        "Original_ID": row["Original_ID"],
        "specimen_key": row["specimen_key"],
        "article_sample_id": article_id,
        "record_role": row["record_role"],
        "is_control": bool(row["is_control"]),
        "is_rerun_or_replicate": bool(row["is_rerun_or_replicate"]),
        "mapping_method": row["mapping_method"],
        "clinical_metadata_join": row["clinical_metadata_join"],
        "spots": row.get("spots"),
        "pcr_results": clinical.get("pcr_results") if clinical else "",
        "mngs_results": clinical.get("mngs_results") if clinical else "",
        "expected_organism_keys": "|".join(expected),
        "expected_from": (
            "mngs_results"
            if expected_mngs
            else ("pcr_results" if expected_pcr else "")
        ),
        "stat_respiratory_organism_keys": "|".join(found_keys),
        "stat_respiratory_hits_json": json.dumps(resp, sort_keys=True),
        "stat_top_virus_hits_json": json.dumps(top_viruses, sort_keys=True),
        "stat_virus_taxa_with_evidence": sum(1 for h in virus_hits if evidence_ok(h)),
        "consistency_class": consistency,
        "evidence_basis": evidence_basis,
        "flags": "|".join(flags),
        "notes": " ".join(notes),
        "article_discrepancy_tables": "|".join(
            sorted({f"T{d['table']}-r{d['table_row']}" for d in disc_rows})
        ),
        "article_final_mngs_by_dtca": " || ".join(
            str(d.get("final_mngs_by_dtca") or "") for d in disc_rows if d.get("final_mngs_by_dtca")
        ),
        "min_evidence_threshold": MIN_EVIDENCE,
        "evidence_rule_note": (
            "STAT self_count/total_count >=3 used as public proxy for "
            "article >=3 non-overlapping reads/contigs criterion; not identical."
        ),
    }


def load_disc_index() -> dict[str, list[dict[str, Any]]]:
    if not DISC_SAMPLE.exists():
        return {}
    df = pd.read_csv(DISC_SAMPLE)
    out: dict[str, list[dict[str, Any]]] = {}
    for rec in df.to_dict(orient="records"):
        sid = rec.get("article_sample_id")
        if sid is None or (isinstance(sid, float) and pd.isna(sid)):
            continue
        sid = str(sid).strip()
        if not sid:
            continue
        out.setdefault(sid, []).append(rec)
    return out


def cached_fetch(acc: str, sess: requests.Session, refresh: bool) -> tuple[str, dict[str, Any] | None, str | None]:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CACHE_DIR / f"{acc}.json"
    if path.exists() and not refresh:
        try:
            return acc, json.loads(path.read_text()), None
        except json.JSONDecodeError:
            pass
    try:
        obj = fetch_taxonomy(sess, acc)
        path.write_text(json.dumps(obj))
        return acc, obj, None
    except Exception as exc:  # noqa: BLE001
        return acc, None, str(exc)


def run_audit(workers: int, refresh: bool, limit: int | None) -> dict[str, Any]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(MANIFEST)
    clinical_df = pd.read_csv(CLINICAL)
    clinical_by_id = {
        str(r["article_sample_id"]): r
        for r in clinical_df.to_dict(orient="records")
    }
    disc_by_sample = load_disc_index()

    if limit is not None:
        manifest = manifest.head(limit)

    runs = manifest["Run"].astype(str).tolist()
    results_by_run: dict[str, tuple[dict[str, Any] | None, str | None]] = {}

    # Thread-local sessions via one session per worker task is fine; create per fetch.
    def _task(acc: str) -> tuple[str, dict[str, Any] | None, str | None]:
        return cached_fetch(acc, session(), refresh)

    print(f"Fetching STAT taxonomy for {len(runs)} Accuracy runs (workers={workers})...")
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_task, acc): acc for acc in runs}
        for fut in as_completed(futs):
            acc, obj, err = fut.result()
            results_by_run[acc] = (obj, err)
            done += 1
            if done % 25 == 0 or done == len(runs):
                print(f"  {done}/{len(runs)}")

    audit_rows: list[dict[str, Any]] = []
    hit_rows: list[dict[str, Any]] = []
    for _, row in manifest.iterrows():
        acc = str(row["Run"])
        obj, err = results_by_run.get(acc, (None, "missing_result"))
        viruses: list[dict[str, Any]] = []
        if obj is not None:
            viruses = virus_rows(obj)
            for h in viruses:
                if evidence_ok(h):
                    hit_rows.append(
                        {
                            "Run": acc,
                            "File_Name": row["File_Name"],
                            "article_sample_id": row.get("article_sample_id"),
                            "record_role": row["record_role"],
                            "tax_id": h["tax_id"],
                            "org": h["org"],
                            "rank": h["rank"],
                            "self_count": h["self_count"],
                            "total_count": h["total_count"],
                            "background": is_background(h["org"]),
                        }
                    )
        article_id = row.get("article_sample_id")
        if isinstance(article_id, float) and pd.isna(article_id):
            article_id = None
        clinical = clinical_by_id.get(str(article_id)) if article_id else None
        audit_rows.append(
            classify_record(row, clinical, viruses, disc_by_sample, err)
        )

    audit_df = pd.DataFrame(audit_rows)
    hits_df = pd.DataFrame(hit_rows)
    audit_df.to_csv(FINDINGS_CSV, index=False)
    hits_df.to_csv(VIRUS_HITS_CSV, index=False)

    class_counts = audit_df["consistency_class"].value_counts().to_dict()
    flag_counts: dict[str, int] = {}
    for flags in audit_df["flags"].fillna(""):
        for f in str(flags).split("|"):
            if f:
                flag_counts[f] = flag_counts.get(f, 0) + 1

    candidate_inconsistencies = audit_df[
        audit_df["consistency_class"].isin(
            [
                "inconsistent_sequence_vs_control_label",
                "possible_inconsistent_sequence_vs_metadata",
            ]
        )
        | audit_df["flags"].fillna("").str.contains(
            "sequence_metadata_conflict|possible_sequence_metadata_conflict",
            regex=True,
        )
    ]

    summary = {
        "generated_at_utc": utc_now(),
        "scope": "audit-sequences via NCBI STAT public run_taxonomy; no FASTQ download",
        "selection_rule": "all Purpose=Accuracy rows from accuracy_manifest.csv",
        "counts": {
            "accuracy_records_audited": int(len(audit_df)),
            "stat_fetch_success": int((audit_df["evidence_basis"] == "sequence").sum()),
            "stat_fetch_failed": int(
                audit_df["consistency_class"].eq("audit_incomplete").sum()
            ),
            "clinical_records": int((audit_df["record_role"] == "clinical_sample").sum()),
            "positive_controls": int((audit_df["record_role"] == "positive_control").sum()),
            "negative_controls": int((audit_df["record_role"] == "negative_control").sum()),
            "rerun_or_replicate_records": int(audit_df["is_rerun_or_replicate"].sum()),
            "records_with_article_discrepancy_flag": int(
                audit_df["flags"].fillna("").str.contains("article_reported_discrepancy").sum()
            ),
            "candidate_sequence_inconsistencies": int(len(candidate_inconsistencies)),
        },
        "consistency_class_counts": class_counts,
        "flag_counts": flag_counts,
        "evidence_rule": {
            "min_stat_count": MIN_EVIDENCE,
            "note": (
                "STAT self_count/total_count >=3 is a public proxy for the article "
                "criterion of >=3 non-overlapping reads/contigs; not identical."
            ),
        },
        "outputs": {
            "audit_csv": str(FINDINGS_CSV.relative_to(ROOT)),
            "audit_json": str(FINDINGS_JSON.relative_to(ROOT)),
            "virus_hits_csv": str(VIRUS_HITS_CSV.relative_to(ROOT)),
            "summary_json": str(SUMMARY_JSON.relative_to(ROOT)),
            "stat_cache_dir": str(CACHE_DIR.relative_to(ROOT)),
        },
        "limitations": [
            "STAT k-mer taxonomy counts are a proxy, not a reimplementation of the paper's read/contig pipeline.",
            "79 NP/NUC identifiers remain without published P-NP crosswalks; those runs are profiled but not clinically compared.",
            "Supplementary Data 1 covers 109 sample IDs; many Accuracy clinical records lack joined organism labels.",
            "Off-panel detections (CMV, papillomavirus, herpesviruses, bocavirus) are reported but not auto-labeled false positives.",
            "No FASTQ re-analysis was performed; ENA/SRA FASTQ fallback was not required for STAT coverage.",
        ],
        "candidate_inconsistency_runs": candidate_inconsistencies[
            ["Run", "File_Name", "Original_ID", "article_sample_id", "consistency_class", "notes"]
        ].to_dict(orient="records"),
    }

    FINDINGS_JSON.write_text(json.dumps({"summary": summary, "records": audit_rows}, indent=2))
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary["counts"], indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--refresh", action="store_true", help="Ignore STAT JSON cache")
    parser.add_argument("--limit", type=int, default=None, help="Audit only first N rows")
    args = parser.parse_args()
    run_audit(workers=args.workers, refresh=args.refresh, limit=args.limit)


if __name__ == "__main__":
    main()
