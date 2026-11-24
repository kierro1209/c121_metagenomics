#!/usr/bin/env python3
"""Inspect STAT virus calls for a known positive and a control."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
manifest = pd.read_csv(ROOT / "data/normalized/accuracy_manifest.csv")
clinical = pd.read_csv(ROOT / "data/normalized/clinical_metadata_samples.csv")

# Prefer a clinical sample with joined article metadata and a clear PCR call.
joined = manifest.merge(
    clinical,
    how="left",
    left_on="article_sample_id",
    right_on="article_sample_id",
)
candidates = joined[
    joined["pcr_results"].notna()
    & joined["is_control"].astype(str).str.lower().eq("false")
].head(5)
print("sample candidates:")
print(
    candidates[
        ["File_Name", "Run", "article_sample_id", "pcr_results", "mngs_results", "spots"]
    ].to_string(index=False)
)

controls = manifest[manifest["record_role"].isin(["positive_control", "negative_control"])].head(4)
print("\ncontrols:")
print(controls[["File_Name", "Run", "Original_ID", "record_role", "spots"]].to_string(index=False))

session = requests.Session()
session.headers.update({"User-Agent": "metagenomics-audit/1.0"})


def fetch(acc: str) -> dict:
    url = (
        "https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy"
        f"?acc={acc}&cluster_name=public"
    )
    r = session.get(url, timeout=180)
    r.raise_for_status()
    data = r.json()
    return data[0] if isinstance(data, list) else data


def virus_rows(obj: dict) -> list[dict]:
    table = obj.get("tax_table") or []
    # Find Viruses node and descendants via parent links.
    by_id = {int(x["tax_id"]): x for x in table if "tax_id" in x}
    virus_ids = set()
    # tax_id 10239 = Viruses
    if 10239 in by_id:
        virus_ids.add(10239)
    # Walk children by parent pointer
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
    rows = []
    for tid in virus_ids:
        row = by_id[tid]
        rows.append(
            {
                "tax_id": tid,
                "org": row.get("org"),
                "rank": row.get("rank"),
                "self_count": row.get("self_count", 0),
                "total_count": row.get("total_count", 0),
            }
        )
    rows.sort(key=lambda x: (-(x["total_count"] or 0), x["org"] or ""))
    return rows


for label, subset in [("clinical", candidates.head(2)), ("control", controls.head(2))]:
    for _, row in subset.iterrows():
        acc = row["Run"]
        print(f"\n=== {label} {acc} {row.get('article_sample_id') or row.get('Original_ID')} ===")
        print("expected pcr:", row.get("pcr_results"))
        print("expected mngs:", row.get("mngs_results"))
        obj = fetch(acc)
        print("tax_totals", obj.get("tax_totals"))
        viruses = virus_rows(obj)
        print("virus taxa", len(viruses))
        for v in viruses[:25]:
            if (v["self_count"] or 0) >= 1 or (v["total_count"] or 0) >= 3:
                print(v)
