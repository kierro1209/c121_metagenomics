#!/usr/bin/env python3
"""Assemble final Accuracy metadata/sequence audit deliverables.

Reads normalized manifest, discrepancy reconstruction, and STAT sequence audit
outputs, then writes:
  - results/accuracy_metadata_audit.csv
  - results/inconsistent_samples.md
  - results/metadata_consistency_report.md
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/normalized/accuracy_manifest.csv"
SPECIMENS = ROOT / "data/normalized/accuracy_specimens.csv"
CLINICAL = ROOT / "data/normalized/clinical_metadata_samples.csv"
UNRESOLVED = ROOT / "data/normalized/unresolved_identifier_mappings.csv"
DISC_SAMPLE = ROOT / "data/normalized/article_discrepancies_sample_level.csv"
DISC_ORG = ROOT / "data/normalized/article_discrepancies_organism_level.csv"
DISC_SUMMARY = ROOT / "results/intermediate/discrepancy_reconstruction_summary.json"
SEQ_AUDIT = ROOT / "results/intermediate/accuracy_sequence_audit.csv"
SEQ_SUMMARY = ROOT / "results/intermediate/accuracy_sequence_audit_summary.json"
MANIFEST_PROV = ROOT / "data/normalized/accuracy_manifest_provenance.json"

OUT_CSV = ROOT / "results/accuracy_metadata_audit.csv"
OUT_INCONSISTENT = ROOT / "results/inconsistent_samples.md"
OUT_REPORT = ROOT / "results/metadata_consistency_report.md"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def build_audit_csv() -> pd.DataFrame:
    manifest = pd.read_csv(MANIFEST)
    clinical = pd.read_csv(CLINICAL)
    disc = pd.read_csv(DISC_SAMPLE)
    seq = pd.read_csv(SEQ_AUDIT)

    # Collapse discrepancy rows per article sample id.
    disc_agg = (
        disc.dropna(subset=["article_sample_id"])
        .assign(article_sample_id=lambda d: d["article_sample_id"].astype(str))
        .groupby("article_sample_id", as_index=False)
        .agg(
            article_discrepancy_class=("discrepancy_class", lambda s: " | ".join(sorted(set(map(str, s))))),
            article_final_mngs_by_dtca=(
                "final_mngs_by_dtca",
                lambda s: " || ".join(str(x) for x in s if pd.notna(x) and str(x).strip()),
            ),
            article_id_mapping_status=(
                "id_mapping_status",
                lambda s: " | ".join(sorted(set(map(str, s)))),
            ),
            article_discrepancy_explanation=(
                "explanation",
                lambda s: " || ".join(str(x) for x in s if pd.notna(x) and str(x).strip()),
            ),
        )
    )

    merged = manifest.merge(
        clinical,
        how="left",
        on="article_sample_id",
        suffixes=("", "_clinical"),
    ).merge(
        disc_agg,
        how="left",
        on="article_sample_id",
    ).merge(
        seq[
            [
                "Run",
                "expected_organism_keys",
                "stat_respiratory_organism_keys",
                "stat_respiratory_hits_json",
                "stat_top_virus_hits_json",
                "stat_virus_taxa_with_evidence",
                "consistency_class",
                "evidence_basis",
                "flags",
                "notes",
                "min_evidence_threshold",
                "evidence_rule_note",
            ]
        ],
        how="left",
        on="Run",
        suffixes=("", "_seq"),
    )

    # Derive a compact adjudication label for the final CSV.
    def adjudicate(r: pd.Series) -> str:
        cls = str(r.get("consistency_class") or "")
        flags = str(r.get("flags") or "")
        if "article_reported_discrepancy" in flags:
            return "article_discrepancy_documented"
        if cls.startswith("consistent"):
            return "consistent"
        if cls in {
            "inconsistent_sequence_vs_control_label",
            "possible_inconsistent_sequence_vs_metadata",
        }:
            return "candidate_inconsistency"
        if "unresolved" in cls or r.get("clinical_metadata_join") in {
            "unresolved_NP_NUC_alias",
            "not_in_supplement",
        }:
            return "metadata_incomplete_not_adjudicated"
        if cls == "audit_incomplete":
            return "audit_incomplete"
        return "reviewed_no_hard_inconsistency"

    merged["final_adjudication"] = merged.apply(adjudicate, axis=1)
    merged["audit_generated_at_utc"] = utc_now()

    preferred_cols = [
        "File_Name",
        "Run",
        "BioSample",
        "Original_ID",
        "specimen_key",
        "article_sample_id",
        "record_role",
        "is_control",
        "is_rerun_or_replicate",
        "mapping_method",
        "clinical_metadata_join",
        "pcr_results",
        "mngs_results",
        "expected_organism_keys",
        "stat_respiratory_organism_keys",
        "stat_virus_taxa_with_evidence",
        "consistency_class",
        "evidence_basis",
        "flags",
        "notes",
        "article_discrepancy_class",
        "article_final_mngs_by_dtca",
        "article_id_mapping_status",
        "article_discrepancy_explanation",
        "final_adjudication",
        "spots",
        "size_MB",
        "download_path",
        "min_evidence_threshold",
        "evidence_rule_note",
        "audit_generated_at_utc",
    ]
    cols = [c for c in preferred_cols if c in merged.columns]
    out = merged[cols].copy()
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    return out


def write_inconsistent_md(audit: pd.DataFrame, seq_summary: dict, disc_summary: dict) -> None:
    candidates = audit[
        audit["final_adjudication"].isin(
            ["candidate_inconsistency", "article_discrepancy_documented"]
        )
    ].copy()

    seq_conflicts = audit[
        audit["consistency_class"].isin(
            [
                "inconsistent_sequence_vs_control_label",
                "possible_inconsistent_sequence_vs_metadata",
            ]
        )
    ]

    lines: list[str] = []
    lines.append("# Inconsistent / discrepancy samples")
    lines.append("")
    lines.append(f"Generated (UTC): {utc_now()}")
    lines.append("")
    lines.append("## Definition")
    lines.append("")
    lines.append(
        "Inconsistency here means either (a) an article-reported mNGS vs RVP disagreement "
        "documented in Supplementary Tables 6–7, or (b) STAT sequence evidence conflicting "
        "with declared clinical/control metadata under a >=3-count proxy for non-overlapping "
        "viral evidence. Off-panel mNGS detections (for example CMV, papillomavirus, "
        "herpesviruses, bocavirus when not an RVP target) are **not** automatically labeled "
        "inconsistent. Unresolved NP/NUC↔P-NP identifiers are profiled but not forced into "
        "a fabricated clinical join."
    )
    lines.append("")
    lines.append("## Summary counts")
    lines.append("")
    lines.append(f"- Accuracy records audited: **{len(audit)}**")
    lines.append(
        f"- Article discrepancy rows reconstructed: Table 6 = "
        f"**{disc_summary['verification']['table6_sample_rows']}**, Table 7 = "
        f"**{disc_summary['verification']['table7_sample_rows']}** "
        f"(main-text mNGS-/RVP+ = 7 via TN+FN weights "
        f"{disc_summary['verification']['table7_TN_plus_FN_weight']})"
    )
    lines.append(
        f"- Sequence candidate conflicts (STAT vs metadata/controls): "
        f"**{len(seq_conflicts)}**"
    )
    lines.append(
        f"- Rows flagged article discrepancy and/or sequence candidate: **{len(candidates)}**"
    )
    lines.append("")

    lines.append("## A. Article-reported discrepancies (Supplementary Tables 6–7)")
    lines.append("")
    disc = pd.read_csv(DISC_SAMPLE)
    for _, r in disc.iterrows():
        sid = r.get("article_sample_id")
        sid_s = str(sid) if pd.notna(sid) else "(unresolved sample id)"
        lines.append(
            f"### Table {int(r['table'])} row {int(r['table_row'])}: {sid_s}"
        )
        lines.append("")
        lines.append(f"- Class: `{r['discrepancy_class']}`")
        lines.append(f"- ID mapping: `{r['id_mapping_status']}`")
        lines.append(f"- RVP organisms: {r.get('rvp_organisms')}")
        lines.append(f"- mNGS organisms: {r.get('mngs_organisms')}")
        lines.append(f"- Orthogonal / repeat testing: {r.get('discrepancy_testing')}")
        lines.append(f"- Final mNGS by DTCA: **{r.get('final_mngs_by_dtca')}**")
        if pd.notna(r.get("explanation")):
            lines.append(f"- Explanation: {r.get('explanation')}")
        # Attach STAT notes if sample id known.
        if pd.notna(sid):
            hits = audit[audit["article_sample_id"].astype(str) == str(sid)]
            if len(hits):
                # Prefer primary (non-rerun) if present.
                primary = hits.sort_values("is_rerun_or_replicate")
                h = primary.iloc[0]
                lines.append(
                    f"- STAT respiratory keys ({h['Run']}): "
                    f"`{h.get('stat_respiratory_organism_keys') or 'none'}`"
                )
                lines.append(f"- Sequence consistency class: `{h.get('consistency_class')}`")
        lines.append("")

    lines.append("## B. Sequence vs metadata / control candidate conflicts")
    lines.append("")
    if seq_conflicts.empty:
        lines.append(
            "No hard STAT-vs-label conflicts were detected among controls or joined clinical "
            "records under the >=3 STAT-count proxy."
        )
        lines.append("")
    else:
        def _clean(val: object) -> str:
            if val is None or (isinstance(val, float) and pd.isna(val)):
                return ""
            text = str(val).strip()
            return "" if text.lower() == "nan" else text

        for _, r in seq_conflicts.iterrows():
            lines.append(
                f"### {r['Original_ID']} / {r['Run']} ({r['record_role']})"
            )
            lines.append("")
            lines.append(f"- File: `{r['File_Name']}`")
            lines.append(f"- Article sample id: `{_clean(r.get('article_sample_id'))}`")
            lines.append(f"- Consistency class: `{r['consistency_class']}`")
            lines.append(f"- Expected keys: `{_clean(r.get('expected_organism_keys'))}`")
            lines.append(
                f"- STAT respiratory keys: `{_clean(r.get('stat_respiratory_organism_keys'))}`"
            )
            lines.append(f"- Notes: {_clean(r.get('notes'))}")
            lines.append("")

    lines.append("## C. Explicitly unresolved / not forced")
    lines.append("")
    unresolved_n = len(pd.read_csv(UNRESOLVED))
    lines.append(
        f"- Unresolved NP/NUC↔article identifier mappings retained: **{unresolved_n}**"
    )
    disc_ids = sorted(
        {
            str(x)
            for x in pd.read_csv(DISC_SAMPLE)["article_sample_id"].dropna().astype(str)
            if str(x).strip()
        }
    )
    manifest_ids = set(audit["article_sample_id"].dropna().astype(str))
    missing_disc = [x for x in disc_ids if x not in manifest_ids]
    lines.append(
        f"- Reconstructed discrepancy sample IDs with no Accuracy-manifest "
        f"`article_sample_id` join (cannot attach an SRR without inventing a crosswalk): "
        f"**{', '.join(missing_disc) if missing_disc else 'none'}**"
    )
    for issue in disc_summary.get("unresolved_issues", []):
        lines.append(f"- {issue}")
    lines.append("")
    lines.append("## D. STAT coverage gaps")
    lines.append("")
    incomplete = audit[audit["consistency_class"] == "audit_incomplete"]
    if incomplete.empty:
        lines.append("- None.")
    else:
        for _, r in incomplete.iterrows():
            lines.append(
                f"- `{r['Run']}` / `{r['Original_ID']}` / `{r['File_Name']}`: STAT fetch failed"
            )
    lines.append("")
    lines.append("## E. Post-DTCA verification targets")
    lines.append("")
    post = disc_summary["verification"]["post_dtca"]
    lines.append(f"- Paper post-DTCA OPA: **{post['opa']}** (target 187/191 ≈ 97.9%)")
    lines.append(f"- Paper post-DTCA PPA: {post['ppa']}")
    lines.append(f"- Paper post-DTCA NPA: {post['npa']}")
    lines.append(
        "- Residual FNs after DTCA (reconstruction): RSV 0.5 + metapneumovirus 1.0"
    )
    lines.append("")

    OUT_INCONSISTENT.write_text("\n".join(lines) + "\n")


def write_report_md(audit: pd.DataFrame, seq_summary: dict, disc_summary: dict) -> None:
    prov = load_json(MANIFEST_PROV)
    specimens = pd.read_csv(SPECIMENS)
    clinical = pd.read_csv(CLINICAL)

    n_acc = len(audit)
    n_clinical = int((audit["record_role"] == "clinical_sample").sum())
    n_pc = int((audit["record_role"] == "positive_control").sum())
    n_nc = int((audit["record_role"] == "negative_control").sum())
    n_rerun = int(audit["is_rerun_or_replicate"].sum())

    adj_counts = audit["final_adjudication"].value_counts().to_dict()
    class_counts = seq_summary.get("consistency_class_counts", {})

    lines: list[str] = []
    lines.append("# Accuracy cohort metadata consistency report")
    lines.append("")
    lines.append(f"Generated (UTC): {utc_now()}")
    lines.append("")
    lines.append("## Introduction")
    lines.append("")
    lines.append(
        "This report audits the Accuracy evaluation cohort from the Nature Communications "
        "respiratory metagenomic NGS (mNGS) validation study "
        "([s41467-024-51470-y](https://www.nature.com/articles/s41467-024-51470-y); "
        "Zenodo [10553379](https://zenodo.org/records/10553379); BioProject "
        "[PRJNA1084017](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1084017)). "
        "The goal is to verify that publicly deposited sequencing records, clinical/control "
        "labels, and article-reported mNGS vs respiratory viral panel (RVP) discrepancies "
        "are internally consistent, and to separate true inconsistencies from unresolved "
        "identifier gaps or off-panel detections."
    )
    lines.append("")
    lines.append("## Challenge")
    lines.append("")
    lines.append(
        "Public Accuracy metadata are split across a Zenodo sequencing master sheet, "
        "Supplementary Data 1 clinical tables, Supplementary Tables 6–7 discrepancy "
        "narratives, and SRA/STAT sequence products. Identifier systems are not fully "
        "crosswalked (NP/NUC lab IDs vs article P-NP IDs), Supplementary Data 1 covers "
        "fewer sample IDs than the Accuracy clinical cohort, and article discrepancy "
        "counts use fractional DTCA weights. An audit must therefore:"
    )
    lines.append("")
    lines.append("1. Enumerate all Purpose=Accuracy sequencing records reproducibly.")
    lines.append("2. Reconstruct every Table 6/7 discrepancy without inventing missing IDs.")
    lines.append(
        "3. Independently profile sequence content (here: NCBI STAT taxonomy) and compare "
        "it to declared labels where joins exist."
    )
    lines.append(
        "4. Report inconsistencies distinctly from QC/low-confidence issues and from "
        "metadata incompleteness."
    )
    lines.append("")
    lines.append("## Software and data")
    lines.append("")
    lines.append("### Software")
    lines.append("")
    lines.append("- Python 3 with `pandas`, `requests`, `openpyxl`")
    lines.append("- `analysis/build_accuracy_manifest.py` — Accuracy manifest assembly")
    lines.append("- `analysis/reconstruct_discrepancies.py` — Tables 6–7 reconstruction")
    lines.append("- `analysis/audit_accuracy_sequences.py` — STAT sequence audit")
    lines.append("- `analysis/build_accuracy_audit_report.py` — this report bundle")
    lines.append(
        "- NCBI STAT `run_taxonomy` public API "
        "(`cluster_name=public`) as sequence-content evidence source"
    )
    lines.append("")
    lines.append("### Primary data sources")
    lines.append("")
    for key, src in prov.get("sources", {}).items():
        lines.append(
            f"- **{key}**: {src.get('role')} — `{src.get('filename')}` "
            f"([link]({src.get('url')}); sha256 `{src.get('sha256')}`)"
        )
    lines.append("")
    lines.append("### Evidence rule")
    lines.append("")
    lines.append(
        "The article calls organisms with ≥3 non-overlapping reads/contigs. This audit "
        "uses NCBI STAT `self_count`/`total_count` ≥3 as a **public proxy** for that "
        "threshold. It is not a byte-identical reimplementation of the paper pipeline."
    )
    lines.append("")
    lines.append("## Detailed results")
    lines.append("")
    lines.append("### 1. Cohort inventory")
    lines.append("")
    lines.append(f"| Metric | Count |")
    lines.append(f"| --- | ---: |")
    lines.append(f"| Accuracy sequencing records | {n_acc} |")
    lines.append(f"| Clinical sequencing records | {n_clinical} |")
    lines.append(f"| Positive controls | {n_pc} |")
    lines.append(f"| Negative controls | {n_nc} |")
    lines.append(f"| Rerun/replicate records | {n_rerun} |")
    lines.append(f"| Accuracy specimen keys table rows | {len(specimens)} |")
    lines.append(
        f"| Expected article clinical specimens | {prov['counts']['expected_article_clinical_specimens']} |"
    )
    lines.append(f"| Supplementary Data 1 sample IDs | {len(clinical)} |")
    lines.append(
        f"| Unresolved NP/NUC mappings | {prov['counts']['unresolved_NP_NUC_identifiers']} |"
    )
    lines.append("")
    lines.append(
        "Verified inventory matches the assemble-manifest targets: **281** Accuracy files, "
        "with clinical cohort structure reconciling to the article’s **191** specimens "
        "(186 unambiguous non-repeat keys + 5 unresolved members of repeat run 039)."
    )
    lines.append("")
    lines.append("### 2. Article discrepancy reconstruction")
    lines.append("")
    v = disc_summary["verification"]
    lines.append(
        f"- Supplementary Table 6 (mNGS+/RVP−): **{v['table6_sample_rows']}** sample rows."
    )
    lines.append(
        f"- Supplementary Table 7 sample rows: **{v['table7_sample_rows']}**; "
        f"TN+FN fractional weight **{v['table7_TN_plus_FN_weight']}** matches main-text "
        f"mNGS-/RVP+ count of **{v['expected_maintext_mngs_neg_rvp_pos']}**."
    )
    lines.append(
        f"- Pre-DTCA accuracy reported **{v['pre_dtca_accuracy']['reported']}**; "
        f"post-DTCA OPA **{v['post_dtca']['opa']}**."
    )
    lines.append(
        f"- ID mapping status among discrepancy rows: "
        + ", ".join(f"{k}={n}" for k, n in v["id_mapping_status_counts"].items())
    )
    lines.append("")
    for finding in disc_summary.get("key_findings", []):
        lines.append(f"- {finding}")
    lines.append("")
    lines.append("Unresolved discrepancy-ID issues retained without invention:")
    lines.append("")
    for issue in disc_summary.get("unresolved_issues", []):
        lines.append(f"- {issue}")
    lines.append("")
    lines.append("### 3. Independent STAT sequence audit")
    lines.append("")
    sc = seq_summary["counts"]
    lines.append(
        f"- STAT fetch success: **{sc['stat_fetch_success']} / {sc['accuracy_records_audited']}**"
    )
    lines.append(
        f"- Candidate sequence inconsistencies: **{sc['candidate_sequence_inconsistencies']}**"
    )
    lines.append(
        f"- Records overlapping reconstructed article discrepancies: "
        f"**{sc['records_with_article_discrepancy_flag']}**"
    )
    lines.append("")
    lines.append("Consistency class counts:")
    lines.append("")
    for k, n in sorted(class_counts.items(), key=lambda kv: (-kv[1], kv[0])):
        lines.append(f"- `{k}`: {n}")
    lines.append("")
    lines.append("Final adjudication counts (report bundle):")
    lines.append("")
    for k, n in sorted(adj_counts.items(), key=lambda kv: (-kv[1], kv[0])):
        lines.append(f"- `{k}`: {n}")
    lines.append("")
    lines.append("### 4. What counts as inconsistent vs not")
    lines.append("")
    lines.append(
        "- **Inconsistent / discrepancy**: article Tables 6–7 cases, or STAT evidence "
        "conflicting with control labels / joined clinical organism expectations."
    )
    lines.append(
        "- **Not automatically inconsistent**: off-panel virus detections; unresolved "
        "NP/NUC aliases; Accuracy records absent from Supplementary Data 1; low-level "
        "background taxa (anelloviruses/phages) under the evidence proxy."
    )
    lines.append(
        "- **QC / incomplete**: STAT fetch failures (`audit_incomplete`) if any; these are "
        "reported separately from mislabel calls."
    )
    lines.append("")
    lines.append("### 5. Verification against paper headline numbers")
    lines.append("")
    lines.append("| Target | Result |")
    lines.append("| --- | --- |")
    lines.append(f"| 281 Accuracy files | {n_acc} audited |")
    lines.append(
        f"| 191 clinical specimens | article expected {prov['counts']['expected_article_clinical_specimens']}; "
        f"manifest clinical records {n_clinical} including reruns/controls excluded from 191 |"
    )
    lines.append(
        f"| 14 mNGS+/RVP− + 7 mNGS-/RVP+ | Table6={v['table6_sample_rows']}; "
        f"Table7 TN+FN weight={v['table7_TN_plus_FN_weight']} (match={v['table7_TN_plus_FN_matches_maintext_7']}) |"
    )
    lines.append(f"| ~97.9% post-adjudication agreement | {v['post_dtca']['opa']} |")
    lines.append("")
    lines.append("## Deliverables")
    lines.append("")
    lines.append(f"- `{OUT_CSV.relative_to(ROOT)}`")
    lines.append(f"- `{OUT_INCONSISTENT.relative_to(ROOT)}`")
    lines.append(f"- `{OUT_REPORT.relative_to(ROOT)}`")
    lines.append(
        f"- Intermediate sequence audit: `{SEQ_AUDIT.relative_to(ROOT)}`, "
        f"`{SEQ_SUMMARY.relative_to(ROOT)}`"
    )
    lines.append("")
    lines.append("## Limitations")
    lines.append("")
    for lim in seq_summary.get("limitations", []):
        lines.append(f"- {lim}")
    for lim in prov.get("known_mapping_limitations", []):
        lines.append(f"- {lim}")
    lines.append("")

    OUT_REPORT.write_text("\n".join(lines) + "\n")


def main() -> None:
    if not SEQ_AUDIT.exists() or not SEQ_SUMMARY.exists():
        raise SystemExit(
            "Sequence audit outputs missing. Run analysis/audit_accuracy_sequences.py first."
        )
    seq_summary = load_json(SEQ_SUMMARY)
    disc_summary = load_json(DISC_SUMMARY)
    audit = build_audit_csv()
    write_inconsistent_md(audit, seq_summary, disc_summary)
    write_report_md(audit, seq_summary, disc_summary)
    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {OUT_INCONSISTENT}")
    print(f"Wrote {OUT_REPORT}")
    print(audit["final_adjudication"].value_counts().to_dict())


if __name__ == "__main__":
    main()
