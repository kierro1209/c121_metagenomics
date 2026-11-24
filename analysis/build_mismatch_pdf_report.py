#!/usr/bin/env python3
"""Build the Accuracy metadata-vs-sequence mismatch PDF report.

Report structure required by the project brief:
  1. Introduction to the research area and the problem
  2. Description of the challenge
  3. Description of the software and data used
  4. Details of the results

Primary finding for this document: which sample IDs have declared metadata
that do not match their actual STAT sequence content. The published article
is used as problem-space context, not as the primary mismatch list.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    KeepTogether,
    ListFlowable,
    ListItem,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
AUDIT_CSV = ROOT / "results/accuracy_metadata_audit.csv"
OUT_PDF = ROOT / "results/metadata_consistency_report.pdf"
OUT_MISMATCH_CSV = ROOT / "results/metadata_sequence_mismatches.csv"


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")


def clean(val: object) -> str:
    if val is None or (isinstance(val, float) and pd.isna(val)):
        return ""
    text = str(val).strip()
    return "" if text.lower() == "nan" else text


def styles():
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "TitleCustom",
            parent=base["Title"],
            fontSize=16,
            leading=20,
            spaceAfter=8,
            alignment=TA_CENTER,
        ),
        "subtitle": ParagraphStyle(
            "SubtitleCustom",
            parent=base["Normal"],
            fontSize=10,
            leading=13,
            textColor=colors.HexColor("#444444"),
            alignment=TA_CENTER,
            spaceAfter=18,
        ),
        "h1": ParagraphStyle(
            "H1Custom",
            parent=base["Heading1"],
            fontSize=13,
            leading=16,
            spaceBefore=14,
            spaceAfter=8,
            textColor=colors.HexColor("#111111"),
        ),
        "h2": ParagraphStyle(
            "H2Custom",
            parent=base["Heading2"],
            fontSize=11,
            leading=14,
            spaceBefore=10,
            spaceAfter=6,
            textColor=colors.HexColor("#222222"),
        ),
        "body": ParagraphStyle(
            "BodyCustom",
            parent=base["Normal"],
            fontSize=9.5,
            leading=13,
            alignment=TA_JUSTIFY,
            spaceAfter=6,
        ),
        "bullet": ParagraphStyle(
            "BulletCustom",
            parent=base["Normal"],
            fontSize=9.5,
            leading=12.5,
            leftIndent=8,
            spaceAfter=3,
        ),
        "small": ParagraphStyle(
            "SmallCustom",
            parent=base["Normal"],
            fontSize=8,
            leading=10,
            textColor=colors.HexColor("#333333"),
        ),
        "tablecell": ParagraphStyle(
            "TableCell",
            parent=base["Normal"],
            fontSize=7.5,
            leading=9.5,
            alignment=TA_LEFT,
        ),
        "tableheader": ParagraphStyle(
            "TableHeader",
            parent=base["Normal"],
            fontSize=8,
            leading=10,
            alignment=TA_LEFT,
            textColor=colors.white,
        ),
    }


def mismatch_frame(audit: pd.DataFrame) -> pd.DataFrame:
    hard = audit[
        audit["consistency_class"].isin(
            [
                "inconsistent_sequence_vs_control_label",
                "possible_inconsistent_sequence_vs_metadata",
            ]
        )
    ].copy()

    def severity(row: pd.Series) -> str:
        cls = str(row["consistency_class"])
        notes = str(row.get("notes") or "")
        if cls == "inconsistent_sequence_vs_control_label":
            if "sars_cov_2" in notes.lower() or "20" in notes:
                return "strong_mismatch"
            return "borderline_mismatch"
        if "Below-threshold near-miss" in notes:
            return "low_confidence_near_miss"
        return "possible_mismatch"

    hard["mismatch_severity"] = hard.apply(severity, axis=1)
    hard["declared_metadata"] = hard.apply(
        lambda r: clean(r.get("mngs_results"))
        or clean(r.get("pcr_results"))
        or ("negative_control" if r["record_role"] == "negative_control" else clean(r.get("record_role"))),
        axis=1,
    )
    hard["observed_sequence_signal"] = hard.apply(
        lambda r: clean(r.get("stat_respiratory_organism_keys"))
        or ("none >=3 STAT counts" if "lacking STAT" in str(r.get("notes") or "") else ""),
        axis=1,
    )
    cols = [
        "Original_ID",
        "article_sample_id",
        "Run",
        "BioSample",
        "File_Name",
        "record_role",
        "declared_metadata",
        "pcr_results",
        "mngs_results",
        "observed_sequence_signal",
        "stat_respiratory_organism_keys",
        "mismatch_severity",
        "consistency_class",
        "notes",
    ]
    # BioSample may be absent depending on merge; keep optional.
    keep = [c for c in cols if c in hard.columns]
    out = hard[keep].sort_values(["mismatch_severity", "Original_ID", "Run"])
    return out


def table(data: list[list[str]], col_widths: list[float], sty) -> Table:
    wrapped = []
    for i, row in enumerate(data):
        style = sty["tableheader"] if i == 0 else sty["tablecell"]
        wrapped.append([Paragraph(str(c), style) for c in row])
    t = Table(wrapped, colWidths=col_widths, repeatRows=1)
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1f2937")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#f8fafc")),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cbd5e1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [colors.HexColor("#ffffff"), colors.HexColor("#f1f5f9")],
                ),
            ]
        )
    )
    return t


def build() -> Path:
    audit = pd.read_csv(AUDIT_CSV)
    mismatches = mismatch_frame(audit)
    OUT_MISMATCH_CSV.parent.mkdir(parents=True, exist_ok=True)
    mismatches.to_csv(OUT_MISMATCH_CSV, index=False)

    sty = styles()
    doc = SimpleDocTemplate(
        str(OUT_PDF),
        pagesize=letter,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.7 * inch,
        bottomMargin=0.7 * inch,
        title="Accuracy Metadata Consistency Report",
        author="c121_metagenomics audit",
    )

    story: list = []
    story.append(Paragraph("Accuracy Cohort Metadata Consistency Report", sty["title"]))
    story.append(
        Paragraph(
            "Which sample metadata do not match their actual sequence content?<br/>"
            f"Generated {utc_now()} &nbsp;|&nbsp; Cohort: Purpose=Accuracy "
            "(Nature Communications s41467-024-51470-y context)",
            sty["subtitle"],
        )
    )

    # ---- 1. Introduction ----
    story.append(Paragraph("1. Introduction", sty["h1"]))
    story.append(
        Paragraph(
            "Respiratory metagenomic next-generation sequencing (mNGS) is used to detect "
            "pathogens directly from clinical specimens. In the Nature Communications "
            "validation study (DOI 10.1038/s41467-024-51470-y; Zenodo 10553379; "
            "BioProject PRJNA1084017), an Accuracy cohort of deposited FASTQ records is "
            "labeled with clinical or control metadata and compared against a commercial "
            "respiratory viral panel (RVP).",
            sty["body"],
        )
    )
    story.append(
        Paragraph(
            "This report answers a narrower reproducibility question: for publicly "
            "available Accuracy runs, which sample IDs have declared metadata that do "
            "not match the organisms actually detectable in the deposited sequence data? "
            "The article provides the problem-space context (cohort design, discrepancy "
            "tables, and adjudication math). The primary evidence here is an independent "
            "NCBI STAT taxonomy audit of each Accuracy SRR.",
            sty["body"],
        )
    )

    # ---- 2. Challenge ----
    story.append(Paragraph("2. Challenge", sty["h1"]))
    story.append(
        Paragraph(
            "Metadata and sequence products live in different public artifacts. Labels "
            "come from a Zenodo sequencing master sheet and Supplementary Data 1; "
            "sequence content is in SRA and NCBI STAT. Identifier systems are incomplete "
            "(79 NP/NUC lab IDs lack a published P-NP crosswalk), and Supplementary Data 1 "
            "covers fewer sample IDs than the Accuracy clinical cohort. An audit must "
            "therefore compare sequence evidence to labels only where joins exist, keep "
            "unresolved IDs flagged rather than invented, and separate true metadata/"
            "sequence conflicts from low-confidence near-misses and QC gaps.",
            sty["body"],
        )
    )
    story.append(
        Paragraph(
            "Inconsistency definition used here: sequence content conflicting with "
            "declared clinical or control metadata under a STAT count ≥3 proxy for the "
            "article’s ≥3 non-overlapping reads/contigs rule. Off-panel detections "
            "(for example CMV, papillomavirus, herpesviruses) are not automatic "
            "mislabels. Article-reported mNGS vs RVP disagreements are contextual "
            "background, not the mismatch list below.",
            sty["body"],
        )
    )

    # ---- 3. Software and data ----
    story.append(Paragraph("3. Software and data", sty["h1"]))
    story.append(Paragraph("Software", sty["h2"]))
    story.append(
        ListFlowable(
            [
                ListItem(Paragraph("Python 3 with pandas, requests, openpyxl, reportlab", sty["bullet"])),
                ListItem(Paragraph("analysis/build_accuracy_manifest.py — Accuracy manifest", sty["bullet"])),
                ListItem(Paragraph("analysis/reconstruct_discrepancies.py — Tables 6–7 context", sty["bullet"])),
                ListItem(Paragraph("analysis/audit_accuracy_sequences.py — STAT sequence audit", sty["bullet"])),
                ListItem(Paragraph("analysis/build_mismatch_pdf_report.py — this PDF", sty["bullet"])),
                ListItem(Paragraph("NCBI STAT run_taxonomy API (cluster_name=public)", sty["bullet"])),
            ],
            bulletType="bullet",
            leftIndent=12,
        )
    )
    story.append(Paragraph("Data sources", sty["h2"]))
    story.append(
        ListFlowable(
            [
                ListItem(
                    Paragraph(
                        "Zenodo sequencing master sheet (Purpose=Accuracy filter): "
                        "respmNGS-sequencing-mastersheet.xlsx",
                        sty["bullet"],
                    )
                ),
                ListItem(
                    Paragraph(
                        "Supplementary Data 1 clinical/lab metadata: "
                        "41467_2024_51470_MOESM4_ESM.xlsx",
                        sty["bullet"],
                    )
                ),
                ListItem(
                    Paragraph(
                        "Supplementary Information (Tables 6–7 context only): "
                        "41467_2024_51470_MOESM1_ESM.pdf",
                        sty["bullet"],
                    )
                ),
                ListItem(
                    Paragraph(
                        "SRA runinfo / BioProject PRJNA1084017 accessions for all Accuracy SRRs",
                        sty["bullet"],
                    )
                ),
            ],
            bulletType="bullet",
            leftIndent=12,
        )
    )
    story.append(
        Paragraph(
            "Evidence rule: STAT self_count/total_count ≥3 is used as a public proxy for "
            "the article criterion of ≥3 non-overlapping reads/contigs. It is not a "
            "byte-identical reimplementation of the paper pipeline.",
            sty["body"],
        )
    )

    # ---- 4. Results ----
    story.append(Paragraph("4. Detailed results", sty["h1"]))
    story.append(Paragraph("4.1 Cohort inventory", sty["h2"]))
    inv = [
        ["Metric", "Count"],
        ["Accuracy sequencing records audited", "281"],
        ["STAT taxonomy fetch success", "279 / 281"],
        ["Clinical sequencing records", "233"],
        ["Positive controls", "24"],
        ["Negative controls", "24"],
        ["Article clinical specimen target", "191"],
        ["Metadata↔sequence mismatch IDs (this report)", str(len(mismatches))],
    ]
    story.append(table(inv, [4.2 * inch, 2.0 * inch], sty))
    story.append(Spacer(1, 8))
    story.append(
        Paragraph(
            "Verification anchors from the project brief: 281 Accuracy files audited; "
            "191 clinical specimens as the article cohort target; article Tables 6–7 "
            "reconcile to 14 mNGS+/RVP− and 7 mNGS−/RVP+ (via TN+FN weight 7.0); "
            "post-DTCA overall agreement 187/191 ≈ 97.9%. Those article figures "
            "contextualize expected residual disagreements after adjudication; they are "
            "not themselves sequence-vs-metadata mismatches.",
            sty["body"],
        )
    )

    story.append(Paragraph("4.2 Sample IDs whose metadata do not match actual sequence data", sty["h2"]))
    story.append(
        Paragraph(
            f"<b>{len(mismatches)} Accuracy records</b> have declared metadata that conflict "
            "with STAT sequence evidence (or lack expected pathogen evidence under the ≥3 "
            "count proxy). IDs:",
            sty["body"],
        )
    )

    # Explicit ID callout list
    id_lines = []
    for _, r in mismatches.iterrows():
        sid = clean(r.get("article_sample_id")) or clean(r.get("Original_ID"))
        id_lines.append(
            f"<b>{sid}</b> &nbsp;|&nbsp; Original_ID={clean(r.get('Original_ID'))} "
            f"&nbsp;|&nbsp; Run={clean(r.get('Run'))} "
            f"&nbsp;|&nbsp; severity={clean(r.get('mismatch_severity'))}"
        )
    story.append(
        ListFlowable(
            [ListItem(Paragraph(x, sty["bullet"])) for x in id_lines],
            bulletType="bullet",
            leftIndent=12,
        )
    )

    story.append(Paragraph("Mismatch detail table", sty["h2"]))
    header = [
        "Sample ID",
        "Run",
        "Role",
        "Declared metadata",
        "Observed STAT signal",
        "Severity",
        "Why inconsistent",
    ]
    rows = [header]
    for _, r in mismatches.iterrows():
        sid = clean(r.get("article_sample_id")) or clean(r.get("Original_ID"))
        rows.append(
            [
                sid,
                clean(r.get("Run")),
                clean(r.get("record_role")),
                clean(r.get("declared_metadata")),
                clean(r.get("observed_sequence_signal")) or "none ≥3",
                clean(r.get("mismatch_severity")),
                clean(r.get("notes")),
            ]
        )
    story.append(
        table(
            rows,
            [
                0.75 * inch,
                0.95 * inch,
                0.85 * inch,
                1.05 * inch,
                0.95 * inch,
                0.95 * inch,
                1.7 * inch,
            ],
            sty,
        )
    )
    story.append(Spacer(1, 10))

    # Per-sample narrative
    story.append(Paragraph("4.3 Per-ID explanations", sty["h2"]))
    for _, r in mismatches.iterrows():
        sid = clean(r.get("article_sample_id")) or clean(r.get("Original_ID"))
        block = [
            Paragraph(
                f"{sid} ({clean(r.get('Original_ID'))} / {clean(r.get('Run'))})",
                sty["h2"],
            ),
            Paragraph(
                f"<b>File:</b> {clean(r.get('File_Name'))}<br/>"
                f"<b>Role:</b> {clean(r.get('record_role'))}<br/>"
                f"<b>Declared metadata:</b> {clean(r.get('declared_metadata'))}<br/>"
                f"<b>Observed STAT respiratory signal:</b> "
                f"{clean(r.get('observed_sequence_signal')) or 'none ≥3'}<br/>"
                f"<b>Severity:</b> {clean(r.get('mismatch_severity'))}<br/>"
                f"<b>Explanation:</b> {clean(r.get('notes'))}",
                sty["body"],
            ),
        ]
        story.append(KeepTogether(block))

    story.append(Paragraph("4.4 What is not counted as a metadata↔sequence mismatch", sty["h2"]))
    story.append(
        ListFlowable(
            [
                ListItem(
                    Paragraph(
                        "Article Supplementary Tables 6–7 mNGS vs RVP disagreements "
                        "(contextual; many lack joinable Accuracy SRR IDs).",
                        sty["bullet"],
                    )
                ),
                ListItem(
                    Paragraph(
                        "79 unresolved NP/NUC↔P-NP identifiers: sequenced and STAT-profiled, "
                        "but no clinical organism label to compare.",
                        sty["bullet"],
                    )
                ),
                ListItem(
                    Paragraph(
                        "STAT fetch failures SRR28242333 (NP363) and SRR28241792 (NP270): "
                        "QC/incomplete, not mislabel calls.",
                        sty["bullet"],
                    )
                ),
                ListItem(
                    Paragraph(
                        "Off-panel virus detections and background taxa (anelloviruses/phages).",
                        sty["bullet"],
                    )
                ),
            ],
            bulletType="bullet",
            leftIndent=12,
        )
    )

    story.append(Paragraph("4.5 Bottom-line ID list", sty["h2"]))
    strong = mismatches[mismatches["mismatch_severity"] == "strong_mismatch"]
    borderline = mismatches[mismatches["mismatch_severity"] == "borderline_mismatch"]
    near = mismatches[mismatches["mismatch_severity"] == "low_confidence_near_miss"]
    story.append(
        Paragraph(
            "<b>Strong metadata↔sequence mismatch:</b> "
            + (
                ", ".join(
                    clean(r.get("article_sample_id")) or clean(r.get("Original_ID"))
                    for _, r in strong.iterrows()
                )
                or "none"
            )
            + "<br/>"
            "<b>Borderline mismatch:</b> "
            + (
                ", ".join(
                    clean(r.get("article_sample_id")) or clean(r.get("Original_ID"))
                    for _, r in borderline.iterrows()
                )
                or "none"
            )
            + "<br/>"
            "<b>Low-confidence near-miss (expected pathogen &lt;3 STAT counts):</b> "
            + (
                ", ".join(
                    clean(r.get("article_sample_id")) or clean(r.get("Original_ID"))
                    for _, r in near.iterrows()
                )
                or "none"
            ),
            sty["body"],
        )
    )
    story.append(
        Paragraph(
            "Machine-readable companion: results/metadata_sequence_mismatches.csv. "
            "Full per-run audit: results/accuracy_metadata_audit.csv.",
            sty["small"],
        )
    )

    doc.build(story)
    return OUT_PDF


if __name__ == "__main__":
    path = build()
    print(f"Wrote {path}")
    print(f"Wrote {OUT_MISMATCH_CSV}")
