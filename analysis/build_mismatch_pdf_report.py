#!/usr/bin/env python3
"""Build a readable PDF of Accuracy samples whose metadata conflict with sequence data.

Keeps the required sections (intro, challenge, software/data, results) but leads with
exact sample-ID discrepancies: what metadata says vs what STAT sequence content shows.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    HRFlowable,
    KeepTogether,
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
ARTIFACT_PDF = Path("/opt/cursor/artifacts/metadata_consistency_report.pdf")


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
            fontSize=18,
            leading=22,
            spaceAfter=4,
            alignment=TA_CENTER,
        ),
        "subtitle": ParagraphStyle(
            "SubtitleCustom",
            parent=base["Normal"],
            fontSize=10,
            leading=13,
            textColor=colors.HexColor("#555555"),
            alignment=TA_CENTER,
            spaceAfter=14,
        ),
        "h1": ParagraphStyle(
            "H1Custom",
            parent=base["Heading1"],
            fontSize=13,
            leading=16,
            spaceBefore=12,
            spaceAfter=6,
        ),
        "h2": ParagraphStyle(
            "H2Custom",
            parent=base["Heading2"],
            fontSize=11,
            leading=14,
            spaceBefore=8,
            spaceAfter=4,
            textColor=colors.HexColor("#111827"),
        ),
        "body": ParagraphStyle(
            "BodyCustom",
            parent=base["Normal"],
            fontSize=10,
            leading=14,
            spaceAfter=6,
        ),
        "id": ParagraphStyle(
            "IdCustom",
            parent=base["Normal"],
            fontSize=12,
            leading=15,
            spaceBefore=2,
            spaceAfter=4,
            textColor=colors.HexColor("#0f172a"),
        ),
        "label": ParagraphStyle(
            "LabelCustom",
            parent=base["Normal"],
            fontSize=9,
            leading=12,
            textColor=colors.HexColor("#334155"),
            spaceAfter=2,
        ),
        "tablecell": ParagraphStyle(
            "TableCell",
            parent=base["Normal"],
            fontSize=8.5,
            leading=11,
            alignment=TA_LEFT,
        ),
        "tableheader": ParagraphStyle(
            "TableHeader",
            parent=base["Normal"],
            fontSize=8.5,
            leading=11,
            textColor=colors.white,
        ),
        "small": ParagraphStyle(
            "SmallCustom",
            parent=base["Normal"],
            fontSize=8.5,
            leading=11,
            textColor=colors.HexColor("#475569"),
            spaceAfter=4,
        ),
    }


def discrepancy_text(row: pd.Series) -> tuple[str, str, str, str]:
    """Return (sample_id, metadata_says, sequence_shows, plain_english)."""
    sid = clean(row.get("article_sample_id")) or clean(row.get("Original_ID"))
    role = clean(row.get("record_role"))
    severity = clean(row.get("mismatch_severity"))

    if role == "negative_control":
        observed = clean(row.get("stat_respiratory_organism_keys"))
        organism = {
            "sars_cov_2": "SARS-CoV-2",
            "adenovirus": "Adenovirus (Mastadenovirus)",
        }.get(observed, observed or "a respiratory pathogen")
        count_note = ""
        notes = clean(row.get("notes"))
        if "SARS" in notes or "sars" in notes.lower():
            count_note = " — 20 STAT counts"
            organism = "SARS-CoV-2"
        elif "adenovirus" in notes.lower() or "Mastadenovirus" in notes:
            count_note = " — 3 STAT counts"
            organism = "Adenovirus (Mastadenovirus)"
        metadata_says = "Negative control (no respiratory pathogen expected)"
        sequence_shows = f"{organism} detected{count_note}"
        plain = (
            f"Metadata labels this run as a negative control, so it should not contain "
            f"respiratory pathogen signal. The deposited sequence instead shows "
            f"{organism}{count_note}."
        )
        return sid, metadata_says, sequence_shows, plain

    # Clinical samples with missing expected organism evidence
    pcr = clean(row.get("pcr_results"))
    mngs = clean(row.get("mngs_results"))
    metadata_says = mngs or pcr or clean(row.get("declared_metadata"))
    if pcr and mngs and pcr != mngs:
        metadata_says = f"PCR: {pcr}; mNGS metadata: {mngs}"

    notes = clean(row.get("notes"))
    if "Rhinovirus A:2/2" in notes:
        sequence_shows = "Rhinovirus A present but only 2 STAT counts (below ≥3 threshold)"
        plain = (
            f"Metadata says this sample is Rhinovirus-positive ({metadata_says}). "
            "STAT finds Rhinovirus A, but only at 2 counts — below the ≥3 evidence cutoff "
            "used for a confident call. So the label and the sequence do not fully agree "
            "at the required evidence level."
        )
    elif "Orthopneumovirus:1/1" in notes:
        sequence_shows = "RSV/Orthopneumovirus present but only 1 STAT count (below ≥3 threshold)"
        plain = (
            f"Metadata says this sample is RSV-positive ({metadata_says}). "
            "STAT finds Orthopneumovirus (RSV lineage) at only 1 count — below the ≥3 "
            "evidence cutoff. The declared RSV call is not supported at threshold."
        )
    else:
        sequence_shows = clean(row.get("observed_sequence_signal")) or "no qualifying pathogen signal"
        plain = (
            f"Metadata says: {metadata_says}. Sequence evidence shows: {sequence_shows}."
        )
    return sid, metadata_says, sequence_shows, plain


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
            if "sars_cov_2" in notes.lower() or ":20" in notes:
                return "strong_mismatch"
            return "borderline_mismatch"
        if "Below-threshold near-miss" in notes:
            return "low_confidence_near_miss"
        return "possible_mismatch"

    hard["mismatch_severity"] = hard.apply(severity, axis=1)
    hard["declared_metadata"] = hard.apply(
        lambda r: clean(r.get("mngs_results"))
        or clean(r.get("pcr_results"))
        or (
            "negative_control"
            if r["record_role"] == "negative_control"
            else clean(r.get("record_role"))
        ),
        axis=1,
    )
    hard["observed_sequence_signal"] = hard.apply(
        lambda r: clean(r.get("stat_respiratory_organism_keys"))
        or (
            "none >=3 STAT counts"
            if "lacking STAT" in str(r.get("notes") or "")
            else ""
        ),
        axis=1,
    )

    # Attach plain-language discrepancy fields
    sid_l, meta_l, seq_l, plain_l = [], [], [], []
    for _, r in hard.iterrows():
        sid, meta, seq, plain = discrepancy_text(r)
        sid_l.append(sid)
        meta_l.append(meta)
        seq_l.append(seq)
        plain_l.append(plain)
    hard["sample_id"] = sid_l
    hard["metadata_says"] = meta_l
    hard["sequence_shows"] = seq_l
    hard["discrepancy"] = plain_l

    severity_order = {
        "strong_mismatch": 0,
        "borderline_mismatch": 1,
        "low_confidence_near_miss": 2,
        "possible_mismatch": 3,
    }
    hard["_ord"] = hard["mismatch_severity"].map(severity_order).fillna(9)
    hard = hard.sort_values(["_ord", "sample_id", "Run"]).drop(columns=["_ord"])
    return hard


def make_table(data: list[list[str]], col_widths: list[float], sty) -> Table:
    wrapped = []
    for i, row in enumerate(data):
        style = sty["tableheader"] if i == 0 else sty["tablecell"]
        wrapped.append([Paragraph(str(c), style) for c in row])
    t = Table(wrapped, colWidths=col_widths, repeatRows=1)
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#111827")),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#cbd5e1")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                (
                    "ROWBACKGROUNDS",
                    (0, 1),
                    (-1, -1),
                    [colors.white, colors.HexColor("#f8fafc")],
                ),
            ]
        )
    )
    return t


def card(row: pd.Series, sty) -> KeepTogether:
    severity = clean(row["mismatch_severity"]).replace("_", " ")
    return KeepTogether(
        [
            Paragraph(
                f"<b>{clean(row['sample_id'])}</b> &nbsp;&nbsp;"
                f"<font size='9' color='#64748b'>{severity}</font>",
                sty["id"],
            ),
            Paragraph(
                f"<b>Run / file:</b> {clean(row['Run'])} &nbsp;·&nbsp; "
                f"{clean(row['File_Name'])}<br/>"
                f"<b>BioSample:</b> {clean(row.get('BioSample'))} &nbsp;·&nbsp; "
                f"<b>Original_ID:</b> {clean(row['Original_ID'])}",
                sty["label"],
            ),
            Paragraph(f"<b>Metadata says:</b> {clean(row['metadata_says'])}", sty["body"]),
            Paragraph(f"<b>Sequence shows:</b> {clean(row['sequence_shows'])}", sty["body"]),
            Paragraph(f"<b>Discrepancy:</b> {clean(row['discrepancy'])}", sty["body"]),
            HRFlowable(
                width="100%",
                thickness=0.6,
                color=colors.HexColor("#e2e8f0"),
                spaceBefore=2,
                spaceAfter=8,
            ),
        ]
    )


def build() -> Path:
    audit = pd.read_csv(AUDIT_CSV)
    mismatches = mismatch_frame(audit)

    export_cols = [
        "sample_id",
        "Original_ID",
        "article_sample_id",
        "Run",
        "BioSample",
        "File_Name",
        "record_role",
        "metadata_says",
        "sequence_shows",
        "discrepancy",
        "mismatch_severity",
        "pcr_results",
        "mngs_results",
        "stat_respiratory_organism_keys",
        "notes",
    ]
    mismatches[[c for c in export_cols if c in mismatches.columns]].to_csv(
        OUT_MISMATCH_CSV, index=False
    )

    sty = styles()
    doc = SimpleDocTemplate(
        str(OUT_PDF),
        pagesize=letter,
        leftMargin=0.7 * inch,
        rightMargin=0.7 * inch,
        topMargin=0.65 * inch,
        bottomMargin=0.65 * inch,
        title="Metadata vs Sequence Discrepancies",
        author="c121_metagenomics audit",
    )

    story: list = []
    story.append(Paragraph("Metadata vs Sequence Discrepancies", sty["title"]))
    story.append(
        Paragraph(
            f"Accuracy cohort &nbsp;·&nbsp; {len(mismatches)} mismatched records "
            f"&nbsp;·&nbsp; {utc_now()}",
            sty["subtitle"],
        )
    )

    # ---- Intro (short) ----
    story.append(Paragraph("1. Introduction", sty["h1"]))
    story.append(
        Paragraph(
            "This report lists Accuracy-cohort sequencing runs where the deposited "
            "sample metadata do not match the organisms detected in the actual "
            "sequence data (NCBI STAT taxonomy, ≥3-count evidence cutoff).",
            sty["body"],
        )
    )

    # ---- Challenge (short, practical) ----
    story.append(Paragraph("2. Challenge", sty["h1"]))
    story.append(
        Paragraph(
            "For each run we compare (a) the declared label — negative control, or "
            "clinical PCR/mNGS organism call — with (b) respiratory pathogens found in "
            "the public STAT taxonomy for that SRR. A mismatch means the label and the "
            "sequence content disagree.",
            sty["body"],
        )
    )

    # ---- Software / data (compact) ----
    story.append(Paragraph("3. Software and data", sty["h1"]))
    story.append(
        Paragraph(
            "<b>Software:</b> Python audit scripts in <font face='Courier'>analysis/</font> "
            "(manifest build, STAT sequence audit, this PDF).<br/>"
            "<b>Data:</b> Zenodo Accuracy master sheet, Supplementary Data 1 clinical "
            "labels, SRA/BioProject PRJNA1084017 runs, NCBI STAT <font face='Courier'>"
            "run_taxonomy</font>.<br/>"
            "<b>Rule:</b> organism call requires STAT count ≥3 (proxy for ≥3 "
            "non-overlapping reads/contigs).",
            sty["body"],
        )
    )

    # ---- Results: ID list first ----
    story.append(Paragraph("4. Results — mismatched sample IDs", sty["h1"]))
    story.append(
        Paragraph(
            "<b>Sample IDs:</b> "
            + ", ".join(f"<b>{clean(s)}</b>" for s in mismatches["sample_id"].tolist()),
            sty["body"],
        )
    )

    summary_rows = [
        ["Sample ID", "Metadata says", "Sequence shows", "Severity"],
    ]
    for _, r in mismatches.iterrows():
        summary_rows.append(
            [
                clean(r["sample_id"]),
                clean(r["metadata_says"]),
                clean(r["sequence_shows"]),
                clean(r["mismatch_severity"]).replace("_", " "),
            ]
        )
    story.append(
        make_table(
            summary_rows,
            [0.9 * inch, 2.15 * inch, 2.55 * inch, 1.2 * inch],
            sty,
        )
    )
    story.append(Spacer(1, 10))

    story.append(Paragraph("Exact discrepancy per sample", sty["h2"]))
    for _, r in mismatches.iterrows():
        story.append(card(r, sty))

    story.append(Paragraph("Not counted as mismatches", sty["h2"]))
    story.append(
        Paragraph(
            "• STAT missing for <b>NP363</b> (<font face='Courier'>SRR28242333</font>) and "
            "<b>NP270</b> (<font face='Courier'>SRR28241792</font>) — incomplete data, not "
            "a label conflict.<br/>"
            "• 79 NP/NUC IDs with no published clinical crosswalk — no metadata organism "
            "to compare.<br/>"
            "• Article mNGS-vs-RVP tables are outside this sequence-vs-metadata list.",
            sty["small"],
        )
    )
    story.append(
        Paragraph(
            "CSV: results/metadata_sequence_mismatches.csv &nbsp;·&nbsp; "
            "Full audit: results/accuracy_metadata_audit.csv",
            sty["small"],
        )
    )

    doc.build(story)
    ARTIFACT_PDF.parent.mkdir(parents=True, exist_ok=True)
    ARTIFACT_PDF.write_bytes(OUT_PDF.read_bytes())
    return OUT_PDF


if __name__ == "__main__":
    path = build()
    print(f"Wrote {path}")
    print(f"Wrote {OUT_MISMATCH_CSV}")
    print(f"Wrote {ARTIFACT_PDF}")
