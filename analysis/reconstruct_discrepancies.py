#!/usr/bin/env python3
"""Reconstruct article-reported mNGS/RVP discrepancies (Supplementary Tables 6-7).

Scope: parse Supplementary Tables 6 and 7, join confident Supplementary Data 1
sample IDs, and emit sample- and organism-level discrepancy artifacts.

Does not download FASTQs, classify reads, or write the final audit report.
Does not invent identifier crosswalks when published files do not support them.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pymupdf
import pandas as pd


ARTICLE_URL = "https://www.nature.com/articles/s41467-024-51470-y"
SUPP_INFO_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2F"
    "s41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM1_ESM.pdf"
)
SUPP_DATA_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2F"
    "s41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM4_ESM.xlsx"
)

EXPECTED_TABLE6_SAMPLES = 14
EXPECTED_TABLE7_SAMPLES = 9
EXPECTED_MAINTEXT_MNGS_NEG_RVP_POS = 7
EXPECTED_PRE_DTCA_SENSITIVITY = (103, 110)
EXPECTED_PRE_DTCA_SPECIFICITY = (76, 81)
EXPECTED_PRE_DTCA_ACCURACY = (179, 191)
EXPECTED_POST_DTCA_PPA = (110.5, 113)
EXPECTED_POST_DTCA_NPA = (76.5, 78)
EXPECTED_POST_DTCA_OPA = (187, 191)


# Curated from Supplementary PDF word coordinates (pages 6-7). Each dict is one
# Table 6/7 sample row. Sample-ID joins are evidence-based only.
TABLE6_ROWS: list[dict[str, Any]] = [
    {
        "table": "6",
        "table_row": 1,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "RSV",
        "mngs_organisms": "RSV; Rhinovirus B",
        "discrepancy_testing": "(+) Rhinovirus",
        "additional_data": "",
        "final_mngs_by_dtca": "TP RSV; TP Rhinovirus",
        "article_sample_id": "P-NP-1",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: RSV+ by PCR; RSV A+ and Rhinovirus B+ by mNGS; "
            "orthogonal Rhinovirus Ct=22.3; comment Flu A/B/RSV PCR only"
        ),
    },
    {
        "table": "6",
        "table_row": 2,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "RSV",
        "mngs_organisms": "RSV; EV-D68",
        "discrepancy_testing": "(+) Enterovirus",
        "additional_data": "",
        "final_mngs_by_dtca": "TP RSV; TP Enterovirus",
        "article_sample_id": "P-NP-6",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: RSV+ by PCR; RSV B+ and EVD68+ by mNGS; "
            "orthogonal Enterovirus Ct=36.7; comment Flu A/B/RSV PCR only"
        ),
    },
    {
        "table": "6",
        "table_row": 3,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "RSV",
        "mngs_organisms": "RSV; Adenovirus C",
        "discrepancy_testing": "(+) Adenovirus",
        "additional_data": "",
        "final_mngs_by_dtca": "TP RSV; TP Adenovirus",
        "article_sample_id": "P-NP-36",
        "id_mapping_status": "conflict",
        "id_mapping_evidence": (
            "Only Supp Data 1 RSV+/Adenovirus-C candidate with Flu A/B/RSV-only "
            "comment, but Supp Data 1 lists mNGS 'RSV not detected' while Table 6 "
            "and DTCA require mNGS RSV detection (TP RSV). Mapping retained as "
            "conflicted candidate only."
        ),
    },
    {
        "table": "6",
        "table_row": 4,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "Adenovirus; Metapneumovirus",
        "mngs_organisms": "Adenovirus C; Metapneumovirus; Bocavirus",
        "discrepancy_testing": "Bocavirus ND",
        "additional_data": (
            "Bocavirus mNGS reads are specific to the viral genus. Patient with "
            "clinical syndrome consistent with acute viral infection."
        ),
        "final_mngs_by_dtca": "TP Adenovirus; TP Metapneumovirus; TP Bocavirus",
        "article_sample_id": "P-NP-113",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: Adenovirus+/MPV+ by PCR; Mastadenovirus C, "
            "Metapneumovirus, and Bocavirus by mNGS; Bocavirus orthogonal "
            "Not Detected"
        ),
    },
    {
        "table": "6",
        "table_row": 5,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Coronavirus OC43",
        "discrepancy_testing": "Coronavirus OC43 ND",
        "additional_data": (
            "Detected by mNGS x2 with high RPM and high coverage. No evidence of "
            "cross-contamination. Patient with clinical syndrome consistent with "
            "acute respiratory viral infection."
        ),
        "final_mngs_by_dtca": "TP Coronavirus OC43",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; Supp Data 1 enumerates virus-positive samples and "
            "does not identify this negative-cohort specimen."
        ),
    },
    {
        "table": "6",
        "table_row": 6,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Coronavirus OC43",
        "discrepancy_testing": "(+) Coronavirus OC43",
        "additional_data": "",
        "final_mngs_by_dtca": "TP Coronavirus OC43",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; no published sample-ID crosswalk in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 7,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Rhinovirus A",
        "discrepancy_testing": "(+) Rhinovirus",
        "additional_data": "",
        "final_mngs_by_dtca": "TP Rhinovirus",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; no published sample-ID crosswalk in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 8,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Rhinovirus A",
        "discrepancy_testing": "(+) Rhinovirus",
        "additional_data": "",
        "final_mngs_by_dtca": "TP Rhinovirus",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; no published sample-ID crosswalk in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 9,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Coronavirus HKU1",
        "discrepancy_testing": "Coronavirus HKU1 ND",
        "additional_data": "",
        "final_mngs_by_dtca": "FP Coronavirus HKU1",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; no published sample-ID crosswalk in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 10,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Coronavirus HKU1",
        "discrepancy_testing": "(+) Coronavirus HKU1",
        "additional_data": "",
        "final_mngs_by_dtca": "TP Coronavirus HKU1",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; no published sample-ID crosswalk in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 11,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Parainfluenza 1; SARS-CoV-2",
        "discrepancy_testing": "(+) Parainfluenza 1; SARS-CoV-2 ND",
        "additional_data": "",
        "final_mngs_by_dtca": "TP Parainfluenza 1; FP SARS-CoV-2",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; no published sample-ID crosswalk in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 12,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Coronavirus NL63",
        "discrepancy_testing": "(+) Coronavirus NL63",
        "additional_data": "",
        "final_mngs_by_dtca": "TP Coronavirus NL63",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; no published sample-ID crosswalk in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 13,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Rhinovirus C",
        "discrepancy_testing": "Rhinovirus ND",
        "additional_data": (
            "Mismatches in the forward primer and probe sequences for "
            "rhinovirus/enterovirus PCR. No evidence of contamination. Patient "
            "with clinical syndrome consistent with acute respiratory viral "
            "infection. Related main-text Fig. 7 (12 non-overlapping reads)."
        ),
        "final_mngs_by_dtca": "TP Rhinovirus",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative rhinovirus C primer-mismatch case described in text/"
            "Fig. 7 without an explicit sample ID in Supp Data 1."
        ),
    },
    {
        "table": "6",
        "table_row": 14,
        "discrepancy_class": "mNGS_positive_RVP_negative",
        "rvp_organisms": "ND",
        "mngs_organisms": "Coronavirus OC43",
        "discrepancy_testing": "NT",
        "additional_data": "",
        "final_mngs_by_dtca": "NA Coronavirus OC43",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "RVP-negative case; orthogonal testing not performed; no sample-ID "
            "crosswalk in Supp Data 1."
        ),
    },
]


TABLE7_ROWS: list[dict[str, Any]] = [
    {
        "table": "7",
        "table_row": 1,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Influenza A; Parainfluenza 3",
        "clinical_rvp_signal": "NA",
        "mngs_organisms": "Influenza A",
        "high_host_background": "Y",
        "discrepancy_testing": "(+) Influenza A; Parainfluenza 3 ND",
        "discrepancy_testing_platform": "MagPix",
        "final_mngs_by_dtca": "TP Influenza A; TN Parainfluenza 3",
        "article_sample_id": "P-NP-19",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: Flu A+ and Parainfluenza 3+ by PCR; mNGS Flu A+(H3N2) "
            "with Parainfluenza 3 not detected; comment notes 1 subthreshold "
            "parainfluenza read"
        ),
        "maintext7_note": (
            "Co-detection sample: TP half for Influenza A plus TN half for "
            "missed Parainfluenza 3; TN/FN halves across Table 7 sum to the "
            "main-text 7."
        ),
    },
    {
        "table": "7",
        "table_row": 2,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Influenza A; SARS-CoV-2",
        "clinical_rvp_signal": "SARS-CoV-2 CT: 36.4",
        "mngs_organisms": "Influenza A",
        "high_host_background": "Y",
        "discrepancy_testing": "Influenza A NT; SARS-CoV-2 ND",
        "discrepancy_testing_platform": "Verigene",
        "final_mngs_by_dtca": "TP Influenza A; TN SARS-CoV-2",
        "article_sample_id": "P-NP-33",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: Flu A+ and SARS-CoV-2+ by PCR; mNGS Flu A+ with "
            "SARS-CoV-2 not detected; Table 7 CT 36.4 / Verigene unique among "
            "Flu A+SARS pairs"
        ),
        "maintext7_note": "Co-detection miss of high-Ct SARS-CoV-2; TN after DTCA",
    },
    {
        "table": "7",
        "table_row": 3,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Influenza A; RSV",
        "clinical_rvp_signal": "NA",
        "mngs_organisms": "ND",
        "high_host_background": "Y",
        "discrepancy_testing": "Influenza A ND; (+) RSV B",
        "discrepancy_testing_platform": "GenMark",
        "final_mngs_by_dtca": "TN Influenza A; FN RSV",
        "article_sample_id": "",
        "id_mapping_status": "unresolved",
        "id_mapping_evidence": (
            "No Supp Data 1 sample with both Influenza A and RSV PCR-positive and "
            "mNGS ND. Not inventing a crosswalk."
        ),
        "maintext7_note": (
            "Paper's dual-adjudication example in one sample (0.5 TN + 0.5 FN)"
        ),
    },
    {
        "table": "7",
        "table_row": 4,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Influenza A; SARS-CoV-2",
        "clinical_rvp_signal": "NA",
        "mngs_organisms": "SARS-CoV-2",
        "high_host_background": "N",
        "discrepancy_testing": "Influenza A ND; SARS-CoV-2 NT",
        "discrepancy_testing_platform": "MagPix",
        "final_mngs_by_dtca": "TP SARS-CoV-2; TN Influenza A",
        "article_sample_id": "P-NP-64",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: Flu A (H3)+ and SARS-CoV-2+ by PCR; mNGS SARS-CoV-2+ "
            "only (Influenza A absent)"
        ),
        "maintext7_note": "Co-detection miss of Influenza A; TN after DTCA",
    },
    {
        "table": "7",
        "table_row": 5,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Rhinovirus/Enterovirus",
        "clinical_rvp_signal": "Rhinovirus/Enterovirus MFI: 52",
        "mngs_organisms": "ND",
        "high_host_background": "Y",
        "discrepancy_testing": "Rhinovirus/Enterovirus ND",
        "discrepancy_testing_platform": "MagPix",
        "final_mngs_by_dtca": "TN Rhinovirus/Enterovirus",
        "article_sample_id": "P-NP-99",
        "id_mapping_status": "probable",
        "id_mapping_evidence": (
            "Supp Data 1 Rhinovirus/Enterovirus+ with mNGS ND and comment "
            "'False Negative'; Table 7 has two such rows (MFI 52 vs NA). "
            "Assigned to P-NP-99/P-NP-100 as an unordered pair; MFI cannot "
            "distinguish them from published files alone."
        ),
        "maintext7_note": "Complete mNGS miss; TN after repeat RVP/DTCA",
    },
    {
        "table": "7",
        "table_row": 6,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Rhinovirus/Enterovirus",
        "clinical_rvp_signal": "NA",
        "mngs_organisms": "ND",
        "high_host_background": "N",
        "discrepancy_testing": "Rhinovirus/Enterovirus ND",
        "discrepancy_testing_platform": "MagPix",
        "final_mngs_by_dtca": "TN Rhinovirus/Enterovirus",
        "article_sample_id": "P-NP-100",
        "id_mapping_status": "probable",
        "id_mapping_evidence": (
            "Paired with Table 7 row 5 / P-NP-99; same Supp Data 1 pattern. "
            "Row-to-ID order is not uniquely determined."
        ),
        "maintext7_note": "Complete mNGS miss; TN after repeat RVP/DTCA",
    },
    {
        "table": "7",
        "table_row": 7,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Rhinovirus/Enterovirus; Parainfluenza 2",
        "clinical_rvp_signal": "NA",
        "mngs_organisms": "Rhinovirus A",
        "high_host_background": "N",
        "discrepancy_testing": "(+) Rhinovirus/Enterovirus; Parainfluenza 2 ND",
        "discrepancy_testing_platform": "MagPix",
        "final_mngs_by_dtca": "TP Rhinovirus/Enterovirus; TN Parainfluenza 2",
        "article_sample_id": "P-NP-101",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1 unique Rhinovirus/Enterovirus + Parainfluenza 2+ with "
            "mNGS Rhinovirus A only"
        ),
        "maintext7_note": "Co-detection miss of Parainfluenza 2; TN after DTCA",
    },
    {
        "table": "7",
        "table_row": 8,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "Metapneumovirus",
        "clinical_rvp_signal": "NA",
        "mngs_organisms": "ND",
        "high_host_background": "Y",
        "discrepancy_testing": "(+) Metapneumovirus",
        "discrepancy_testing_platform": "GenMark",
        "final_mngs_by_dtca": "FN Metapneumovirus",
        "article_sample_id": "P-NP-109",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: Metapneumovirus+; mNGS respiratory virus not detected; "
            "comment False Negative"
        ),
        "maintext7_note": "Confirmed mNGS false negative after DTCA",
    },
    {
        "table": "7",
        "table_row": 9,
        "discrepancy_class": "mNGS_negative_RVP_positive",
        "rvp_organisms": "RSV A",
        "clinical_rvp_signal": "NA",
        "mngs_organisms": "ND",
        "high_host_background": "Y",
        "discrepancy_testing": "RSV ND",
        "discrepancy_testing_platform": "MagPix",
        "final_mngs_by_dtca": "TN RSV A",
        "article_sample_id": "BAL105",
        "id_mapping_status": "confident",
        "id_mapping_evidence": (
            "Supp Data 1: RSV A+; mNGS respiratory virus not detected (CMV only); "
            "comment False Negative; unique BAL RSV A ND pattern"
        ),
        "maintext7_note": (
            "Table 7 n=9 sample rows; fractional TN+FN weights across all rows "
            "equal the main-text 7 (5.5 TN + 1.5 FN)."
        ),
    },
]


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return re.sub(r"\s+", " ", str(value)).strip()


def split_final_calls(final_text: str) -> list[tuple[str, str]]:
    """Parse 'TP RSV; FP SARS-CoV-2' into (classification, organism) pairs."""
    calls: list[tuple[str, str]] = []
    for part in re.split(r"\s*;\s*", final_text):
        part = part.strip()
        if not part:
            continue
        match = re.match(r"^(TP|FP|TN|FN|NA)\s+(.+)$", part)
        if not match:
            raise ValueError(f"Unparseable DTCA call: {part!r} in {final_text!r}")
        calls.append((match.group(1), match.group(2).strip()))
    return calls


def concise_explanation(row: dict[str, Any]) -> str:
    bits = [
        f"Table {row['table']} row {row['table_row']}: RVP[{row['rvp_organisms']}] "
        f"vs mNGS[{row['mngs_organisms']}]."
    ]
    testing = clean(row.get("discrepancy_testing"))
    if testing:
        bits.append(f"Orthogonal/repeat testing: {testing}.")
    platform = clean(row.get("discrepancy_testing_platform"))
    if platform:
        bits.append(f"Platform: {platform}.")
    additional = clean(row.get("additional_data"))
    if additional:
        bits.append(additional)
    bits.append(f"Final mNGS by DTCA: {row['final_mngs_by_dtca']}.")
    if row.get("id_mapping_status") != "confident":
        bits.append(f"ID mapping {row['id_mapping_status']}: {row['id_mapping_evidence']}")
    return " ".join(bits)


def chart_adjudication_summary(row: dict[str, Any]) -> str:
    additional = clean(row.get("additional_data"))
    finals = row["final_mngs_by_dtca"]
    if additional:
        return f"{additional} Final classification(s): {finals}."
    testing = clean(row.get("discrepancy_testing"))
    if "ND" in testing and any(x in finals for x in ("TN", "TP", "FN", "FP")):
        return (
            "DTCA used orthogonal/repeat testing and clinical review; "
            f"final classification(s): {finals}."
        )
    return f"Final classification(s) after DTCA: {finals}."


def load_clinical_lookup(clinical_samples_csv: Path) -> dict[str, dict[str, str]]:
    samples = pd.read_csv(clinical_samples_csv, dtype=str).fillna("")
    lookup: dict[str, dict[str, str]] = {}
    for _, sample in samples.iterrows():
        lookup[clean(sample["article_sample_id"])] = {
            "sample_type": clean(sample["sample_type"]),
            "supp_data1_pcr_results": clean(sample["pcr_results"]),
            "supp_data1_mngs_results": clean(sample["mngs_results"]),
            "supp_data1_viral_load_copies_ml": clean(sample["viral_load_copies_ml"]),
            "supp_data1_clinical_diagnosis": clean(sample["clinical_diagnosis"]),
            "supp_data1_disease_severity": clean(sample["disease_severity"]),
            "supp_data1_symptoms": clean(sample["symptoms"]),
            "supp_data1_discrepancy_testing_ct": clean(sample["discrepancy_testing_ct"]),
            "supp_data1_comments": clean(sample["comments"]),
        }
    return lookup


def attach_clinical(row: dict[str, Any], lookup: dict[str, dict[str, str]]) -> dict[str, Any]:
    sample_id = clean(row.get("article_sample_id"))
    clinical = lookup.get(sample_id, {})
    empty = {
        "sample_type": "",
        "supp_data1_pcr_results": "",
        "supp_data1_mngs_results": "",
        "supp_data1_viral_load_copies_ml": "",
        "supp_data1_clinical_diagnosis": "",
        "supp_data1_disease_severity": "",
        "supp_data1_symptoms": "",
        "supp_data1_discrepancy_testing_ct": "",
        "supp_data1_comments": "",
    }
    merged = {**row, **empty, **clinical}
    merged["chart_adjudication"] = chart_adjudication_summary(row)
    merged["explanation"] = concise_explanation(row)
    merged["initial_rvp_call"] = row["rvp_organisms"]
    merged["initial_mngs_call"] = row["mngs_organisms"]
    return merged


def expand_organism_rows(sample_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    organism_rows: list[dict[str, Any]] = []
    for sample in sample_rows:
        calls = split_final_calls(sample["final_mngs_by_dtca"])
        n_calls = len(calls)
        fractional_weight = 1.0 / n_calls if n_calls else 1.0
        for call_index, (classification, organism) in enumerate(calls, start=1):
            organism_rows.append(
                {
                    "discrepancy_uid": (
                        f"T{sample['table']}-R{sample['table_row']}-C{call_index}"
                    ),
                    "table": sample["table"],
                    "table_row": sample["table_row"],
                    "organism_call_index": call_index,
                    "discrepancy_class": sample["discrepancy_class"],
                    "article_sample_id": sample.get("article_sample_id", ""),
                    "id_mapping_status": sample.get("id_mapping_status", ""),
                    "organism": organism,
                    "final_classification": classification,
                    "fractional_sample_weight": fractional_weight,
                    "initial_rvp_call": sample["rvp_organisms"],
                    "initial_mngs_call": sample["mngs_organisms"],
                    "orthogonal_or_repeat_testing": sample.get("discrepancy_testing", ""),
                    "discrepancy_testing_platform": sample.get(
                        "discrepancy_testing_platform", ""
                    ),
                    "clinical_rvp_signal": sample.get("clinical_rvp_signal", ""),
                    "high_host_background": sample.get("high_host_background", ""),
                    "additional_data": sample.get("additional_data", ""),
                    "chart_adjudication": sample.get("chart_adjudication", ""),
                    "explanation": sample.get("explanation", ""),
                }
            )
    return organism_rows


def verify_pdf_titles(supp_pdf: Path) -> dict[str, Any]:
    doc = pymupdf.open(supp_pdf)
    page6 = doc[5].get_text("text")
    page7 = doc[6].get_text("text")
    checks = {
        "table6_title_has_n14": bool(
            re.search(r"Supplementary Table 6[\s\S]*?\(n\s*=\s*14 Samples", page6)
        ),
        "table7_title_has_n9": bool(
            re.search(r"Supplementary Table 7[\s\S]*?\(n\s*=\s*9 Samples", page7)
        ),
        "table6_mentions_tp_rhinovirus": "TP" in page6 and "Rhinovirus" in page6,
        "table7_mentions_fn_metapneumovirus": "FN" in page7
        and "Metapneumovirus" in page7,
        "pdf_sha256": digest(supp_pdf),
        "pdf_pages": doc.page_count,
    }
    if not checks["table6_title_has_n14"]:
        raise ValueError("Could not verify Supplementary Table 6 n=14 in PDF")
    if not checks["table7_title_has_n9"]:
        raise ValueError("Could not verify Supplementary Table 7 n=9 in PDF")
    return checks


REQUIRED_SAMPLE_FIELDS = (
    "discrepancy_class",
    "rvp_organisms",
    "mngs_organisms",
    "discrepancy_testing",
    "final_mngs_by_dtca",
    "chart_adjudication",
    "explanation",
    "id_mapping_status",
    "id_mapping_evidence",
)


def verify_sample_field_completeness(sample_df: pd.DataFrame) -> dict[str, Any]:
    """Ensure every Table 6/7 sample has adjudication fields filled."""
    missing: list[dict[str, Any]] = []
    for _, row in sample_df.iterrows():
        for field in REQUIRED_SAMPLE_FIELDS:
            if not clean(row.get(field)):
                missing.append(
                    {
                        "table": clean(row.get("table")),
                        "table_row": int(row.get("table_row")),
                        "field": field,
                    }
                )
        # Final classification must parse into organism calls.
        try:
            calls = split_final_calls(clean(row.get("final_mngs_by_dtca")))
        except ValueError as exc:
            missing.append(
                {
                    "table": clean(row.get("table")),
                    "table_row": int(row.get("table_row")),
                    "field": "final_mngs_by_dtca",
                    "error": str(exc),
                }
            )
            continue
        if not calls:
            missing.append(
                {
                    "table": clean(row.get("table")),
                    "table_row": int(row.get("table_row")),
                    "field": "final_mngs_by_dtca",
                    "error": "no organism classifications parsed",
                }
            )
    if missing:
        raise ValueError(f"Incomplete sample-level adjudication fields: {missing}")
    return {
        "required_fields": list(REQUIRED_SAMPLE_FIELDS),
        "samples_checked": int(len(sample_df)),
        "all_required_fields_present": True,
    }


def verify_counts(sample_df: pd.DataFrame, organism_df: pd.DataFrame) -> dict[str, Any]:
    table6 = sample_df[sample_df["table"] == "6"]
    table7 = sample_df[sample_df["table"] == "7"]
    if len(table6) != EXPECTED_TABLE6_SAMPLES:
        raise ValueError(f"Expected {EXPECTED_TABLE6_SAMPLES} Table 6 rows, got {len(table6)}")
    if len(table7) != EXPECTED_TABLE7_SAMPLES:
        raise ValueError(f"Expected {EXPECTED_TABLE7_SAMPLES} Table 7 rows, got {len(table7)}")

    t7_org = organism_df[organism_df["table"] == "7"]
    tn_weight = float(
        t7_org.loc[t7_org["final_classification"] == "TN", "fractional_sample_weight"].sum()
    )
    fn_weight = float(
        t7_org.loc[t7_org["final_classification"] == "FN", "fractional_sample_weight"].sum()
    )
    tp_weight = float(
        t7_org.loc[t7_org["final_classification"] == "TP", "fractional_sample_weight"].sum()
    )
    # Main-text "7" reconciles as Table 7 TN+FN fractional weights (5.5+1.5),
    # while the table itself lists 9 sample rows (co-detections also carry TP halves).
    maintext7_weight = tn_weight + fn_weight
    completeness = verify_sample_field_completeness(sample_df)

    return {
        "table6_sample_rows": int(len(table6)),
        "table6_mngs_pos_rvp_neg_count": int(len(table6)),
        "table7_sample_rows": int(len(table7)),
        "maintext_mngs_neg_rvp_pos_count": EXPECTED_MAINTEXT_MNGS_NEG_RVP_POS,
        "expected_maintext_mngs_neg_rvp_pos": EXPECTED_MAINTEXT_MNGS_NEG_RVP_POS,
        "table7_TN_weight": tn_weight,
        "table7_FN_weight": fn_weight,
        "table7_TP_weight": tp_weight,
        "table7_TN_plus_FN_weight": maintext7_weight,
        "paper_maintext_TN_reclassified": 5.5,
        "paper_maintext_FN_remaining_implied": 1.5,
        "table7_TN_matches_paper_5_5": bool(abs(tn_weight - 5.5) < 1e-9),
        "table7_TN_plus_FN_matches_maintext_7": bool(
            abs(maintext7_weight - EXPECTED_MAINTEXT_MNGS_NEG_RVP_POS) < 1e-9
        ),
        "exact_counts": {
            "mngs_positive_rvp_negative_samples": EXPECTED_TABLE6_SAMPLES,
            "mngs_negative_rvp_positive_maintext_equivalents": EXPECTED_MAINTEXT_MNGS_NEG_RVP_POS,
            "table7_published_sample_rows": EXPECTED_TABLE7_SAMPLES,
            "reconciliation": (
                "Main-text 7 = Table 7 TN fractional weight 5.5 + FN fractional "
                "weight 1.5; Table 7 lists 9 sample rows because co-detections "
                "also contribute TP halves that are not part of the 7."
            ),
        },
        "sample_field_completeness": completeness,
        "pre_dtca_sensitivity": {
            "reported": f"{EXPECTED_PRE_DTCA_SENSITIVITY[0]} of {EXPECTED_PRE_DTCA_SENSITIVITY[1]}",
            "implied_fn_samples": EXPECTED_PRE_DTCA_SENSITIVITY[1]
            - EXPECTED_PRE_DTCA_SENSITIVITY[0],
        },
        "pre_dtca_specificity": {
            "reported": f"{EXPECTED_PRE_DTCA_SPECIFICITY[0]} of {EXPECTED_PRE_DTCA_SPECIFICITY[1]}",
            "implied_fp_samples": EXPECTED_PRE_DTCA_SPECIFICITY[1]
            - EXPECTED_PRE_DTCA_SPECIFICITY[0],
        },
        "pre_dtca_accuracy": {
            "reported": f"{EXPECTED_PRE_DTCA_ACCURACY[0]} of {EXPECTED_PRE_DTCA_ACCURACY[1]}"
        },
        "post_dtca": {
            "ppa": f"{EXPECTED_POST_DTCA_PPA[0]} of {EXPECTED_POST_DTCA_PPA[1]} = "
            f"{EXPECTED_POST_DTCA_PPA[0] / EXPECTED_POST_DTCA_PPA[1]:.3%}",
            "npa": f"{EXPECTED_POST_DTCA_NPA[0]} of {EXPECTED_POST_DTCA_NPA[1]} = "
            f"{EXPECTED_POST_DTCA_NPA[0] / EXPECTED_POST_DTCA_NPA[1]:.3%}",
            "opa": f"{EXPECTED_POST_DTCA_OPA[0]} of {EXPECTED_POST_DTCA_OPA[1]} = "
            f"{EXPECTED_POST_DTCA_OPA[0] / EXPECTED_POST_DTCA_OPA[1]:.3%}",
        },
        "id_mapping_status_counts": {
            str(key): int(value)
            for key, value in sample_df["id_mapping_status"]
            .value_counts(dropna=False)
            .items()
        },
        "confident_or_probable_mapped_sample_ids": sorted(
            {
                clean(x)
                for x in sample_df.loc[
                    sample_df["id_mapping_status"].isin(["confident", "probable"]),
                    "article_sample_id",
                ]
                if clean(x)
            }
        ),
    }


def build_unresolved_id_ledger(sample_df: pd.DataFrame) -> pd.DataFrame:
    """Machine-readable ledger of non-confident sample-ID joins."""
    rows: list[dict[str, Any]] = []
    for _, row in sample_df.iterrows():
        status = clean(row.get("id_mapping_status"))
        if status == "confident":
            continue
        rows.append(
            {
                "table": clean(row.get("table")),
                "table_row": int(row.get("table_row")),
                "discrepancy_class": clean(row.get("discrepancy_class")),
                "article_sample_id": clean(row.get("article_sample_id")),
                "id_mapping_status": status,
                "id_mapping_evidence": clean(row.get("id_mapping_evidence")),
                "rvp_organisms": clean(row.get("rvp_organisms")),
                "mngs_organisms": clean(row.get("mngs_organisms")),
                "final_mngs_by_dtca": clean(row.get("final_mngs_by_dtca")),
                "np_nuc_crosswalk_note": (
                    "No published NP/NUC↔P-NP crosswalk; do not invent lab-ID "
                    "aliases for this discrepancy row."
                ),
            }
        )
    return pd.DataFrame(rows)


def build_maintext7_equivalents(organism_df: pd.DataFrame) -> pd.DataFrame:
    """Organism calls whose fractional weights sum to the main-text 7."""
    t7 = organism_df[organism_df["table"] == "7"].copy()
    equivalents = t7[t7["final_classification"].isin(["TN", "FN"])].copy()
    equivalents.insert(0, "contributes_to_maintext_7", True)
    equivalents["maintext7_component"] = equivalents["final_classification"].map(
        {"TN": "reclassified_true_negative", "FN": "residual_false_negative"}
    )
    return equivalents.reset_index(drop=True)


def build_case_cards(sample_df: pd.DataFrame, organism_df: pd.DataFrame) -> list[dict[str, Any]]:
    """Compact per-sample cards for reporting agents."""
    cards: list[dict[str, Any]] = []
    for _, sample in sample_df.iterrows():
        table = clean(sample.get("table"))
        table_row = int(sample.get("table_row"))
        org_rows = organism_df[
            (organism_df["table"] == table) & (organism_df["table_row"] == table_row)
        ]
        organism_calls = [
            {
                "organism": clean(org.get("organism")),
                "final_classification": clean(org.get("final_classification")),
                "fractional_sample_weight": float(org.get("fractional_sample_weight")),
                "orthogonal_or_repeat_testing": clean(org.get("orthogonal_or_repeat_testing")),
            }
            for _, org in org_rows.iterrows()
        ]
        cards.append(
            {
                "case_uid": f"T{table}-R{table_row}",
                "table": table,
                "table_row": table_row,
                "discrepancy_class": clean(sample.get("discrepancy_class")),
                "article_sample_id": clean(sample.get("article_sample_id")),
                "id_mapping_status": clean(sample.get("id_mapping_status")),
                "id_mapping_evidence": clean(sample.get("id_mapping_evidence")),
                "organisms_rvp": clean(sample.get("rvp_organisms")),
                "organisms_mngs": clean(sample.get("mngs_organisms")),
                "orthogonal_or_repeat_pcr": clean(sample.get("discrepancy_testing")),
                "discrepancy_testing_platform": clean(
                    sample.get("discrepancy_testing_platform")
                ),
                "chart_adjudication": clean(sample.get("chart_adjudication")),
                "final_classification": clean(sample.get("final_mngs_by_dtca")),
                "concise_explanation": clean(sample.get("explanation")),
                "organism_calls": organism_calls,
                "contributes_to_maintext_7_weight": float(
                    org_rows.loc[
                        org_rows["final_classification"].isin(["TN", "FN"]),
                        "fractional_sample_weight",
                    ].sum()
                )
                if table == "7"
                else 0.0,
            }
        )
    return cards


def relpath(path: Path, project_root: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def write_readme(path: Path, verification: dict[str, Any]) -> None:
    lines = [
        "# Article discrepancy reconstruction",
        "",
        "Scope: Supplementary Tables 6 and 7 plus Supplementary Data 1 joins only.",
        "No FASTQ classification and no final inconsistency report.",
        "",
        "## Sources",
        f"- Article: {ARTICLE_URL}",
        f"- Supplementary Information (Tables 6-7): {SUPP_INFO_URL}",
        f"- Supplementary Data 1: {SUPP_DATA_URL}",
        "",
        "## Exact counts",
        f"- mNGS+/RVP- samples (Table 6): "
        f"{verification['exact_counts']['mngs_positive_rvp_negative_samples']}",
        f"- mNGS-/RVP+ main-text equivalents: "
        f"{verification['exact_counts']['mngs_negative_rvp_positive_maintext_equivalents']}",
        f"- Table 7 published sample rows: "
        f"{verification['exact_counts']['table7_published_sample_rows']} "
        f"(title n=9; TN+FN weight "
        f"{verification['table7_TN_plus_FN_weight']} = "
        f"TN {verification['table7_TN_weight']} + FN {verification['table7_FN_weight']}; "
        f"match={verification['table7_TN_plus_FN_matches_maintext_7']})",
        f"- Pre-DTCA sensitivity misses: "
        f"{verification['pre_dtca_sensitivity']['implied_fn_samples']} = 110-103",
        f"- Post-DTCA OPA: {verification['post_dtca']['opa']}",
        "",
        "## Outputs",
        "- `supplementary_table6_discrepancies.csv`",
        "- `supplementary_table7_discrepancies.csv`",
        "- `article_discrepancies_sample_level.csv`",
        "- `article_discrepancies_organism_level.csv`",
        "- `article_discrepancies_maintext7_tn_fn.csv`",
        "- `article_discrepancies_unresolved_ids.csv`",
        "- `../results/intermediate/discrepancy_case_cards.json`",
        "- `../results/intermediate/discrepancy_reconstruction_summary.json`",
        "",
        "## Unresolved issues preserved",
        "- Most Table 6 RVP-negative rows lack sample IDs in Supplementary Data 1.",
        "- Table 6 row 3 vs P-NP-36 RSV call conflicts between Table 6 and Supp Data 1.",
        "- Table 7 title n=9 vs main-text n=7: reconciled as TN+FN fractional weights 5.5+1.5=7.",
        "- P-NP-99/P-NP-100 cannot be ordered against the two rhinovirus Table 7 rows.",
        "- Influenza A + RSV / mNGS ND Table 7 row has no Supp Data 1 match.",
        "- NP/NUC↔P-NP lab-ID crosswalks are unpublished; not invented here "
        "(see `unresolved_identifier_mappings.csv`).",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_artifacts(project_root: Path) -> dict[str, Any]:
    raw_dir = project_root / "data" / "raw"
    norm_dir = project_root / "data" / "normalized"
    inter_dir = project_root / "results" / "intermediate"
    norm_dir.mkdir(parents=True, exist_ok=True)
    inter_dir.mkdir(parents=True, exist_ok=True)

    supp_pdf = raw_dir / "41467_2024_51470_MOESM1_ESM.pdf"
    supp_xlsx = raw_dir / "41467_2024_51470_MOESM4_ESM.xlsx"
    clinical_samples = norm_dir / "clinical_metadata_samples.csv"
    if not clinical_samples.exists():
        raise FileNotFoundError(
            f"Missing {clinical_samples}; run analysis/build_accuracy_manifest.py first"
        )
    if not supp_pdf.exists() or not supp_xlsx.exists():
        raise FileNotFoundError("Missing supplementary PDF and/or Data 1 workbook")

    pdf_checks = verify_pdf_titles(supp_pdf)
    lookup = load_clinical_lookup(clinical_samples)

    sample_rows = [attach_clinical(row, lookup) for row in TABLE6_ROWS + TABLE7_ROWS]
    sample_df = pd.DataFrame(sample_rows)
    organism_df = pd.DataFrame(expand_organism_rows(sample_rows))
    unresolved_df = build_unresolved_id_ledger(sample_df)
    maintext7_df = build_maintext7_equivalents(organism_df)
    case_cards = build_case_cards(sample_df, organism_df)

    table6_df = sample_df[sample_df["table"] == "6"].copy()
    table7_df = sample_df[sample_df["table"] == "7"].copy()

    table6_path = norm_dir / "supplementary_table6_discrepancies.csv"
    table7_path = norm_dir / "supplementary_table7_discrepancies.csv"
    sample_path = norm_dir / "article_discrepancies_sample_level.csv"
    organism_path = norm_dir / "article_discrepancies_organism_level.csv"
    unresolved_path = norm_dir / "article_discrepancies_unresolved_ids.csv"
    maintext7_path = norm_dir / "article_discrepancies_maintext7_tn_fn.csv"
    readme_path = norm_dir / "README_discrepancy_reconstruction.md"
    case_cards_path = inter_dir / "discrepancy_case_cards.json"
    summary_path = inter_dir / "discrepancy_reconstruction_summary.json"

    table6_df.to_csv(table6_path, index=False)
    table7_df.to_csv(table7_path, index=False)
    sample_df.to_csv(sample_path, index=False)
    organism_df.to_csv(organism_path, index=False)
    unresolved_df.to_csv(unresolved_path, index=False)
    maintext7_df.to_csv(maintext7_path, index=False)
    case_cards_path.write_text(json.dumps(case_cards, indent=2), encoding="utf-8")

    verification = verify_counts(sample_df, organism_df)
    write_readme(readme_path, verification)

    maintext7_weight = float(maintext7_df["fractional_sample_weight"].sum())
    if abs(maintext7_weight - EXPECTED_MAINTEXT_MNGS_NEG_RVP_POS) >= 1e-9:
        raise ValueError(
            f"Main-text 7 equivalents sum to {maintext7_weight}, expected "
            f"{EXPECTED_MAINTEXT_MNGS_NEG_RVP_POS}"
        )

    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": (
            "reconstruct-discrepancies only; no FASTQ audit and no final report"
        ),
        "exact_counts": {
            "mngs_positive_rvp_negative": 14,
            "mngs_negative_rvp_positive_maintext": 7,
            "table7_sample_rows_published": 9,
            "maintext7_tn_fn_weight": maintext7_weight,
        },
        "sources": {
            "article": ARTICLE_URL,
            "supplementary_information": {
                "url": SUPP_INFO_URL,
                "local_path": relpath(supp_pdf, project_root),
                "sha256": digest(supp_pdf),
            },
            "supplementary_data_1": {
                "url": SUPP_DATA_URL,
                "local_path": relpath(supp_xlsx, project_root),
                "sha256": digest(supp_xlsx),
            },
            "clinical_metadata_samples": relpath(clinical_samples, project_root),
        },
        "pdf_verification": pdf_checks,
        "verification": verification,
        "outputs": {
            "table6": relpath(table6_path, project_root),
            "table7": relpath(table7_path, project_root),
            "sample_level": relpath(sample_path, project_root),
            "organism_level": relpath(organism_path, project_root),
            "maintext7_tn_fn": relpath(maintext7_path, project_root),
            "unresolved_ids": relpath(unresolved_path, project_root),
            "case_cards": relpath(case_cards_path, project_root),
            "readme": relpath(readme_path, project_root),
            "summary": relpath(summary_path, project_root),
        },
        "key_findings": [
            "Exact counts verified: 14 mNGS+/RVP- (Table 6) and 7 mNGS-/RVP+ main-text equivalents (Table 7 TN 5.5 + FN 1.5).",
            "Supplementary Table 7 publishes 9 sample rows; co-detection TP halves are excluded from the main-text 7.",
            "All 23 Table 6+7 sample rows have organism, orthogonal/repeat PCR, chart adjudication, final classification, and concise explanation fields.",
            "Pre-DTCA sensitivity 103/110 implies 7 FN sample-equivalents; post-DTCA OPA is 187/191 (97.9%).",
            "Confident Supp Data 1 joins: P-NP-1, P-NP-6, P-NP-113, P-NP-19, P-NP-33, P-NP-64, P-NP-101, P-NP-109, BAL105.",
            "Probable unordered pair: P-NP-99 and P-NP-100 for the two rhinovirus/enterovirus Table 7 ND rows.",
            "Confirmed residual mNGS false negatives after DTCA: RSV (0.5) and metapneumovirus (1.0).",
        ],
        "unresolved_issues": [
            "10 of 14 Table 6 rows (mostly RVP ND) have no sample ID in Supplementary Data 1.",
            "Table 6 row 3 candidate P-NP-36 conflicts on whether mNGS detected RSV.",
            "Table 7 Influenza A+RSV / mNGS ND row has no Supp Data 1 match.",
            "P-NP-99 vs P-NP-100 cannot be uniquely ordered to Table 7 rows 5 vs 6.",
            "Published Zenodo/article files provide no NP/NUC↔P-NP crosswalk; unresolved lab IDs remain flagged in unresolved_identifier_mappings.csv and were not invented.",
        ],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Project root containing data/ and analysis/",
    )
    args = parser.parse_args()
    summary = build_artifacts(args.project_root.resolve())
    print(json.dumps(summary["verification"], indent=2))
    print("Wrote:")
    for key, path in summary["outputs"].items():
        print(f"  {key}: {path}")
    verification = summary["verification"]
    if not verification["table7_TN_matches_paper_5_5"]:
        raise SystemExit("Table 7 TN weight does not match paper 5.5")
    if not verification["table7_TN_plus_FN_matches_maintext_7"]:
        raise SystemExit("Table 7 TN+FN weight does not match main-text 7")


if __name__ == "__main__":
    main()
