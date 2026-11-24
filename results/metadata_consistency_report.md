# Accuracy cohort metadata consistency report

Generated (UTC): 2026-11-24T04:19:37.556577+00:00

## Introduction

This report audits the Accuracy evaluation cohort from the Nature Communications respiratory metagenomic NGS (mNGS) validation study ([s41467-024-51470-y](https://www.nature.com/articles/s41467-024-51470-y); Zenodo [10553379](https://zenodo.org/records/10553379); BioProject [PRJNA1084017](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1084017)). The goal is to verify that publicly deposited sequencing records, clinical/control labels, and article-reported mNGS vs respiratory viral panel (RVP) discrepancies are internally consistent, and to separate true inconsistencies from unresolved identifier gaps or off-panel detections.

## Challenge

Public Accuracy metadata are split across a Zenodo sequencing master sheet, Supplementary Data 1 clinical tables, Supplementary Tables 6–7 discrepancy narratives, and SRA/STAT sequence products. Identifier systems are not fully crosswalked (NP/NUC lab IDs vs article P-NP IDs), Supplementary Data 1 covers fewer sample IDs than the Accuracy clinical cohort, and article discrepancy counts use fractional DTCA weights. An audit must therefore:

1. Enumerate all Purpose=Accuracy sequencing records reproducibly.
2. Reconstruct every Table 6/7 discrepancy without inventing missing IDs.
3. Independently profile sequence content (here: NCBI STAT taxonomy) and compare it to declared labels where joins exist.
4. Report inconsistencies distinctly from QC/low-confidence issues and from metadata incompleteness.

## Software and data

### Software

- Python 3 with `pandas`, `requests`, `openpyxl`
- `analysis/build_accuracy_manifest.py` — Accuracy manifest assembly
- `analysis/reconstruct_discrepancies.py` — Tables 6–7 reconstruction
- `analysis/audit_accuracy_sequences.py` — STAT sequence audit
- `analysis/build_accuracy_audit_report.py` — this report bundle
- NCBI STAT `run_taxonomy` public API (`cluster_name=public`) as sequence-content evidence source

### Primary data sources

- **zenodo_master**: Sequencing record master sheet; Purpose filter and run identifiers — `respmNGS-sequencing-mastersheet.xlsx` ([link](https://zenodo.org/api/records/10553379/files/respmNGS-sequencing-mastersheet.xlsx/content); sha256 `10bfb73ada4fd9a87f8335d532c8a036674f82300527c0f351ff7acd560d8762`)
- **supplementary_data_1**: Clinical and laboratory metadata (Supplementary Data 1) — `41467_2024_51470_MOESM4_ESM.xlsx` ([link](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM4_ESM.xlsx); sha256 `ed87b95b307c0a0e167d5a4779e1a09a190939bad11eb3abcb4fd3afbbd0297f`)
- **supplementary_information**: Supplementary Tables 1-8, including discrepancy Tables 6-7 — `41467_2024_51470_MOESM1_ESM.pdf` ([link](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM1_ESM.pdf); sha256 `3ddd17f680954a6514a7e652f2dec2b06817042f775158812f5af8238a2d0e78`)
- **figure_source_data**: Figure source data; contextual, not complete clinical metadata — `41467_2024_51470_MOESM6_ESM.xlsx` ([link](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM6_ESM.xlsx); sha256 `8ceeaa5aa47cf90df0ec0ac0c6764cc6bdf1cd1d34a586b99ea38c50839bbd3a`)
- **sra_runinfo**: SRA/BioSample accessions and read metadata for sequence retrieval — `PRJNA1084017_runinfo.csv` ([link](https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1084017); sha256 `9ca53d667e8b5dfe5e0492073cc720295d533367e9b9d58f09621e1c24487d2f`)

### Evidence rule

The article calls organisms with ≥3 non-overlapping reads/contigs. This audit uses NCBI STAT `self_count`/`total_count` ≥3 as a **public proxy** for that threshold. It is not a byte-identical reimplementation of the paper pipeline.

## Detailed results

### 1. Cohort inventory

| Metric | Count |
| --- | ---: |
| Accuracy sequencing records | 281 |
| Clinical sequencing records | 233 |
| Positive controls | 24 |
| Negative controls | 24 |
| Rerun/replicate records | 27 |
| Accuracy specimen keys table rows | 206 |
| Expected article clinical specimens | 191 |
| Supplementary Data 1 sample IDs | 109 |
| Unresolved NP/NUC mappings | 79 |

Verified inventory matches the assemble-manifest targets: **281** Accuracy files, with clinical cohort structure reconciling to the article’s **191** specimens (186 unambiguous non-repeat keys + 5 unresolved members of repeat run 039).

### 2. Article discrepancy reconstruction

- Supplementary Table 6 (mNGS+/RVP−): **14** sample rows.
- Supplementary Table 7 sample rows: **9**; TN+FN fractional weight **7.0** matches main-text mNGS-/RVP+ count of **7**.
- Pre-DTCA accuracy reported **179 of 191**; post-DTCA OPA **187 of 191 = 97.906%**.
- ID mapping status among discrepancy rows: unresolved=11, confident=9, probable=2, conflict=1

- Supplementary Table 6 contains exactly 14 mNGS-positive/RVP-negative samples.
- Supplementary Table 7 lists 9 sample rows; fractional TN+FN weights are 5.5+1.5=7, matching the main-text mNGS-/RVP+ count and DTCA summary.
- Pre-DTCA sensitivity 103/110 implies 7 FN sample-equivalents; post-DTCA OPA is 187/191 (97.9%).
- Confident Supp Data 1 joins: P-NP-1, P-NP-6, P-NP-113, P-NP-19, P-NP-33, P-NP-64, P-NP-101, P-NP-109, BAL105.
- Probable unordered pair: P-NP-99 and P-NP-100 for the two rhinovirus/enterovirus Table 7 ND rows.
- Confirmed residual mNGS false negatives after DTCA: RSV (0.5) and metapneumovirus (1.0).

Unresolved discrepancy-ID issues retained without invention:

- 10 of 14 Table 6 rows (mostly RVP ND) have no sample ID in Supplementary Data 1.
- Table 6 row 3 candidate P-NP-36 conflicts on whether mNGS detected RSV.
- Table 7 Influenza A+RSV / mNGS ND row has no Supp Data 1 match.
- P-NP-99 vs P-NP-100 cannot be uniquely ordered to Table 7 rows 5 vs 6.

### 3. Independent STAT sequence audit

- STAT fetch success: **279 / 281**
- Candidate sequence inconsistencies: **4**
- Records overlapping reconstructed article discrepancies: **4**

Consistency class counts:

- `sequence_profiled_no_supp_data1_join`: 105
- `sequence_profiled_metadata_unresolved`: 77
- `consistent_sequence_supports_declared_organisms`: 47
- `consistent_positive_control_virus_detected`: 24
- `consistent_negative_control`: 22
- `audit_incomplete`: 2
- `inconsistent_sequence_vs_control_label`: 2
- `possible_inconsistent_sequence_vs_metadata`: 2

Final adjudication counts (report bundle):

- `metadata_incomplete_not_adjudicated`: 184
- `consistent`: 89
- `article_discrepancy_documented`: 4
- `candidate_inconsistency`: 4

### 4. What counts as inconsistent vs not

- **Inconsistent / discrepancy**: article Tables 6–7 cases, or STAT evidence conflicting with control labels / joined clinical organism expectations.
- **Not automatically inconsistent**: off-panel virus detections; unresolved NP/NUC aliases; Accuracy records absent from Supplementary Data 1; low-level background taxa (anelloviruses/phages) under the evidence proxy.
- **QC / incomplete**: STAT fetch failures (`audit_incomplete`) if any; these are reported separately from mislabel calls.

### 5. Verification against paper headline numbers

| Target | Result |
| --- | --- |
| 281 Accuracy files | 281 audited |
| 191 clinical specimens | article expected 191; manifest clinical records 233 including reruns/controls excluded from 191 |
| 14 mNGS+/RVP− + 7 mNGS-/RVP+ | Table6=14; Table7 TN+FN weight=7.0 (match=True) |
| ~97.9% post-adjudication agreement | 187 of 191 = 97.906% |

## Deliverables

- `results/accuracy_metadata_audit.csv`
- `results/inconsistent_samples.md`
- `results/metadata_consistency_report.md`
- Intermediate sequence audit: `results/intermediate/accuracy_sequence_audit.csv`, `results/intermediate/accuracy_sequence_audit_summary.json`

## Limitations

- STAT k-mer taxonomy counts are a proxy, not a reimplementation of the paper's read/contig pipeline.
- 79 NP/NUC identifiers remain without published P-NP crosswalks; those runs are profiled but not clinically compared.
- Supplementary Data 1 covers 109 sample IDs; many Accuracy clinical records lack joined organism labels.
- Off-panel detections (CMV, papillomavirus, herpesviruses, bocavirus) are reported but not auto-labeled false positives.
- No FASTQ re-analysis was performed; ENA/SRA FASTQ fallback was not required for STAT coverage.
- Supplementary Data 1 contains 109 explicit sample IDs across 134 result lines, while the paper describes 110 virus-positive samples and the normalized sequencing master contains 191 clinical specimens.
- 79 NP/NUC identifiers cannot be paired to the article's P-NP identifiers from the published files alone. They remain stable specimen keys and are not orphaned FASTQ rows.
- The 191 article-cohort count structurally reconciles as 186 unambiguous non-repeat-run specimen keys plus 5 members of the 20 clinical records in explicit repeat run 039. The published metadata do not identify which 5; all 20 are retained and flagged.

