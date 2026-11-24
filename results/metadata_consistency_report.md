# Metadata consistency audit: UCSF viral respiratory mNGS Accuracy cohort

**Paper:** Tan et al., *Nature Communications* (2024) — [Laboratory validation of a clinical metagenomic next-generation sequencing assay for respiratory virus detection and discovery](https://www.nature.com/articles/s41467-024-51470-y) ([doi:10.1038/s41467-024-51470-y](https://doi.org/10.1038/s41467-024-51470-y))  
**BioProject:** [PRJNA1084017](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1084017)  
**Zenodo sequencing master sheet:** [10.5281/zenodo.10553379](https://zenodo.org/records/10553379)  
**Generated (UTC):** 2026-11-24T04:16:17.988612+00:00

## 1. Introduction

Respiratory infections remain a leading cause of morbidity, yet conventional multiplex RT-PCR panels (RVP) only interrogate a fixed pathogen list and can fail when primers drift or when uncommon/novel viruses are present. Metagenomic next-generation sequencing (mNGS) offers agnostic detection of known and divergent viruses from the same nucleic-acid extract. Tan et al. validated a largely automated clinical respiratory viral mNGS assay (SURPI+) on residual clinical specimens and reported strong agreement with RVP that further improved after discrepancy testing and clinical adjudication (DTCA).

This audit revisits the Accuracy evaluation cohort as deposited in public repositories. The goal is **metadata consistency**: whether each deposited sequencing record’s declared clinical/control identity and article-reported mNGS/RVP calls are coherent with (i) published supplementary discrepancy tables and (ii) independent sequence-derived taxonomy evidence.

## 2. The challenge

The Accuracy set mixes clinical specimens, positive/negative controls, and reruns under lab identifiers (NP/NUC) that only partially crosswalk to article sample IDs (P-NP/BAL). Supplementary Data 1 lists virus-positive clinical metadata for 109 explicit IDs, while the article Accuracy denominator is **191 clinical specimens**. Supplementary Tables 6–7 enumerate **14 mNGS+/RVP−** and (main-text) **7 mNGS−/RVP+** discrepancy equivalents, but many Table 6 RVP-negative rows lack publishable sample IDs. Independently confirming biological content for all **281** Purpose=Accuracy FASTQs (~11 GB) is constrained without local Kraken/SURPI databases; NCBI’s public STAT taxonomy API provides sequence-derived k-mer taxonomy for nearly all SRRs and is used here as an independent check.

**Inconsistency definition used here**
1. **(a)** article-reported mNGS vs RVP disagreement (Tables 6–7 / DTCA), or  
2. **(b)** sequence content conflicting with declared clinical/control metadata.  

QC failures, STAT-unavailable runs, unresolved aliases, and low-confidence/partial hits are separated from true mislabels. mNGS-only detections of viruses **not targeted by RVP** are not counted as inconsistencies.

## 3. Software and data

### Data sources
| Source | Role |
| --- | --- |
| Zenodo `respmNGS-sequencing-mastersheet.xlsx` | Purpose=Accuracy filter; 281 sequencing records |
| Supplementary Data 1 (`41467_2024_51470_MOESM4_ESM.xlsx`) | Clinical PCR/mNGS/viral-load metadata |
| Supplementary Information PDF (`...MOESM1_ESM.pdf`) | Tables 6–7 discrepancy adjudication |
| SRA RunInfo for PRJNA1084017 | SRR/BioSample linkage and read metrics |
| NCBI STAT `run_taxonomy` API | Independent taxonomy for each SRR |

### Software / scripts in this repository
- `analysis/build_accuracy_manifest.py` — Accuracy manifest assembly  
- `analysis/reconstruct_discrepancies.py` — Tables 6–7 reconstruction + DTCA case cards  
- `analysis/audit_sequences_stat.py` — STAT harvest/classify for all Accuracy SRRs  
- `analysis/audit_accuracy_fastq.py` — fuller FASTQ/STAT audit scaffold (cloud handoff)  
- Outputs under `data/normalized/` and `results/`

### Sequence-evidence rule
Paper positive criterion: **≥3 non-overlapping viral reads/contigs**. This audit uses NCBI STAT **`total_count ≥ 3`** as a **proxy** (k-mer/spot taxonomy counts, not SURPI+ non-overlapping coverage). Evidence type is labeled `sequence_STAT_kmer_counts`. Deposited libraries are human-depleted (`nohuman`).

## 4. Results

### 4.1 Cohort inventory (verified)
| Quantity | Count |
| --- | ---: |
| Accuracy sequencing files/records | **281** |
| Clinical sequencing records | 233 |
| Positive controls | 24 |
| Negative controls | 24 |
| Rerun/replicate records | 27 |
| Article clinical specimens (Accuracy denominator) | **191** (186 unambiguous keys + 5 unresolved members of repeat run 039) |
| STAT taxonomy available | 279 / 281 |
| STAT fetch failures | 2 |

Manifest provenance counts match the assembled Accuracy filter (`Purpose == "Accuracy"`).

### 4.2 Article discrepancy reconstruction (definition a)
Exact reconstructed counts:
- **14** mNGS+/RVP− samples (Supplementary Table 6)
- **7** mNGS−/RVP+ **main-text equivalents** (Table 7 publishes 9 sample rows; TN 5.5 + FN 1.5 = 7.0)
- Pre-DTCA accuracy: **179 / 191 = 93.7%** (sensitivity 103/110; specificity 76/81)
- Post-DTCA overall predictive agreement: **187 / 191 = 97.9%** (PPA 110.5/113; NPA 76.5/78)

DTCA residual mNGS false negatives after adjudication: RSV (0.5) + metapneumovirus (1.0).  
ID mapping from Tables 6–7 into Supp Data 1: confident 9, probable 2, conflict 1 (P-NP-36), unresolved 11 — **no NP/NUC↔P-NP crosswalks invented**.

Structured artifacts: `supplementary_table6_discrepancies.csv`, `supplementary_table7_discrepancies.csv`, `article_discrepancies_*.csv`, `results/intermediate/discrepancy_case_cards.json`.

### 4.3 Independent STAT sequence audit (definition b)
All **281** Accuracy SRRs were queried via NCBI STAT. **279** returned taxonomy; **2** returned empty payloads (`SRR28242333`, `SRR28241792`).

Final consistency status among sequencing records:
- `unresolved_identifier_or_incomplete_metadata`: **182**
- `consistent`: **87**
- `inconsistent_sequence_vs_metadata`: **6**
- `article_discrepancy_adjudicated`: **4**
- `stat_unavailable`: **2**

Controls:
- Positive controls with all four Accuplex targets at threshold: **22 / 24**
- Negative controls without threshold respiratory-pathogen hits: **22 / 24**

Sequence-vs-metadata conflicts (status `inconsistent_sequence_vs_metadata`): **6** records (see `results/inconsistent_samples.md` section C). These are candidate contamination/mislabels or missing expected pathogen signal under the STAT proxy and should be reviewed separately from article DTCA cases.

Many clinical records remain `unresolved_identifier_or_incomplete_metadata` because 79 NP/NUC lab IDs lack a published article crosswalk; presence/absence of STAT pathogen signal in those rows is reported but **not** forced into a mislabel call.

### 4.4 Reconciliation to ~97.9% post-adjudication agreement
The article’s Accuracy denominator is specimen-level (**191**), not FASTQ-level (**281**). Using DTCA as the composite reference:
- Concordant after adjudication: **187**
- Discordant after adjudication: **4** specimen-equivalents (191−187)
- OPA = 187/191 ≈ **97.906% ≈ 97.9%**

This audit reproduces those DTCA counts from Supplementary Tables 6–7 and does not invent additional specimen-level concordance beyond published adjudication. Sequence-level STAT conflicts are additional **repository metadata/sequence** findings and are not automatically subtracted from the article’s 97.9% clinical OPA.

### 4.5 Key inconsistent samples (summary)
See **`results/inconsistent_samples.md`** for the full list. Highlights:
- Mapped Table 6/7 examples: P-NP-1, P-NP-6, P-NP-113 (mNGS+/RVP−); P-NP-19, P-NP-33, P-NP-64, P-NP-99/100, P-NP-101, P-NP-109, BAL105 (mNGS−/RVP+ pathway).
- Residual confirmed FN after DTCA: **P-NP-109** (metapneumovirus) and the unresolved Influenza A+RSV dual-adjudication sample (0.5 FN RSV).
- STAT control anomalies: two PC records missing expected Accuplex targets; two NC records with threshold SARS-CoV-2 or Adenovirus signal.
- Clinical STAT mismatches without article discrepancy join: e.g. P-NP-14 (RSV declared, no STAT RSV), BAL139 (rhinovirus declared, no STAT rhinovirus/enterovirus).

## 5. Limitations
1. STAT k-mer `total_count` ≠ SURPI+ non-overlapping genomic coverage; threshold ≥3 is a proxy.  
2. Two SRRs lack STAT payloads.  
3. 79 NP/NUC identifiers remain without published P-NP crosswalks.  
4. Most Table 6 RVP-negative discrepancy rows have no Supp Data 1 sample ID.  
5. Absorbed base-branch partial STAT intermediates (`results/intermediate/sequence_audit/`, commit bf5fb30) and completed a full 281-run STAT API audit on this branch. Sibling cloud sequence-audit agent `bc-c112c72e-7a68-4840-978f-0fdb051d5d1c` had not yet published a completion branch at PR time.  
6. No fabricated identifier mappings.

## 6. Deliverables
- `results/accuracy_metadata_audit.csv` — per-record audit table (281 rows)
- `results/inconsistent_samples.md` — sample IDs and why inconsistent
- `results/metadata_consistency_report.md` — this document
- Supporting discrepancy + STAT intermediates under `data/normalized/` and `results/intermediate/`
