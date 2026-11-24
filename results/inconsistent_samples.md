# Inconsistent samples

Definition used:
- **(a)** article-reported mNGS vs RVP disagreement (Supplementary Tables 6–7), or
- **(b)** sequence content (NCBI STAT) conflicting with declared clinical/control metadata.

QC failures, STAT-unavailable runs, unresolved NP/NUC aliases, and low-confidence/partial STAT support are listed separately and are **not** counted as true mislabels unless they also meet (a) or (b).
mNGS-only detections of viruses not targeted by RVP are **not** treated as inconsistencies.

## A. Article-reported mNGS+/RVP− (n=14)

- **P-NP-1** — Table 6 row 1; RVP[RSV] vs mNGS[RSV; Rhinovirus B]; DTCA: TP RSV; TP Rhinovirus; mapping=confident.
  - Why: Table 6 row 1: RVP[RSV] vs mNGS[RSV; Rhinovirus B]. Orthogonal/repeat testing: (+) Rhinovirus. Final mNGS by DTCA: TP RSV; TP Rhinovirus.
- **P-NP-6** — Table 6 row 2; RVP[RSV] vs mNGS[RSV; EV-D68]; DTCA: TP RSV; TP Enterovirus; mapping=confident.
  - Why: Table 6 row 2: RVP[RSV] vs mNGS[RSV; EV-D68]. Orthogonal/repeat testing: (+) Enterovirus. Final mNGS by DTCA: TP RSV; TP Enterovirus.
- **P-NP-36** — Table 6 row 3; RVP[RSV] vs mNGS[RSV; Adenovirus C]; DTCA: TP RSV; TP Adenovirus; mapping=conflict.
  - Why: Table 6 row 3: RVP[RSV] vs mNGS[RSV; Adenovirus C]. Orthogonal/repeat testing: (+) Adenovirus. Final mNGS by DTCA: TP RSV; TP Adenovirus. ID mapping conflict: Only Supp Data 1 RSV+/Adenovirus-C candidate with Flu A/B/RSV-only comment, but Supp Data 1 lists mNGS 'RSV not detected' while Table 6 and DTCA require mNGS RSV detection (TP RSV). Mapping retained as conflicted candidate only.
- **P-NP-113** — Table 6 row 4; RVP[Adenovirus; Metapneumovirus] vs mNGS[Adenovirus C; Metapneumovirus; Bocavirus]; DTCA: TP Adenovirus; TP Metapneumovirus; TP Bocavirus; mapping=confident.
  - Why: Table 6 row 4: RVP[Adenovirus; Metapneumovirus] vs mNGS[Adenovirus C; Metapneumovirus; Bocavirus]. Orthogonal/repeat testing: Bocavirus ND. Bocavirus mNGS reads are specific to the viral genus. Patient with clinical syndrome consistent with acute viral infection. Final mNGS by DTCA: TP Adenovirus; TP Metapneumovirus; TP Bocavirus.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 5; RVP[ND] vs mNGS[Coronavirus OC43]; DTCA: TP Coronavirus OC43; mapping=unresolved.
  - Why: Table 6 row 5: RVP[ND] vs mNGS[Coronavirus OC43]. Orthogonal/repeat testing: Coronavirus OC43 ND. Detected by mNGS x2 with high RPM and high coverage. No evidence of cross-contamination. Patient with clinical syndrome consistent with acute respiratory viral infection. Final mNGS by DTCA: TP Coronavirus OC43. ID mapping unresolved: RVP-negative case; Supp Data 1 enumerates virus-positive samples and does not identify this negative-cohort specimen.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 6; RVP[ND] vs mNGS[Coronavirus OC43]; DTCA: TP Coronavirus OC43; mapping=unresolved.
  - Why: Table 6 row 6: RVP[ND] vs mNGS[Coronavirus OC43]. Orthogonal/repeat testing: (+) Coronavirus OC43. Final mNGS by DTCA: TP Coronavirus OC43. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 7; RVP[ND] vs mNGS[Rhinovirus A]; DTCA: TP Rhinovirus; mapping=unresolved.
  - Why: Table 6 row 7: RVP[ND] vs mNGS[Rhinovirus A]. Orthogonal/repeat testing: (+) Rhinovirus. Final mNGS by DTCA: TP Rhinovirus. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 8; RVP[ND] vs mNGS[Rhinovirus A]; DTCA: TP Rhinovirus; mapping=unresolved.
  - Why: Table 6 row 8: RVP[ND] vs mNGS[Rhinovirus A]. Orthogonal/repeat testing: (+) Rhinovirus. Final mNGS by DTCA: TP Rhinovirus. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 9; RVP[ND] vs mNGS[Coronavirus HKU1]; DTCA: FP Coronavirus HKU1; mapping=unresolved.
  - Why: Table 6 row 9: RVP[ND] vs mNGS[Coronavirus HKU1]. Orthogonal/repeat testing: Coronavirus HKU1 ND. Final mNGS by DTCA: FP Coronavirus HKU1. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 10; RVP[ND] vs mNGS[Coronavirus HKU1]; DTCA: TP Coronavirus HKU1; mapping=unresolved.
  - Why: Table 6 row 10: RVP[ND] vs mNGS[Coronavirus HKU1]. Orthogonal/repeat testing: (+) Coronavirus HKU1. Final mNGS by DTCA: TP Coronavirus HKU1. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 11; RVP[ND] vs mNGS[Parainfluenza 1; SARS-CoV-2]; DTCA: TP Parainfluenza 1; FP SARS-CoV-2; mapping=unresolved.
  - Why: Table 6 row 11: RVP[ND] vs mNGS[Parainfluenza 1; SARS-CoV-2]. Orthogonal/repeat testing: (+) Parainfluenza 1; SARS-CoV-2 ND. Final mNGS by DTCA: TP Parainfluenza 1; FP SARS-CoV-2. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 12; RVP[ND] vs mNGS[Coronavirus NL63]; DTCA: TP Coronavirus NL63; mapping=unresolved.
  - Why: Table 6 row 12: RVP[ND] vs mNGS[Coronavirus NL63]. Orthogonal/repeat testing: (+) Coronavirus NL63. Final mNGS by DTCA: TP Coronavirus NL63. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 13; RVP[ND] vs mNGS[Rhinovirus C]; DTCA: TP Rhinovirus; mapping=unresolved.
  - Why: Table 6 row 13: RVP[ND] vs mNGS[Rhinovirus C]. Orthogonal/repeat testing: Rhinovirus ND. Mismatches in the forward primer and probe sequences for rhinovirus/enterovirus PCR. No evidence of contamination. Patient with clinical syndrome consistent with acute respiratory viral infection. Related main-text Fig. 7 (12 non-overlapping reads). Final mNGS by DTCA: TP Rhinovirus. ID mapping unresolved: RVP-negative rhinovirus C primer-mismatch case described in text/Fig. 7 without an explicit sample ID in Supp Data 1.
- **(sample ID unresolved in Supp Data 1)** — Table 6 row 14; RVP[ND] vs mNGS[Coronavirus OC43]; DTCA: NA Coronavirus OC43; mapping=unresolved.
  - Why: Table 6 row 14: RVP[ND] vs mNGS[Coronavirus OC43]. Orthogonal/repeat testing: NT. Final mNGS by DTCA: NA Coronavirus OC43. ID mapping unresolved: RVP-negative case; orthogonal testing not performed; no sample-ID crosswalk in Supp Data 1.

## B. Article-reported mNGS−/RVP+ main-text equivalents (n=7; Table 7 has 9 sample rows)

Main-text 7 = Table 7 TN weight 5.5 + FN weight 1.5.

- **P-NP-19** — Table 7 row 1; RVP[Influenza A; Parainfluenza 3] vs mNGS[Influenza A]; DTCA: TP Influenza A; TN Parainfluenza 3; mapping=confident.
  - Why: Table 7 row 1: RVP[Influenza A; Parainfluenza 3] vs mNGS[Influenza A]. Orthogonal/repeat testing: (+) Influenza A; Parainfluenza 3 ND. Platform: MagPix. Final mNGS by DTCA: TP Influenza A; TN Parainfluenza 3.
- **P-NP-33** — Table 7 row 2; RVP[Influenza A; SARS-CoV-2] vs mNGS[Influenza A]; DTCA: TP Influenza A; TN SARS-CoV-2; mapping=confident.
  - Why: Table 7 row 2: RVP[Influenza A; SARS-CoV-2] vs mNGS[Influenza A]. Orthogonal/repeat testing: Influenza A NT; SARS-CoV-2 ND. Platform: Verigene. Final mNGS by DTCA: TP Influenza A; TN SARS-CoV-2.
- **(sample ID unresolved in Supp Data 1)** — Table 7 row 3; RVP[Influenza A; RSV] vs mNGS[ND]; DTCA: TN Influenza A; FN RSV; mapping=unresolved.
  - Why: Table 7 row 3: RVP[Influenza A; RSV] vs mNGS[ND]. Orthogonal/repeat testing: Influenza A ND; (+) RSV B. Platform: GenMark. Final mNGS by DTCA: TN Influenza A; FN RSV. ID mapping unresolved: No Supp Data 1 sample with both Influenza A and RSV PCR-positive and mNGS ND. Not inventing a crosswalk.
- **P-NP-64** — Table 7 row 4; RVP[Influenza A; SARS-CoV-2] vs mNGS[SARS-CoV-2]; DTCA: TP SARS-CoV-2; TN Influenza A; mapping=confident.
  - Why: Table 7 row 4: RVP[Influenza A; SARS-CoV-2] vs mNGS[SARS-CoV-2]. Orthogonal/repeat testing: Influenza A ND; SARS-CoV-2 NT. Platform: MagPix. Final mNGS by DTCA: TP SARS-CoV-2; TN Influenza A.
- **P-NP-99** — Table 7 row 5; RVP[Rhinovirus/Enterovirus] vs mNGS[ND]; DTCA: TN Rhinovirus/Enterovirus; mapping=probable.
  - Why: Table 7 row 5: RVP[Rhinovirus/Enterovirus] vs mNGS[ND]. Orthogonal/repeat testing: Rhinovirus/Enterovirus ND. Platform: MagPix. Final mNGS by DTCA: TN Rhinovirus/Enterovirus. ID mapping probable: Supp Data 1 Rhinovirus/Enterovirus+ with mNGS ND and comment 'False Negative'; Table 7 has two such rows (MFI 52 vs NA). Assigned to P-NP-99/P-NP-100 as an unordered pair; MFI cannot distinguish them from published files alone.
- **P-NP-100** — Table 7 row 6; RVP[Rhinovirus/Enterovirus] vs mNGS[ND]; DTCA: TN Rhinovirus/Enterovirus; mapping=probable.
  - Why: Table 7 row 6: RVP[Rhinovirus/Enterovirus] vs mNGS[ND]. Orthogonal/repeat testing: Rhinovirus/Enterovirus ND. Platform: MagPix. Final mNGS by DTCA: TN Rhinovirus/Enterovirus. ID mapping probable: Paired with Table 7 row 5 / P-NP-99; same Supp Data 1 pattern. Row-to-ID order is not uniquely determined.
- **P-NP-101** — Table 7 row 7; RVP[Rhinovirus/Enterovirus; Parainfluenza 2] vs mNGS[Rhinovirus A]; DTCA: TP Rhinovirus/Enterovirus; TN Parainfluenza 2; mapping=confident.
  - Why: Table 7 row 7: RVP[Rhinovirus/Enterovirus; Parainfluenza 2] vs mNGS[Rhinovirus A]. Orthogonal/repeat testing: (+) Rhinovirus/Enterovirus; Parainfluenza 2 ND. Platform: MagPix. Final mNGS by DTCA: TP Rhinovirus/Enterovirus; TN Parainfluenza 2.
- **P-NP-109** — Table 7 row 8; RVP[Metapneumovirus] vs mNGS[ND]; DTCA: FN Metapneumovirus; mapping=confident.
  - Why: Table 7 row 8: RVP[Metapneumovirus] vs mNGS[ND]. Orthogonal/repeat testing: (+) Metapneumovirus. Platform: GenMark. Final mNGS by DTCA: FN Metapneumovirus.
- **BAL105** — Table 7 row 9; RVP[RSV A] vs mNGS[ND]; DTCA: TN RSV A; mapping=confident.
  - Why: Table 7 row 9: RVP[RSV A] vs mNGS[ND]. Orthogonal/repeat testing: RSV ND. Platform: MagPix. Final mNGS by DTCA: TN RSV A.

### Main-text 7 TN/FN components

- P-NP-19: Parainfluenza 3 → TN (weight 0.5; reclassified_true_negative)
- P-NP-33: SARS-CoV-2 → TN (weight 0.5; reclassified_true_negative)
- (unresolved): Influenza A → TN (weight 0.5; reclassified_true_negative)
- (unresolved): RSV → FN (weight 0.5; residual_false_negative)
- P-NP-64: Influenza A → TN (weight 0.5; reclassified_true_negative)
- P-NP-99: Rhinovirus/Enterovirus → TN (weight 1.0; reclassified_true_negative)
- P-NP-100: Rhinovirus/Enterovirus → TN (weight 1.0; reclassified_true_negative)
- P-NP-101: Parainfluenza 2 → TN (weight 0.5; reclassified_true_negative)
- P-NP-109: Metapneumovirus → FN (weight 1.0; residual_false_negative)
- BAL105: RSV A → TN (weight 1.0; reclassified_true_negative)

## C. Sequence-vs-metadata conflicts (STAT audit)

- **NC-2** / `UCSF_VRP_mNGS_039_0891_nohuman.fastq.gz` / SRR28234621 (negative_control)
  - Why: Negative control has threshold STAT hits for: SARS_CoV_2. Interpreted as possible contamination or mislabel; not treated as QC failure alone.
  - STAT passing targets: SARS_CoV_2; flags: nc_pathogen_signal
- **PC-1** / `UCSF_VRP_mNGS_048_1085_nohuman.fastq.gz` / SRR28242212 (positive_control)
  - Why: Positive control missing threshold STAT for: Influenza_A, Influenza_B, RSV, SARS_CoV_2.
  - STAT passing targets: nan; flags: pc_missing_target
- **PC-2** / `UCSF_VRP_mNGS_048_1086_nohuman.fastq.gz` / SRR28242201 (positive_control)
  - Why: Positive control missing threshold STAT for: Influenza_A, Influenza_B, RSV.
  - STAT passing targets: SARS_CoV_2; flags: pc_missing_target
- **P-NP-14** / `UCSF_VRP_mNGS_050_1109_nohuman.fastq.gz` / SRR28241899 (clinical_sample)
  - Why: Declared targets ['RSV'] not supported by STAT; STAT passing=[].
  - STAT passing targets: nan; flags: declared_vs_stat_mismatch
- **BAL139** / `UCSF_VRP_mNGS_067_1335_nohuman.fastq.gz` / SRR28242341 (clinical_sample)
  - Why: Declared targets ['Rhinovirus_Enterovirus'] not supported by STAT; STAT passing=[].
  - STAT passing targets: nan; flags: declared_vs_stat_mismatch
- **NC-2** / `UCSF_VRP_mNGS_074_1500_nohuman.fastq.gz` / SRR28242480 (negative_control)
  - Why: Negative control has threshold STAT hits for: Adenovirus. Interpreted as possible contamination or mislabel; not treated as QC failure alone.
  - STAT passing targets: Adenovirus; flags: nc_pathogen_signal

## D. Separated: low-confidence / partial / unavailable (not automatic mislabels)

- **NP363** / SRR28242333 — status=stat_unavailable; No NCBI STAT taxonomy payload for this SRR.
- **NP270** / SRR28241792 — status=stat_unavailable; No NCBI STAT taxonomy payload for this SRR.

## E. Unresolved identifier limitations

- Unresolved NP/NUC↔article crosswalks in Accuracy cohort: see `data/normalized/unresolved_identifier_mappings.csv`.
- Discrepancy rows lacking publishable sample IDs: 14 (see `article_discrepancies_unresolved_ids.csv`).
- No NP/NUC↔P-NP crosswalks were invented.
