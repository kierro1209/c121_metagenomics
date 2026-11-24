# Inconsistent / discrepancy samples

Generated (UTC): 2026-11-24T04:19:37.533401+00:00

## Definition

Inconsistency here means either (a) an article-reported mNGS vs RVP disagreement documented in Supplementary Tables 6–7, or (b) STAT sequence evidence conflicting with declared clinical/control metadata under a >=3-count proxy for non-overlapping viral evidence. Off-panel mNGS detections (for example CMV, papillomavirus, herpesviruses, bocavirus when not an RVP target) are **not** automatically labeled inconsistent. Unresolved NP/NUC↔P-NP identifiers are profiled but not forced into a fabricated clinical join.

## Summary counts

- Accuracy records audited: **281**
- Article discrepancy rows reconstructed: Table 6 = **14**, Table 7 = **9** (main-text mNGS-/RVP+ = 7 via TN+FN weights 7.0)
- Sequence candidate conflicts (STAT vs metadata/controls): **4**
- Rows flagged article discrepancy and/or sequence candidate: **8**

## A. Article-reported discrepancies (Supplementary Tables 6–7)

### Table 6 row 1: P-NP-1

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `confident`
- RVP organisms: RSV
- mNGS organisms: RSV; Rhinovirus B
- Orthogonal / repeat testing: (+) Rhinovirus
- Final mNGS by DTCA: **TP RSV; TP Rhinovirus**
- Explanation: Table 6 row 1: RVP[RSV] vs mNGS[RSV; Rhinovirus B]. Orthogonal/repeat testing: (+) Rhinovirus. Final mNGS by DTCA: TP RSV; TP Rhinovirus.

### Table 6 row 2: P-NP-6

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `confident`
- RVP organisms: RSV
- mNGS organisms: RSV; EV-D68
- Orthogonal / repeat testing: (+) Enterovirus
- Final mNGS by DTCA: **TP RSV; TP Enterovirus**
- Explanation: Table 6 row 2: RVP[RSV] vs mNGS[RSV; EV-D68]. Orthogonal/repeat testing: (+) Enterovirus. Final mNGS by DTCA: TP RSV; TP Enterovirus.

### Table 6 row 3: P-NP-36

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `conflict`
- RVP organisms: RSV
- mNGS organisms: RSV; Adenovirus C
- Orthogonal / repeat testing: (+) Adenovirus
- Final mNGS by DTCA: **TP RSV; TP Adenovirus**
- Explanation: Table 6 row 3: RVP[RSV] vs mNGS[RSV; Adenovirus C]. Orthogonal/repeat testing: (+) Adenovirus. Final mNGS by DTCA: TP RSV; TP Adenovirus. ID mapping conflict: Only Supp Data 1 RSV+/Adenovirus-C candidate with Flu A/B/RSV-only comment, but Supp Data 1 lists mNGS 'RSV not detected' while Table 6 and DTCA require mNGS RSV detection (TP RSV). Mapping retained as conflicted candidate only.
- STAT respiratory keys (SRR28241870): `adenovirus`
- Sequence consistency class: `consistent_sequence_supports_declared_organisms`

### Table 6 row 4: P-NP-113

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `confident`
- RVP organisms: Adenovirus; Metapneumovirus
- mNGS organisms: Adenovirus C; Metapneumovirus; Bocavirus
- Orthogonal / repeat testing: Bocavirus ND
- Final mNGS by DTCA: **TP Adenovirus; TP Metapneumovirus; TP Bocavirus**
- Explanation: Table 6 row 4: RVP[Adenovirus; Metapneumovirus] vs mNGS[Adenovirus C; Metapneumovirus; Bocavirus]. Orthogonal/repeat testing: Bocavirus ND. Bocavirus mNGS reads are specific to the viral genus. Patient with clinical syndrome consistent with acute viral infection. Final mNGS by DTCA: TP Adenovirus; TP Metapneumovirus; TP Bocavirus.

### Table 6 row 5: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Coronavirus OC43
- Orthogonal / repeat testing: Coronavirus OC43 ND
- Final mNGS by DTCA: **TP Coronavirus OC43**
- Explanation: Table 6 row 5: RVP[ND] vs mNGS[Coronavirus OC43]. Orthogonal/repeat testing: Coronavirus OC43 ND. Detected by mNGS x2 with high RPM and high coverage. No evidence of cross-contamination. Patient with clinical syndrome consistent with acute respiratory viral infection. Final mNGS by DTCA: TP Coronavirus OC43. ID mapping unresolved: RVP-negative case; Supp Data 1 enumerates virus-positive samples and does not identify this negative-cohort specimen.

### Table 6 row 6: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Coronavirus OC43
- Orthogonal / repeat testing: (+) Coronavirus OC43
- Final mNGS by DTCA: **TP Coronavirus OC43**
- Explanation: Table 6 row 6: RVP[ND] vs mNGS[Coronavirus OC43]. Orthogonal/repeat testing: (+) Coronavirus OC43. Final mNGS by DTCA: TP Coronavirus OC43. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.

### Table 6 row 7: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Rhinovirus A
- Orthogonal / repeat testing: (+) Rhinovirus
- Final mNGS by DTCA: **TP Rhinovirus**
- Explanation: Table 6 row 7: RVP[ND] vs mNGS[Rhinovirus A]. Orthogonal/repeat testing: (+) Rhinovirus. Final mNGS by DTCA: TP Rhinovirus. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.

### Table 6 row 8: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Rhinovirus A
- Orthogonal / repeat testing: (+) Rhinovirus
- Final mNGS by DTCA: **TP Rhinovirus**
- Explanation: Table 6 row 8: RVP[ND] vs mNGS[Rhinovirus A]. Orthogonal/repeat testing: (+) Rhinovirus. Final mNGS by DTCA: TP Rhinovirus. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.

### Table 6 row 9: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Coronavirus HKU1
- Orthogonal / repeat testing: Coronavirus HKU1 ND
- Final mNGS by DTCA: **FP Coronavirus HKU1**
- Explanation: Table 6 row 9: RVP[ND] vs mNGS[Coronavirus HKU1]. Orthogonal/repeat testing: Coronavirus HKU1 ND. Final mNGS by DTCA: FP Coronavirus HKU1. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.

### Table 6 row 10: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Coronavirus HKU1
- Orthogonal / repeat testing: (+) Coronavirus HKU1
- Final mNGS by DTCA: **TP Coronavirus HKU1**
- Explanation: Table 6 row 10: RVP[ND] vs mNGS[Coronavirus HKU1]. Orthogonal/repeat testing: (+) Coronavirus HKU1. Final mNGS by DTCA: TP Coronavirus HKU1. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.

### Table 6 row 11: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Parainfluenza 1; SARS-CoV-2
- Orthogonal / repeat testing: (+) Parainfluenza 1; SARS-CoV-2 ND
- Final mNGS by DTCA: **TP Parainfluenza 1; FP SARS-CoV-2**
- Explanation: Table 6 row 11: RVP[ND] vs mNGS[Parainfluenza 1; SARS-CoV-2]. Orthogonal/repeat testing: (+) Parainfluenza 1; SARS-CoV-2 ND. Final mNGS by DTCA: TP Parainfluenza 1; FP SARS-CoV-2. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.

### Table 6 row 12: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Coronavirus NL63
- Orthogonal / repeat testing: (+) Coronavirus NL63
- Final mNGS by DTCA: **TP Coronavirus NL63**
- Explanation: Table 6 row 12: RVP[ND] vs mNGS[Coronavirus NL63]. Orthogonal/repeat testing: (+) Coronavirus NL63. Final mNGS by DTCA: TP Coronavirus NL63. ID mapping unresolved: RVP-negative case; no published sample-ID crosswalk in Supp Data 1.

### Table 6 row 13: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Rhinovirus C
- Orthogonal / repeat testing: Rhinovirus ND
- Final mNGS by DTCA: **TP Rhinovirus**
- Explanation: Table 6 row 13: RVP[ND] vs mNGS[Rhinovirus C]. Orthogonal/repeat testing: Rhinovirus ND. Mismatches in the forward primer and probe sequences for rhinovirus/enterovirus PCR. No evidence of contamination. Patient with clinical syndrome consistent with acute respiratory viral infection. Related main-text Fig. 7 (12 non-overlapping reads). Final mNGS by DTCA: TP Rhinovirus. ID mapping unresolved: RVP-negative rhinovirus C primer-mismatch case described in text/Fig. 7 without an explicit sample ID in Supp Data 1.

### Table 6 row 14: (unresolved sample id)

- Class: `mNGS_positive_RVP_negative`
- ID mapping: `unresolved`
- RVP organisms: ND
- mNGS organisms: Coronavirus OC43
- Orthogonal / repeat testing: NT
- Final mNGS by DTCA: **NA Coronavirus OC43**
- Explanation: Table 6 row 14: RVP[ND] vs mNGS[Coronavirus OC43]. Orthogonal/repeat testing: NT. Final mNGS by DTCA: NA Coronavirus OC43. ID mapping unresolved: RVP-negative case; orthogonal testing not performed; no sample-ID crosswalk in Supp Data 1.

### Table 7 row 1: P-NP-19

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `confident`
- RVP organisms: Influenza A; Parainfluenza 3
- mNGS organisms: Influenza A
- Orthogonal / repeat testing: (+) Influenza A; Parainfluenza 3 ND
- Final mNGS by DTCA: **TP Influenza A; TN Parainfluenza 3**
- Explanation: Table 7 row 1: RVP[Influenza A; Parainfluenza 3] vs mNGS[Influenza A]. Orthogonal/repeat testing: (+) Influenza A; Parainfluenza 3 ND. Platform: MagPix. Final mNGS by DTCA: TP Influenza A; TN Parainfluenza 3.
- STAT respiratory keys (SRR28241893): `influenza_a`
- Sequence consistency class: `consistent_sequence_supports_declared_organisms`

### Table 7 row 2: P-NP-33

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `confident`
- RVP organisms: Influenza A; SARS-CoV-2
- mNGS organisms: Influenza A
- Orthogonal / repeat testing: Influenza A NT; SARS-CoV-2 ND
- Final mNGS by DTCA: **TP Influenza A; TN SARS-CoV-2**
- Explanation: Table 7 row 2: RVP[Influenza A; SARS-CoV-2] vs mNGS[Influenza A]. Orthogonal/repeat testing: Influenza A NT; SARS-CoV-2 ND. Platform: Verigene. Final mNGS by DTCA: TP Influenza A; TN SARS-CoV-2.
- STAT respiratory keys (SRR28241873): `influenza_a`
- Sequence consistency class: `consistent_sequence_supports_declared_organisms`

### Table 7 row 3: (unresolved sample id)

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `unresolved`
- RVP organisms: Influenza A; RSV
- mNGS organisms: ND
- Orthogonal / repeat testing: Influenza A ND; (+) RSV B
- Final mNGS by DTCA: **TN Influenza A; FN RSV**
- Explanation: Table 7 row 3: RVP[Influenza A; RSV] vs mNGS[ND]. Orthogonal/repeat testing: Influenza A ND; (+) RSV B. Platform: GenMark. Final mNGS by DTCA: TN Influenza A; FN RSV. ID mapping unresolved: No Supp Data 1 sample with both Influenza A and RSV PCR-positive and mNGS ND. Not inventing a crosswalk.

### Table 7 row 4: P-NP-64

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `confident`
- RVP organisms: Influenza A; SARS-CoV-2
- mNGS organisms: SARS-CoV-2
- Orthogonal / repeat testing: Influenza A ND; SARS-CoV-2 NT
- Final mNGS by DTCA: **TP SARS-CoV-2; TN Influenza A**
- Explanation: Table 7 row 4: RVP[Influenza A; SARS-CoV-2] vs mNGS[SARS-CoV-2]. Orthogonal/repeat testing: Influenza A ND; SARS-CoV-2 NT. Platform: MagPix. Final mNGS by DTCA: TP SARS-CoV-2; TN Influenza A.

### Table 7 row 5: P-NP-99

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `probable`
- RVP organisms: Rhinovirus/Enterovirus
- mNGS organisms: ND
- Orthogonal / repeat testing: Rhinovirus/Enterovirus ND
- Final mNGS by DTCA: **TN Rhinovirus/Enterovirus**
- Explanation: Table 7 row 5: RVP[Rhinovirus/Enterovirus] vs mNGS[ND]. Orthogonal/repeat testing: Rhinovirus/Enterovirus ND. Platform: MagPix. Final mNGS by DTCA: TN Rhinovirus/Enterovirus. ID mapping probable: Supp Data 1 Rhinovirus/Enterovirus+ with mNGS ND and comment 'False Negative'; Table 7 has two such rows (MFI 52 vs NA). Assigned to P-NP-99/P-NP-100 as an unordered pair; MFI cannot distinguish them from published files alone.

### Table 7 row 6: P-NP-100

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `probable`
- RVP organisms: Rhinovirus/Enterovirus
- mNGS organisms: ND
- Orthogonal / repeat testing: Rhinovirus/Enterovirus ND
- Final mNGS by DTCA: **TN Rhinovirus/Enterovirus**
- Explanation: Table 7 row 6: RVP[Rhinovirus/Enterovirus] vs mNGS[ND]. Orthogonal/repeat testing: Rhinovirus/Enterovirus ND. Platform: MagPix. Final mNGS by DTCA: TN Rhinovirus/Enterovirus. ID mapping probable: Paired with Table 7 row 5 / P-NP-99; same Supp Data 1 pattern. Row-to-ID order is not uniquely determined.

### Table 7 row 7: P-NP-101

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `confident`
- RVP organisms: Rhinovirus/Enterovirus; Parainfluenza 2
- mNGS organisms: Rhinovirus A
- Orthogonal / repeat testing: (+) Rhinovirus/Enterovirus; Parainfluenza 2 ND
- Final mNGS by DTCA: **TP Rhinovirus/Enterovirus; TN Parainfluenza 2**
- Explanation: Table 7 row 7: RVP[Rhinovirus/Enterovirus; Parainfluenza 2] vs mNGS[Rhinovirus A]. Orthogonal/repeat testing: (+) Rhinovirus/Enterovirus; Parainfluenza 2 ND. Platform: MagPix. Final mNGS by DTCA: TP Rhinovirus/Enterovirus; TN Parainfluenza 2.

### Table 7 row 8: P-NP-109

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `confident`
- RVP organisms: Metapneumovirus
- mNGS organisms: ND
- Orthogonal / repeat testing: (+) Metapneumovirus
- Final mNGS by DTCA: **FN Metapneumovirus**
- Explanation: Table 7 row 8: RVP[Metapneumovirus] vs mNGS[ND]. Orthogonal/repeat testing: (+) Metapneumovirus. Platform: GenMark. Final mNGS by DTCA: FN Metapneumovirus.

### Table 7 row 9: BAL105

- Class: `mNGS_negative_RVP_positive`
- ID mapping: `confident`
- RVP organisms: RSV A
- mNGS organisms: ND
- Orthogonal / repeat testing: RSV ND
- Final mNGS by DTCA: **TN RSV A**
- Explanation: Table 7 row 9: RVP[RSV A] vs mNGS[ND]. Orthogonal/repeat testing: RSV ND. Platform: MagPix. Final mNGS by DTCA: TN RSV A.
- STAT respiratory keys (SRR28242349): `nan`
- Sequence consistency class: `consistent_sequence_supports_declared_organisms`

## B. Sequence vs metadata / control candidate conflicts

### NC-2 / SRR28234621 (negative_control)

- File: `UCSF_VRP_mNGS_039_0891_nohuman.fastq.gz`
- Article sample id: ``
- Consistency class: `inconsistent_sequence_vs_control_label`
- Expected keys: ``
- STAT respiratory keys: `sars_cov_2`
- Notes: Negative control has respiratory pathogen STAT evidence >=3: sars_cov_2(Severe acute respiratory syndrome coronavirus 2:20)

### P-NP-14 / SRR28241899 (clinical_sample)

- File: `UCSF_VRP_mNGS_050_1109_nohuman.fastq.gz`
- Article sample id: `P-NP-14`
- Consistency class: `possible_inconsistent_sequence_vs_metadata`
- Expected keys: `rsv`
- STAT respiratory keys: ``
- Notes: Declared organisms lacking STAT evidence >=3: rsv; found=none Below-threshold near-miss taxa: Orthopneumovirus:1/1

### BAL139 / SRR28242341 (clinical_sample)

- File: `UCSF_VRP_mNGS_067_1335_nohuman.fastq.gz`
- Article sample id: `BAL139`
- Consistency class: `possible_inconsistent_sequence_vs_metadata`
- Expected keys: `rhinovirus`
- STAT respiratory keys: ``
- Notes: Declared organisms lacking STAT evidence >=3: rhinovirus; found=none Below-threshold near-miss taxa: Rhinovirus A:2/2

### NC-2 / SRR28242480 (negative_control)

- File: `UCSF_VRP_mNGS_074_1500_nohuman.fastq.gz`
- Article sample id: ``
- Consistency class: `inconsistent_sequence_vs_control_label`
- Expected keys: ``
- STAT respiratory keys: `adenovirus`
- Notes: Negative control has respiratory pathogen STAT evidence >=3: adenovirus(Mastadenovirus:3)

## C. Explicitly unresolved / not forced

- Unresolved NP/NUC↔article identifier mappings retained: **79**
- Reconstructed discrepancy sample IDs with no Accuracy-manifest `article_sample_id` join (cannot attach an SRR without inventing a crosswalk): **P-NP-1, P-NP-100, P-NP-101, P-NP-109, P-NP-113, P-NP-6, P-NP-64, P-NP-99**
- 10 of 14 Table 6 rows (mostly RVP ND) have no sample ID in Supplementary Data 1.
- Table 6 row 3 candidate P-NP-36 conflicts on whether mNGS detected RSV.
- Table 7 Influenza A+RSV / mNGS ND row has no Supp Data 1 match.
- P-NP-99 vs P-NP-100 cannot be uniquely ordered to Table 7 rows 5 vs 6.

## D. STAT coverage gaps

- `SRR28242333` / `NP363` / `UCSF_VRP_mNGS_068_1364_nohuman.fastq.gz`: STAT fetch failed
- `SRR28241792` / `NP270` / `UCSF_VRP_mNGS_074_1487_nohuman.fastq.gz`: STAT fetch failed

## E. Post-DTCA verification targets

- Paper post-DTCA OPA: **187 of 191 = 97.906%** (target 187/191 ≈ 97.9%)
- Paper post-DTCA PPA: 110.5 of 113 = 97.788%
- Paper post-DTCA NPA: 76.5 of 78 = 98.077%
- Residual FNs after DTCA (reconstruction): RSV 0.5 + metapneumovirus 1.0

