# Article discrepancy reconstruction

Scope: Supplementary Tables 6 and 7 plus Supplementary Data 1 joins only.
No FASTQ classification and no final inconsistency report.

## Sources
- Article: https://www.nature.com/articles/s41467-024-51470-y
- Supplementary Information (Tables 6-7): https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM1_ESM.pdf
- Supplementary Data 1: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM4_ESM.xlsx

## Exact counts
- mNGS+/RVP- samples (Table 6): 14
- mNGS-/RVP+ main-text equivalents: 7
- Table 7 published sample rows: 9 (title n=9; TN+FN weight 7.0 = TN 5.5 + FN 1.5; match=True)
- Pre-DTCA sensitivity misses: 7 = 110-103
- Post-DTCA OPA: 187 of 191 = 97.906%

## Outputs
- `supplementary_table6_discrepancies.csv`
- `supplementary_table7_discrepancies.csv`
- `article_discrepancies_sample_level.csv`
- `article_discrepancies_organism_level.csv`
- `article_discrepancies_maintext7_tn_fn.csv`
- `article_discrepancies_unresolved_ids.csv`
- `../results/intermediate/discrepancy_case_cards.json`
- `../results/intermediate/discrepancy_reconstruction_summary.json`

## Unresolved issues preserved
- Most Table 6 RVP-negative rows lack sample IDs in Supplementary Data 1.
- Table 6 row 3 vs P-NP-36 RSV call conflicts between Table 6 and Supp Data 1.
- Table 7 title n=9 vs main-text n=7: reconciled as TN+FN fractional weights 5.5+1.5=7.
- P-NP-99/P-NP-100 cannot be ordered against the two rhinovirus Table 7 rows.
- Influenza A + RSV / mNGS ND Table 7 row has no Supp Data 1 match.
- NP/NUC↔P-NP lab-ID crosswalks are unpublished; not invented here (see `unresolved_identifier_mappings.csv`).
