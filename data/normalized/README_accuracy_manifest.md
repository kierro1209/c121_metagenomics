# Accuracy manifest provenance

Scope: metadata assembly only. No sequence audit or final adjudication was run.

Selection: `Purpose == "Accuracy"` in the Zenodo sequencing master sheet.

## Source URLs

- Article: https://www.nature.com/articles/s41467-024-51470-y
- Zenodo record: https://zenodo.org/records/10553379
- Zenodo API record: https://zenodo.org/api/records/10553379
- Sequencing record master sheet; Purpose filter and run identifiers: https://zenodo.org/api/records/10553379/files/respmNGS-sequencing-mastersheet.xlsx/content
- Clinical and laboratory metadata (Supplementary Data 1): https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM4_ESM.xlsx
- Supplementary Tables 1-8, including discrepancy Tables 6-7: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM1_ESM.pdf
- Figure source data; contextual, not complete clinical metadata: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51470-y/MediaObjects/41467_2024_51470_MOESM6_ESM.xlsx
- SRA/BioSample accessions and read metadata for sequence retrieval: https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1084017

## Verified counts

- Accuracy sequencing records: 281
- Clinical sequencing records: 233
- Control sequencing records: 48
- Observed clinical specimen keys (before unresolved alias collapse): 206
- Article clinical cohort: 191 (186 unambiguous keys + 5 unresolved members of repeat run 039)
- Rerun/replicate records beyond each specimen's first record: 27
- Unresolved NP/NUC-to-article identifier crosswalks: 79

See `accuracy_manifest_provenance.json` for hashes and machine-readable details.
