# STAT sequence audit artifacts

Independent taxonomy audit of all 281 Accuracy SRRs via NCBI STAT
`run_taxonomy` API (`https://www.ncbi.nlm.nih.gov/Traces/sra-db-be/run_taxonomy`).

## Rule
- Threshold: `total_count >= 3` as a proxy for the paper’s ≥3 non-overlapping reads/contigs.
- Evidence type: sequence-derived STAT k-mer/spot counts (not SURPI+ coverage).

## Key files
- `accuracy_stat_audit_detail.csv` — full per-run STAT + classification
- `accuracy_stat_audit_summary.csv` — compact summary
- `accuracy_stat_audit_provenance.json` — counts and limitations
- `accuracy_stat_fetch_errors.csv` — SRRs without STAT payloads
- `tax_analysis_info_accuracy.csv` — S3 tax_analysis_info rows for Accuracy SRRs
- `run_taxonomy_cache/` — raw API JSON (gitignored)

## Related
- Scaffold for extended FASTQ confirmation: `analysis/audit_accuracy_fastq.py`
- Final merged deliverables: `results/accuracy_metadata_audit.csv`, `results/inconsistent_samples.md`, `results/metadata_consistency_report.md`
