# Accuracy sequence audit (STAT)

Scope: independent sequence-content profiling of all 281 Purpose=Accuracy runs.
Uses NCBI STAT public `run_taxonomy` as a proxy for the article ≥3 non-overlapping
reads/contigs evidence rule. Does not download FASTQs and does not invent
NP/NUC↔P-NP crosswalks.

## Script

- `analysis/audit_accuracy_sequences.py`
- Cache: `results/intermediate/stat_cache/` (gitignored)
- Outputs:
  - `results/intermediate/accuracy_sequence_audit.csv`
  - `results/intermediate/accuracy_sequence_audit.json`
  - `results/intermediate/accuracy_stat_virus_hits.csv`
  - `results/intermediate/accuracy_sequence_audit_summary.json`

## Final report bundle

Built by `analysis/build_accuracy_audit_report.py`:

- `results/accuracy_metadata_audit.csv`
- `results/inconsistent_samples.md`
- `results/metadata_consistency_report.md`
