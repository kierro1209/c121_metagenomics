#!/usr/bin/env python3
"""Harvest NCBI STAT tax_analysis rows for Accuracy SRRs from public S3 parquet."""
from __future__ import annotations
import json
import sys
from pathlib import Path
import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/normalized/accuracy_manifest.csv"
OUT_DIR = ROOT / "data/normalized/stat"
OUT_DIR.mkdir(parents=True, exist_ok=True)
FILE_LIST = Path("/tmp/stat_tax_files.txt")
PROGRESS = OUT_DIR / "harvest_progress.json"
HITS = OUT_DIR / "tax_analysis_hits.parquet"
SUMMARY = OUT_DIR / "tax_analysis_info_filtered.csv"

def main() -> None:
    manifest = pd.read_csv(MANIFEST)
    srrs = sorted(manifest["Run"].dropna().unique().tolist())
    srr_sql = ",".join(repr(s) for s in srrs)
    files = [ln.strip() for ln in FILE_LIST.read_text().splitlines() if ln.strip()]
    done = set()
    if PROGRESS.exists():
        done = set(json.loads(PROGRESS.read_text()).get("done_files", []))
    con = duckdb.connect(str(OUT_DIR / "harvest.duckdb"))
    con.execute("INSTALL httpfs; LOAD httpfs;")
    con.execute("SET s3_region='us-east-1';")
    con.execute("SET s3_access_key_id=''; SET s3_secret_access_key='';")
    con.execute(
        """
        CREATE TABLE IF NOT EXISTS tax_hits AS
        SELECT * FROM read_parquet(?)
        WHERE 1=0
        """,
        [
            f"s3://sra-pub-metadata-us-east-1/sra_tax_analysis_tool/tax_analysis/{files[0]}"
        ],
    )
    # recreate empty with schema from first file if empty
    n = con.execute("SELECT count(*) FROM tax_hits").fetchone()[0]
    if n == 0 and not done:
        con.execute("DROP TABLE tax_hits")
        con.execute(
            f"""
            CREATE TABLE tax_hits AS
            SELECT * FROM read_parquet('s3://sra-pub-metadata-us-east-1/sra_tax_analysis_tool/tax_analysis/{files[0]}')
            WHERE acc IN ({srr_sql})
            """
        )
        done.add(files[0])
        PROGRESS.write_text(json.dumps({"done_files": sorted(done), "n_files": len(files)}))
        print(f"seeded with {files[0]} rows={con.execute('select count(*) from tax_hits').fetchone()[0]}", flush=True)

    for i, fn in enumerate(files, 1):
        if fn in done:
            continue
        uri = f"s3://sra-pub-metadata-us-east-1/sra_tax_analysis_tool/tax_analysis/{fn}"
        try:
            con.execute(
                f"""
                INSERT INTO tax_hits
                SELECT * FROM read_parquet('{uri}')
                WHERE acc IN ({srr_sql})
                """
            )
        except Exception as e:
            print(f"ERROR {fn}: {e}", flush=True)
            PROGRESS.write_text(json.dumps({"done_files": sorted(done), "n_files": len(files), "last_error": str(e)}))
            raise
        done.add(fn)
        cnt = con.execute("SELECT count(*) FROM tax_hits").fetchone()[0]
        PROGRESS.write_text(json.dumps({"done_files": sorted(done), "n_files": len(files), "rows": cnt, "last_file": fn}))
        print(f"[{i}/{len(files)}] {fn} cumulative_rows={cnt}", flush=True)

    con.execute(f"COPY tax_hits TO '{HITS}' (FORMAT PARQUET)")
    # export viral-focused summary
    viral = con.execute(
        """
        SELECT * FROM tax_hits
        WHERE lower(name) LIKE '%virus%'
           OR lower(name) LIKE '%viridae%'
           OR lower(name) LIKE '%virinae%'
           OR lower(name) LIKE 'influenza%'
           OR lower(name) LIKE '%coronavirus%'
           OR lower(name) LIKE '%rhinovirus%'
           OR lower(name) LIKE '%enterovirus%'
           OR lower(name) LIKE '%adenovirus%'
           OR lower(name) LIKE '%metapneumovirus%'
           OR lower(name) LIKE '%bocavirus%'
           OR lower(name) LIKE '%respiratory syncytial%'
           OR lower(name) LIKE '%orthopneumovirus%'
           OR lower(name) LIKE '%respirovirus%'
           OR lower(name) LIKE '%rubulavirus%'
           OR lower(name) LIKE '%mastadenovirus%'
           OR lower(name) LIKE '%sarbecovirus%'
           OR lower(name) LIKE '%betacoronavirus%'
           OR lower(name) LIKE '%alphacoronavirus%'
           OR lower(name) LIKE '%pneumoviridae%'
           OR lower(name) LIKE '%paramyxoviridae%'
           OR lower(name) LIKE '%picornaviridae%'
           OR name IN ('Viruses','Orthomyxoviridae','Coronaviridae')
        ORDER BY acc, ileft, ilevel
        """
    ).fetchdf()
    viral.to_csv(OUT_DIR / "stat_viral_taxa.csv", index=False)
    print("wrote", HITS, "viral_rows", len(viral), flush=True)

if __name__ == "__main__":
    main()
