"""Script 1: download TCGA-BRCA gene count files from the GDC.

Usage:
    python scripts/01_download_expression.py --limit 60      # small test run
    python scripts/01_download_expression.py                 # everything
"""
import argparse
import json
import time
from pathlib import Path

import pandas as pd
import requests

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"
GDC_STATUS = "https://api.gdc.cancer.gov/status"


def build_filters():
    """Describe which files we want: open-access TCGA-BRCA STAR gene counts."""

    def clause(field, value):
        return {"op": "in", "content": {"field": field, "value": [value]}}

    return {
        "op": "and",
        "content": [
            clause("cases.project.project_id", "TCGA-BRCA"),
            clause("data_category", "Transcriptome Profiling"),
            clause("data_type", "Gene Expression Quantification"),
            clause("analysis.workflow_type", "STAR - Counts"),
            clause("access", "open"),
        ],
    }


def list_files():
    """Ask the GDC for the list of matching files and return it as a table."""
    body = {
        "filters": build_filters(),
        "fields": ("file_id,file_name,cases.submitter_id,"
                   "cases.samples.submitter_id,cases.samples.sample_type"),
        "format": "JSON",
        "size": 5000,
    }
    response = requests.post(GDC_FILES, json=body, timeout=120)
    response.raise_for_status()
    rows = []
    for hit in response.json()["data"]["hits"]:
        case = hit["cases"][0]
        sample = case["samples"][0]
        rows.append({
            "file_id": hit["file_id"],
            "file_name": hit["file_name"],
            "patient_id": case["submitter_id"],
            "sample_barcode": sample["submitter_id"],
            "sample_type": sample["sample_type"],
        })
    return pd.DataFrame(rows)


def download_one(file_id, destination, tries=3):
    """Download a single file, retrying a few times if the network hiccups."""
    for attempt in range(1, tries + 1):
        try:
            with requests.get(f"{GDC_DATA}/{file_id}", stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(destination, "wb") as handle:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        handle.write(chunk)
            return
        except requests.RequestException as error:
            print(f"  attempt {attempt} failed for {file_id}: {error}")
            time.sleep(2 * attempt)
    raise RuntimeError(f"Could not download {file_id} after {tries} tries")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=None,
                        help="only download this many samples (random, fixed seed)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--out-dir", default="data/raw")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    expression_dir = out_dir / "expression"
    expression_dir.mkdir(parents=True, exist_ok=True)

    status = requests.get(GDC_STATUS, timeout=60).json()
    (out_dir / "gdc_status.json").write_text(json.dumps(status, indent=2))
    print("GDC data release:", status.get("data_release"))

    files = list_files()
    files = files[files["sample_type"] == "Primary Tumor"]
    files = files.sort_values("sample_barcode").drop_duplicates("patient_id")
    print(f"Found {len(files)} primary tumor samples (one per patient).")

    if args.limit is not None and args.limit < len(files):
        files = files.sample(n=args.limit, random_state=args.seed)
        print(f"Keeping a random subset of {len(files)} samples.")
    files = files.sort_values("sample_barcode").reset_index(drop=True)

    for number, row in enumerate(files.itertuples(), start=1):
        destination = expression_dir / f"{row.file_id}.tsv"
        if destination.exists() and destination.stat().st_size > 0:
            continue
        print(f"[{number}/{len(files)}] downloading {row.patient_id}")
        download_one(row.file_id, destination)
        time.sleep(0.2)

    # The manifest is written last, so its existence means downloading finished.
    files.to_csv(out_dir / "manifest.tsv", sep="\t", index=False)
    print("Done. Wrote", out_dir / "manifest.tsv")


if __name__ == "__main__":
    main()
