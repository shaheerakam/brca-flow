"""Script 2: combine expression files and clinical data into clean tables.

Reads:   data/raw/manifest.tsv and the downloaded expression files
Fetches: clinical data from the public cBioPortal API (then saves a copy)
Writes:  data/processed/counts.tsv.gz   (genes x patients, raw counts)
         data/processed/genes.tsv       (gene_id -> gene_name)
         data/processed/sample_sheet.tsv (one row per patient)
"""
import argparse
import sys
from pathlib import Path

import pandas as pd
import requests

from brca_flow.expression import read_star_counts
from brca_flow.tcga import (er_from_subtype, find_er_column,
                            normalize_er_status, parse_os_status)

CBIOPORTAL_API = "https://www.cbioportal.org/api"


def fetch_clinical_table(study_id, data_type, cache_path):
    """Get clinical data from cBioPortal as a wide table (one row per patient).

    data_type is 'PATIENT' or 'SAMPLE'. The result is cached on disk so the
    internet is only needed the first time.
    """
    if cache_path.exists():
        print("Using saved copy:", cache_path)
        return pd.read_csv(cache_path, sep="\t")

    url = f"{CBIOPORTAL_API}/studies/{study_id}/clinical-data"
    params = {"clinicalDataType": data_type, "projection": "SUMMARY"}
    response = requests.get(url, params=params, timeout=300,
                            headers={"Accept": "application/json"})
    response.raise_for_status()
    long_table = pd.DataFrame(response.json())
    wide = long_table.pivot_table(index="patientId", columns="clinicalAttributeId",
                                  values="value", aggfunc="first").reset_index()
    wide.to_csv(cache_path, sep="\t", index=False)
    return wide


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="data/raw/manifest.tsv")
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--out-dir", default="data/processed")
    parser.add_argument("--study-id", default="brca_tcga_pan_can_atlas_2018")
    parser.add_argument("--er-column", default=None,
                        help="name of the clinical column holding ER status")
    parser.add_argument("--er-from-subtype", action="store_true",
                        help="approximate ER status from the PAM50 SUBTYPE column")
    args = parser.parse_args()

    raw_dir, out_dir = Path(args.raw_dir), Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t")

    # ---- Part A: expression counts -------------------------------------
    series, genes = [], None
    for row in manifest.itertuples():
        table = read_star_counts(raw_dir / "expression" / f"{row.file_id}.tsv")
        table = table[table["gene_type"] == "protein_coding"].set_index("gene_id")
        if genes is None:
            genes = table[["gene_name"]]
        series.append(table["unstranded"].rename(row.patient_id))
    counts = pd.concat(series, axis=1)
    counts.to_csv(out_dir / "counts.tsv.gz", sep="\t")
    genes.reset_index().to_csv(out_dir / "genes.tsv", sep="\t", index=False)
    print(f"Count table: {counts.shape[0]} genes x {counts.shape[1]} patients")

    # ---- Part B: clinical data -----------------------------------------
    patient = fetch_clinical_table(args.study_id, "PATIENT",
                                   raw_dir / "clinical_patient.tsv")
    sample = fetch_clinical_table(args.study_id, "SAMPLE",
                                  raw_dir / "clinical_sample.tsv")
    sample = sample.drop_duplicates("patientId")
    clinical = patient.merge(sample, on="patientId", how="outer",
                             suffixes=("", "_sample"))

    # ---- Part C: the sample sheet --------------------------------------
    sheet = manifest[["patient_id", "sample_barcode", "file_id"]].merge(
        clinical, left_on="patient_id", right_on="patientId", how="left")

    if args.er_from_subtype:
        sheet["er_status"] = sheet["SUBTYPE"].map(er_from_subtype)
        print("WARNING: ER status was approximated from PAM50 subtype.")
    else:
        er_column = args.er_column or find_er_column(list(clinical.columns))
        if er_column is None or er_column not in clinical.columns:
            print("Could not find an ER status column. Columns available:")
            for name in sorted(clinical.columns):
                print("  ", name)
            print("Re-run with --er-column NAME, or use --er-from-subtype.")
            sys.exit(2)
        print("Using ER column:", er_column)
        sheet["er_status"] = sheet[er_column].map(normalize_er_status)

    sheet["os_months"] = pd.to_numeric(sheet["OS_MONTHS"], errors="coerce")
    sheet["os_event"] = sheet["OS_STATUS"].map(parse_os_status)
    keep = ["patient_id", "sample_barcode", "file_id", "er_status",
            "os_months", "os_event"]
    sheet[keep].to_csv(out_dir / "sample_sheet.tsv", sep="\t", index=False)

    print(sheet["er_status"].value_counts(dropna=False).to_string())
    has_survival = sheet["os_months"].notna() & sheet["os_event"].notna()
    print("Patients with survival data:", int(has_survival.sum()))


if __name__ == "__main__":
    main()
