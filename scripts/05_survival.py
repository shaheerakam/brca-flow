"""Script 5: Kaplan-Meier survival curves for one gene (high vs low expression).

Usage:
    python scripts/05_survival.py               # uses the top gene from script 4
    python scripts/05_survival.py --gene GATA3  # or pick a gene yourself

Writes: results/figures/km_gene.png and results/tables/survival_summary.txt
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from brca_flow.expression import log2_cpm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", default="data/processed/counts.tsv.gz")
    parser.add_argument("--sample-sheet", default="data/processed/sample_sheet.tsv")
    parser.add_argument("--de-table",
                        default="results/tables/de_er_positive_vs_negative.tsv")
    parser.add_argument("--gene", default=None, help="gene symbol, e.g. GATA3")
    parser.add_argument("--fig-dir", default="results/figures")
    parser.add_argument("--table-dir", default="results/tables")
    args = parser.parse_args()

    fig_dir, table_dir = Path(args.fig_dir), Path(args.table_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    de = pd.read_csv(args.de_table, sep="\t", index_col=0).dropna(subset=["padj"])
    if args.gene:
        matches = de[de["gene_name"] == args.gene]
        if matches.empty:
            raise SystemExit(f"Gene {args.gene!r} was not found in {args.de_table}")
        gene_id = matches["baseMean"].idxmax()
    else:
        gene_id = de.index[0]  # table is sorted by adjusted p-value
    gene_name = de.loc[gene_id, "gene_name"]

    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    sheet = pd.read_csv(args.sample_sheet, sep="\t").set_index("patient_id")
    sheet = sheet.dropna(subset=["os_months", "os_event"])
    sheet = sheet[sheet["os_months"] > 0]
    expression = log2_cpm(counts)[sheet.index].loc[gene_id]

    high = expression >= expression.median()
    kmf_high, kmf_low = KaplanMeierFitter(), KaplanMeierFitter()
    fig, ax = plt.subplots(figsize=(6.5, 5))
    kmf_high.fit(sheet.loc[high, "os_months"], sheet.loc[high, "os_event"],
                 label=f"{gene_name} high (n={int(high.sum())})")
    kmf_low.fit(sheet.loc[~high, "os_months"], sheet.loc[~high, "os_event"],
                label=f"{gene_name} low (n={int((~high).sum())})")
    kmf_high.plot_survival_function(ax=ax, ci_show=False, color="#c0392b")
    kmf_low.plot_survival_function(ax=ax, ci_show=False, color="#2c7fb8")

    test = logrank_test(sheet.loc[high, "os_months"], sheet.loc[~high, "os_months"],
                        sheet.loc[high, "os_event"], sheet.loc[~high, "os_event"])
    ax.set_xlabel("Months since diagnosis")
    ax.set_ylabel("Overall survival probability")
    ax.set_title(f"{gene_name}: log-rank p = {test.p_value:.3g}")
    fig.tight_layout()
    fig.savefig(fig_dir / "km_gene.png", dpi=200)
    plt.close(fig)

    summary = (f"Gene: {gene_name} ({gene_id})\n"
               f"Patients analysed: {len(sheet)}\n"
               f"Deaths (events): {int(sheet['os_event'].sum())}\n"
               f"Log-rank p-value: {test.p_value:.4g}\n")
    (table_dir / "survival_summary.txt").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
