"""Script 4: find genes that differ between ER-positive and ER-negative tumors.

Uses PyDESeq2, a Python version of the widely used DESeq2 method.

Writes: results/tables/de_er_positive_vs_negative.tsv
        results/tables/positive_controls.tsv
        results/figures/volcano_er.png
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

from brca_flow.expression import filter_low_counts

# Genes whose direction we already know from published breast cancer biology.
# Expected sign of log2 fold change for ER_positive vs ER_negative.
POSITIVE_CONTROLS = {"ESR1": +1, "GATA3": +1, "FOXA1": +1, "KRT5": -1, "KRT17": -1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", default="data/processed/counts.tsv.gz")
    parser.add_argument("--genes", default="data/processed/genes.tsv")
    parser.add_argument("--sample-sheet", default="data/processed/sample_sheet.tsv")
    parser.add_argument("--fig-dir", default="results/figures")
    parser.add_argument("--table-dir", default="results/tables")
    args = parser.parse_args()

    fig_dir, table_dir = Path(args.fig_dir), Path(args.table_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    genes = pd.read_csv(args.genes, sep="\t").set_index("gene_id")
    sheet = pd.read_csv(args.sample_sheet, sep="\t").set_index("patient_id")
    sheet = sheet[sheet["er_status"].isin(["ER_positive", "ER_negative"])]
    counts = filter_low_counts(counts[sheet.index], min_count=10, min_samples=10)
    print(f"Testing {counts.shape[0]} genes in {counts.shape[1]} tumors")

    # PyDESeq2 wants samples as rows and genes as columns.
    dds = DeseqDataSet(counts=counts.T, metadata=sheet[["er_status"]],
                       design="~er_status", quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["er_status", "ER_positive", "ER_negative"],
                       quiet=True)
    stats.summary()
    results = stats.results_df.join(genes).sort_values("padj")
    results.index.name = "gene_id"
    results.to_csv(table_dir / "de_er_positive_vs_negative.tsv", sep="\t")

    # Sanity check against known biology.
    rows = []
    for name, expected in POSITIVE_CONTROLS.items():
        hit = results[results["gene_name"] == name]
        if hit.empty:
            rows.append({"gene": name, "expected": expected, "log2FC": np.nan,
                         "padj": np.nan, "result": "NOT TESTED"})
            continue
        best = hit.iloc[0]
        ok = np.sign(best["log2FoldChange"]) == expected
        rows.append({"gene": name, "expected": expected,
                     "log2FC": round(best["log2FoldChange"], 2),
                     "padj": best["padj"], "result": "PASS" if ok else "FAIL"})
    controls = pd.DataFrame(rows)
    controls.to_csv(table_dir / "positive_controls.tsv", sep="\t", index=False)
    print(controls.to_string(index=False))

    # Volcano plot
    plot = results.dropna(subset=["padj"])
    significant = (plot["padj"] < 0.05) & (plot["log2FoldChange"].abs() > 1)
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(plot["log2FoldChange"], -np.log10(plot["padj"]), s=6,
               color=np.where(significant, "#c0392b", "#b0b0b0"))
    ax.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=0.6)
    ax.axvline(-1, linestyle="--", color="black", linewidth=0.6)
    ax.axvline(1, linestyle="--", color="black", linewidth=0.6)
    ax.set_xlabel("log2 fold change (ER positive vs ER negative)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Differential expression by ER status")
    fig.tight_layout()
    fig.savefig(fig_dir / "volcano_er.png", dpi=200)
    plt.close(fig)
    print(f"{int(significant.sum())} genes with padj < 0.05 and |log2FC| > 1")


if __name__ == "__main__":
    main()
