"""Script 3: quality control - library sizes and a PCA plot.

Writes: results/figures/qc_library_sizes.png
        results/figures/pca_er_status.png
        results/tables/qc_summary.txt
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # draw to files, not to a window
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA

from brca_flow.expression import filter_low_counts, log2_cpm


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", default="data/processed/counts.tsv.gz")
    parser.add_argument("--sample-sheet", default="data/processed/sample_sheet.tsv")
    parser.add_argument("--fig-dir", default="results/figures")
    parser.add_argument("--table-dir", default="results/tables")
    parser.add_argument("--min-samples", type=int, default=10)
    args = parser.parse_args()

    fig_dir, table_dir = Path(args.fig_dir), Path(args.table_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    sheet = pd.read_csv(args.sample_sheet, sep="\t").set_index("patient_id")
    sheet = sheet[sheet["er_status"].notna()]
    counts = counts[sheet.index]

    # Figure 1: how many reads does each sample have?
    library_sizes = counts.sum(axis=0) / 1e6
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(library_sizes, bins=30, color="#3b6ea5")
    ax.set_xlabel("Total reads per sample (millions)")
    ax.set_ylabel("Number of samples")
    ax.set_title("Library sizes")
    fig.tight_layout()
    fig.savefig(fig_dir / "qc_library_sizes.png", dpi=200)
    plt.close(fig)

    # Figure 2: PCA on the 2,000 most variable genes
    filtered = filter_low_counts(counts, min_count=10, min_samples=args.min_samples)
    logged = log2_cpm(filtered)
    top_genes = logged.var(axis=1).nlargest(2000).index
    matrix = logged.loc[top_genes].T
    matrix = matrix - matrix.mean(axis=0)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(matrix)

    fig, ax = plt.subplots(figsize=(6, 5))
    colors = {"ER_positive": "#c0392b", "ER_negative": "#2c7fb8"}
    for group, color in colors.items():
        mask = (sheet["er_status"] == group).values
        ax.scatter(coords[mask, 0], coords[mask, 1], s=18, alpha=0.8,
                   color=color, label=group)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.0%} of variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.0%} of variance)")
    ax.set_title("PCA of tumor expression, colored by ER status")
    ax.legend()
    fig.tight_layout()
    fig.savefig(fig_dir / "pca_er_status.png", dpi=200)
    plt.close(fig)

    summary = [
        f"Samples with known ER status: {counts.shape[1]}",
        f"ER_positive: {(sheet['er_status'] == 'ER_positive').sum()}",
        f"ER_negative: {(sheet['er_status'] == 'ER_negative').sum()}",
        f"Genes before filtering: {counts.shape[0]}",
        f"Genes after filtering:  {filtered.shape[0]}",
        f"Smallest library (millions of reads): {library_sizes.min():.1f}",
        f"Largest library (millions of reads):  {library_sizes.max():.1f}",
    ]
    (table_dir / "qc_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
