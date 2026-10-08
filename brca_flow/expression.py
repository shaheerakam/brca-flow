"""Helpers for reading and transforming gene expression counts."""
import numpy as np
import pandas as pd


def read_star_counts(path):
    """Read one GDC 'STAR - Counts' file into a table.

    The file starts with a comment line, then a header, then four summary
    rows whose gene_id starts with 'N_'. We drop the comment and summary rows.
    """
    table = pd.read_csv(path, sep="\t", comment="#")
    return table[~table["gene_id"].str.startswith("N_")]


def filter_low_counts(counts, min_count=10, min_samples=10):
    """Keep genes with at least min_count reads in at least min_samples samples.

    counts is a table with genes as rows and samples as columns.
    """
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def log2_cpm(counts, pseudocount=1.0):
    """Convert raw counts to log2 counts-per-million (genes x samples)."""
    library_sizes = counts.sum(axis=0)
    cpm = counts.div(library_sizes, axis=1) * 1_000_000
    return np.log2(cpm + pseudocount)
