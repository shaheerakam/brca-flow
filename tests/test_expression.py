"""Tests for the expression helper functions. Run them with:  pytest"""
import numpy as np
import pandas as pd
import pytest

from brca_flow.expression import filter_low_counts, log2_cpm


def test_filter_low_counts_drops_quiet_genes():
    counts = pd.DataFrame({"s1": [100, 0], "s2": [120, 1], "s3": [90, 0]},
                          index=["busy_gene", "quiet_gene"])
    kept = filter_low_counts(counts, min_count=10, min_samples=2)
    assert list(kept.index) == ["busy_gene"]


def test_log2_cpm_is_independent_of_library_size():
    # Sample s2 has exactly twice the reads of s1, so after scaling
    # to counts-per-million the two samples must look identical.
    counts = pd.DataFrame({"s1": [500, 500], "s2": [1000, 1000]},
                          index=["gene_a", "gene_b"])
    logged = log2_cpm(counts)
    assert logged.loc["gene_a", "s1"] == pytest.approx(logged.loc["gene_a", "s2"])
    # Each gene is half the library = 500,000 per million.
    assert logged.loc["gene_a", "s1"] == pytest.approx(np.log2(500_001))
