# Analysis Notes

## Part 6: Quality control (60-tumor test run)

### Library sizes (results/figures/qc_library_sizes.png)
Total reads per sample ranged from about 28 to 91 million. No sample had very few reads, so I saw no clear low-quality outliers and did not remove any samples. The histogram has several humps, which may reflect samples sequenced in different batches. I did not test for batch effects, so I note it as a limitation.

### PCA (results/figures/pca_er_status.png)
PC1 explains about 20% of the variance and PC2 about 9%. ER-negative tumors fall mostly on the left of PC1 and ER-positive tumors mostly on the right, so the groups partly separate. A few tumors sit in the middle and overlap, so ER status is not a perfect divider in this data.

### Caveat
ER status here was approximated from PAM50 subtype (--er-from-subtype), and PAM50 subtypes are defined from gene expression. The separation is therefore partly expected and is not independent confirmation that the pipeline is correct.

### Numbers from qc_summary.txt
Samples with known ER status: 56
ER_positive: 40
ER_negative: 16
Genes before filtering: 19962
Genes after filtering:  16509
Smallest library (millions of reads): 27.4
Largest library (millions of reads):  91.3
