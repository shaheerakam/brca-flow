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

## Part 7: Validation by positive controls (60-tumor test run)

Check: five genes with a known direction in breast cancer should differ in the expected direction between ER-positive and ER-negative tumors. ESR1, GATA3 and FOXA1 should be higher in ER-positive tumors. KRT5 and KRT17 should be lower.

Result: all five matched the expected direction (5 of 5 PASS), shown below from results/tables/positive_controls.tsv.

What this supports: the ER-positive vs ER-negative labels and the contrast direction are the right way round, and the counts, filtering and PyDESeq2 steps are behaving sensibly. A flipped comparison or mislabeled groups would have failed most of these.

Caveats:
- ER status was approximated from PAM50 subtype. I believe several of these genes (for example ESR1, FOXA1, KRT5 and KRT17) are among the genes used to define PAM50 subtypes, so the check is partly circular and is weaker than an independent test. I still need to confirm this against the PAM50 gene list.
- This is a 60-tumor test run. The check needs repeating on the larger or full dataset.
- A direction match does not show that the effect sizes or p-values are accurate. That is why I will also compare against UALCAN, GEPIA2 and Kaplan-Meier Plotter (not done yet).

Positive control table:

gene	expected	log2FC	padj	result
ESR1	1	5.13	9.62078233436037e-34	PASS
GATA3	1	2.94	4.361028419746399e-14	PASS
FOXA1	1	2.16	1.3547016754235446e-06	PASS
KRT5	-1	-3.1	5.203413166325104e-05	PASS
KRT17	-1	-1.26	0.08483221444381164	PASS

## Part 7: Validation by positive controls (60-tumor test run)

Check: five genes with a known direction in breast cancer should differ in the expected direction between ER-positive and ER-negative tumors. ESR1, GATA3 and FOXA1 should be higher in ER-positive tumors. KRT5 and KRT17 should be lower.

Result: all five matched the expected direction (5 of 5 PASS), shown below from results/tables/positive_controls.tsv.

What this supports: the ER-positive vs ER-negative labels and the contrast direction are the right way round, and the counts, filtering and PyDESeq2 steps are behaving sensibly. A flipped comparison or mislabeled groups would have failed most of these.

Caveats:
- ER status was approximated from PAM50 subtype. I believe several of these genes (for example ESR1, FOXA1, KRT5 and KRT17) are among the genes used to define PAM50 subtypes, so the check is partly circular and is weaker than an independent test. I still need to confirm this against the PAM50 gene list.
- This is a 60-tumor test run. The check needs repeating on the larger or full dataset.
- A direction match does not show that the effect sizes or p-values are accurate. That is why I will also compare against UALCAN, GEPIA2 and Kaplan-Meier Plotter (not done yet).

Positive control table:

gene	expected	log2FC	padj	result
ESR1	1	5.13	9.62078233436037e-34	PASS
GATA3	1	2.94	4.361028419746399e-14	PASS
FOXA1	1	2.16	1.3547016754235446e-06	PASS
KRT5	-1	-3.1	5.203413166325104e-05	PASS
KRT17	-1	-1.26	0.08483221444381164	PASS

### Volcano plot (results/figures/volcano_er.png)
Of 16509 genes tested in the 60-tumor run, 1567 were higher in ER-positive tumors and 1736 were higher in ER-negative tumors, using padj below 0.05 and an absolute log2 fold change above 1 (red dots). The plot has the expected two-wing shape, with significant genes on both sides. Gray dots near the center are genes with small changes or weak evidence.

The most significant genes reach a -log10 adjusted p-value of about 45. The largest fold changes extend further on the ER-positive side (up to about +12) than on the ER-negative side (about -7.5). Some of the most extreme fold changes may come from genes with low counts, where estimates are less stable. I have not checked which genes these are.

Caveats: the number of significant genes is large, which is partly because ER status was approximated from PAM50 subtype, a label defined from expression, so many genes track it. With only 60 tumors, this is a test run and the counts will change on the larger dataset. Thresholds are the usual conventions and not tuned.

Top genes by adjusted p-value:

gene_id	log2FoldChange	padj	gene_name
ENSG00000173467.9	6.748380546474649	4.500217830314659e-47	AGR3
ENSG00000134830.6	4.352635221121162	2.1910789562612763e-45	C5AR2
ENSG00000166535.20	-7.55105138326054	2.6009052125441374e-39	A2ML1
ENSG00000120262.10	4.086931976670416	4.0069974715730915e-36	CCDC170
ENSG00000082175.15	6.769998340633373	9.39212878004479e-34	PGR

## Part 8: Survival analysis (60-tumor test run)

### Kaplan-Meier curve (results/figures/km_gene.png)
Gene: GATA3, chosen before looking at survival results, not selected from the data. Patients were split at the median GATA3 expression into a high group (n=30) and a low group (n=29). The log-rank test gave p = 0.577, so there is no significant difference in overall survival between the groups in this subset.

The curves cross. The GATA3-high group drops earlier (steps at roughly 5, 10 and 33 months) and then stays flat at about 0.85. The GATA3-low group stays higher until about 56 months and then drops to about 0.69 by roughly 75 months. The low-expression curve ends near 107 months, and the high-expression curve extends to about 170 months. Because the curves cross and the test is not significant, I do not read this as evidence that GATA3 level affects survival.

### Caveats
- Very few deaths: each curve falls in only a handful of steps, so one patient moving changes the curve a lot. The later drop in the low group rests on few patients still at risk.
- This is a 60-tumor test run. The result needs repeating on the larger or full dataset.
- ER status was approximated from PAM50 subtype, and I did not adjust for ER status, stage or age, which all relate to survival.
- Exploratory analysis. A non-significant result is still a result, and I did not try other genes to find a significant one.

### Numbers from survival_summary.txt
Gene: GATA3 (ENSG00000107485.18)
Patients analysed: 59
Deaths (events): 6
Log-rank p-value: 0.5765
