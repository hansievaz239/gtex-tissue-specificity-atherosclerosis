# GTEx Multi-Tissue Gene Expression: Tissue Specificity and Atherosclerosis Gene Panel

A bulk RNA-seq tissue-specificity analysis using public GTEx data, comparing
gene expression profiles across nine human tissues and examining where a
panel of atherosclerosis genetic-risk genes is actually expressed in the body.
Implemented in both R (`data.table`/`ggplot2`) and Python (pandas/matplotlib).

## Motivation

Genetic association studies identify disease-risk genes, but association
alone doesn't say where in the body a gene is biologically active. This
project connects back to an earlier project (atherosclerosis drug-target
landscape, via Open Targets) by asking a complementary question: of the
genes with strong genetic evidence for atherosclerosis risk, which are
actually expressed in arterial or lipid-metabolizing tissue specifically,
versus expressed broadly or in unrelated tissues? Tissue-specificity analysis
is a core technique in functional genomics for interpreting GWAS hits and
prioritizing targets.

## Data source

- **GTEx (Genotype-Tissue Expression project)**, v8 gene-level TPM data,
  downloaded per-tissue from the GTEx Portal (no login required for these
  bulk per-tissue files).
- Nine tissues, chosen to span distinct biology including the two organs
  central to lipid metabolism: **Artery - Aorta** (432 samples),
  **Artery - Coronary** (240 samples), **Brain - Hypothalamus**
  (202 samples), **Brain - Cerebellum** (241 samples), **Adrenal Gland**
  (258 samples), **Cervix - Endocervix** (10 samples), **Ovary**
  (180 samples), **Pancreas** (328 samples), **Liver** (226 samples).
- Each file contains per-sample TPM (transcripts per million) for 56,200
  genes; this analysis uses the per-gene **median** TPM across samples
  within each tissue as the tissue's representative expression value.

## Method

1. **Load and summarize**: for each tissue, computed the median TPM per gene
   across all samples in that tissue, producing a genes x 9-tissues matrix.
2. **Expression filter**: retained genes with median TPM > 1 in at least one
   of the nine tissues (23,043 / 56,200 genes).
3. **Tissue-tissue correlation**: computed Pearson correlation between
   tissues on log2(TPM+1) values across all filtered genes, to assess overall
   transcriptomic similarity.
4. **Tissue-specificity ranking**: for each tissue, ranked genes by
   log2 fold-change of that tissue's expression versus the mean of the other
   eight tissues, to identify the most tissue-enriched genes.
5. **Atherosclerosis gene panel**: pulled median expression across all nine
   tissues for a panel of genetic-risk genes (LDLR, APOE, PCSK9, HMGCR,
   SORT1, ABCG8, LPL, SH2B3, CETP, LPA) carried over from the earlier Open
   Targets atherosclerosis project, to see which tissues drive their
   biological activity.

## Files

- `step1_median_per_tissue.py`, `step2_analysis.py` - Python analysis scripts
- `gtex_tissue_analysis.R` - equivalent full pipeline in R
- `tissue_correlation_matrix.csv` / `tissue_correlation_heatmap.png` -
  pairwise tissue transcriptome correlation
- `atherosclerosis_gene_panel_by_tissue.csv` /
  `atherosclerosis_panel_heatmap.png` - the atherosclerosis-risk gene panel's
  expression across all nine tissues
- `top_tissue_specific_genes.csv` - top 20 most tissue-enriched genes per
  tissue, ranked by log2 fold-change vs. the other eight tissues

## Key findings

- **Artery - Aorta and Artery - Coronary remain each other's closest match**
  (Pearson r = 0.985 on log2 TPM), confirming a strong shared arterial
  transcriptomic signature. **Liver is the most transcriptomically distinct
  tissue in the panel**, showing its lowest correlation with Brain -
  Cerebellum (r = 0.56) - consistent with liver's highly specialized
  metabolic gene expression program.
- **Top tissue-specific genes recover well-established tissue biology**: MGP,
  ITGA8, PLN, MYH11, ACTA2, BGN, and COL8A1 (vascular smooth-muscle
  structural/contractile genes) top the two artery lists; GFAP, AVP, OXT,
  SST, and MBP top hypothalamus; PVALB and GABRD top cerebellum; CYP11B1,
  CYP17A1, and STAR (steroidogenesis enzymes) top adrenal gland; PRSS1/PRSS2
  and CPA1 (digestive enzymes) top pancreas; and ALB, FGB, and
  APOC3/APOA2 (plasma proteins and apolipoproteins) top liver - a strong
  sanity check that the pipeline correctly recovers known biology across very
  different organ systems.
- **PCSK9 is essentially undetectable in all eight non-liver tissues but
  clearly expressed in liver** (median TPM ~25.5 vs. <3 everywhere else),
  directly confirming PCSK9's well-known liver-dominant expression.
- **ABCG8 and LPA show the same pattern**: both essentially zero across all
  eight other tissues but clearly expressed in liver (ABCG8 ~48 TPM, LPA
  ~15 TPM), consistent with their known roles in hepatic biliary cholesterol
  excretion (ABCG8) and hepatic apolipoprotein(a)/Lp(a) synthesis (LPA).
- **SORT1 remains clearly elevated in both artery tissues** relative to the
  other seven tissues, independently supporting the arterial relevance of
  this well-established atherosclerosis/LDL-cholesterol GWAS locus, and is
  notably *low* in liver despite SORT1's established role in hepatic LDL
  metabolism in the literature - a reminder that bulk tissue median
  expression doesn't always track mechanistic importance.
- **APOE is dramatically highest in liver and adrenal gland** (median TPM
  >3,000 and >3,500 respectively), consistent with both liver's central role
  in lipoprotein production and adrenal steroidogenic tissue's reliance on
  lipoprotein-derived cholesterol.

## Limitations

- Only nine of GTEx's ~54 tissues are included; a full tissue-specificity
  analysis (as in the GTEx publications) would use all available tissues for
  a complete picture.
- Analysis uses the median across samples per tissue as a single
  representative value, discarding inter-individual variability within each
  tissue; a more complete analysis would model that variability directly
  (e.g., via differential expression testing).
- Correlation and tissue-specificity were computed only on genes passing a
  simple TPM > 1 filter in at least one tissue; results for lowly-expressed
  genes near this threshold should be interpreted cautiously.
- Bulk tissue expression reflects the average across all cell types in a
  tissue sample; genes with cell-type-restricted but functionally important
  expression (e.g., SORT1 in hepatocytes specifically) can appear lower than
  expected when diluted across a whole-tissue sample.

## Author

Hansie Vaz - PhD Researcher, Computational Biology, Dublin City University
