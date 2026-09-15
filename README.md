# Differential Expression Analysis: Skeletal Muscle Ageing

DESeq2-based transcriptomic analysis of human skeletal muscle across age groups, using publicly available RNA-seq data (GSE164471, Tumasian et al., *Nature Communications* 2021).

> **TL;DR** — DESeq2 differential expression comparing young (20–34 yrs) vs. old (65+ yrs) human skeletal muscle RNA-seq, with GO Biological Process enrichment on the resulting gene sets. Full results in [`results/`](./results); figures below.

## Motivation

Skeletal muscle ageing (sarcopenia) is a major driver of frailty, falls, and loss of independence in older adults, but the underlying transcriptional changes aren't fully mapped. Public RNA-seq datasets like GSE164471 make it possible to re-analyze this question independently, practice a complete differential expression workflow end to end, and generate a candidate gene/pathway list that could inform which targets are worth following up experimentally — a first step toward the kind of translational, biomarker-oriented work this project is meant to build toward.

## Hypothesis

Aged skeletal muscle will show a distinct transcriptional signature from young muscle, with differentially expressed genes enriched for processes already implicated in muscle ageing biology — e.g. mitochondrial/metabolic function, protein turnover, inflammation, and extracellular matrix remodeling. If age is a genuine biological signal rather than noise, young and old samples should also separate cleanly on PCA before any differential expression testing is done (Figure 1).

## Biological Question

Which genes and pathways are differentially expressed between young (20–34 years) and old (65+ years) human skeletal muscle?

## Dataset

| | |
|---|---|
| **Source** | GEO accession [GSE164471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164471) |
| **Tissue** | Human skeletal muscle biopsies |
| **Groups** | Young (20–34 yrs) vs. Old (65–79 yrs and 80+) |
| **Technology** | Total RNA-seq, Illumina HiSeq 2500 |

## Methods

- **Quality control:** DESeq2 variance-stabilising transformation + PCA
- **Differential expression:** DESeq2 (padj < 0.05, |log2FC| > 1)
- **Visualization:** Volcano plot (EnhancedVolcano), heatmap of top DEGs (pheatmap)
- **Pathway enrichment:** clusterProfiler GO Biological Process, run separately on the up- and downregulated gene sets

> GO enrichment plots are only generated when a given gene set (up or down) actually returns significant enriched terms — if `figures/04_GO_upregulated.png` or `05_GO_downregulated.png` is missing, that run found nothing significant for that direction, rather than the step having failed.

## Repository Structure

```
deseq2-aging-analysis/
├── deseq2-aging-analysis.Rmd   # full analysis, knit-able in RStudio
├── data/                       # input count matrices / metadata
├── figures/                    # exported plots (PCA, volcano, heatmap, GO)
├── results/                    # DE gene tables (significant + full)
└── README.md
```

## Key Figures

### Figure 1 — PCA: Sample Separation by Age Group
![PCA](https://github.com/VishCipher/deseq2-aging-analysis/raw/main/figures/01_pca_plot.png)

**What it is:** A principal component analysis (PCA) plot. Each point is one RNA-seq sample, positioned by its overall gene expression profile rather than any single gene.

**What it means:** If young and old samples form two separate clusters, age is a major driver of transcriptional variation in this dataset — the basic justification for running differential expression at all. A sample sitting far outside its group's cluster is a possible outlier worth a second look.

### Figure 2 — Volcano Plot: Differentially Expressed Genes
![Volcano](https://github.com/VishCipher/deseq2-aging-analysis/raw/main/figures/02_volcano_plot.png)

**What it is:** Every gene tested, plotted by effect size (x-axis: log2 fold change, old vs. young) against statistical confidence (y-axis: −log10 adjusted p-value).

**What it means:** Points high on the plot are the genes we're most confident actually changed; points far left or right changed the most. Genes in the upper-right increase with age, genes in the upper-left decrease — together, these are the significant DEGs reported in `results/significant_DEGs.csv`.

### Figure 3 — Heatmap: Top 40 Differentially Expressed Genes
![Heatmap](https://github.com/VishCipher/deseq2-aging-analysis/raw/main/figures/03_heatmap.png)

**What it is:** The top 40 DEGs by adjusted p-value (rows) across every individual sample (columns), colored by relative expression (row-scaled z-score, not raw counts).

**What it means:** This is a visual cross-check on the volcano plot — if the young and old samples split cleanly into two color blocks, it confirms these genes genuinely separate the two groups, rather than being driven by one or two unusual samples.

### Figure 4 — GO Enrichment: Pathways Upregulated in Aged Muscle
![GO Up](https://github.com/VishCipher/deseq2-aging-analysis/raw/main/figures/04_GO_upregulated.png)

**What it is:** GO Biological Process terms significantly enriched among the upregulated DEGs (clusterProfiler, BH-adjusted p < 0.05), shown as a dot plot.

**What it means:** Instead of looking at genes one at a time, this asks "what biological processes, as a group, are turning on with age?" — the terms shown are the ones the upregulated gene set is most strongly associated with.

*(Only present if this run's upregulated gene set returned significant enriched terms — see note above.)*

### Figure 5 — GO Enrichment: Pathways Downregulated in Aged Muscle
![GO Down](https://github.com/VishCipher/deseq2-aging-analysis/raw/main/figures/05_GO_downregulated.png)

**What it is:** The same enrichment analysis, restricted to the downregulated gene set.

**What it means:** The "what's turning off" counterpart to Figure 4. Read together, Figures 4 and 5 give the pathway-level summary of aged vs. young muscle, beyond the individual genes in Figures 2–3.

*(Only present if this run's downregulated gene set returned significant enriched terms — see note above.)*

## Results

- [`results/significant_DEGs.csv`](./results/significant_DEGs.csv) — genes passing padj < 0.05 and |log2FC| > 1
- [`results/full_DESeq2_results.csv`](./results/full_DESeq2_results.csv) — full DESeq2 output for every tested gene

## How to Reproduce

1. Clone this repo and open `deseq2-aging-analysis.Rmd` in RStudio.
2. Install the required packages (below).
3. Click **Knit** — the full analysis, from QC through GO enrichment, regenerates end to end.

**Required R packages:**

```r
BiocManager::install(c("DESeq2", "clusterProfiler",
                        "org.Hs.eg.db", "EnhancedVolcano"))
install.packages(c("tidyverse", "pheatmap"))
```

## Possible Extension

The `.Rmd` includes a commented-out section for cross-referencing DEGs against the [GenAge](https://genomics.senescence.info/genes/human.zip) human ageing-gene database. It's not run by default — enabling it requires downloading `genage_human.csv` separately and uncommenting that chunk.

## Reference

Tumasian RA 3rd, Harish A, Kundu G, Yang JH et al. Skeletal muscle transcriptome in healthy aging. *Nat Commun* 2021 Apr 1;12(1):2014. PMID: 33795677

## Author

**Vishishtaa Pandit** — [GitHub](https://github.com/VishCipher) · [LinkedIn](https://www.linkedin.com/in/vishishtaa-pandit28/)
