# Toy GWAS Workshop

**Population genetics, statistical genetics, and GWAS — a hands-on tutorial**

A fully simulated case/control GWAS dataset (10 fake chromosomes, 50,000 SNPs,
~4,600 samples across a majority ancestry and two minority ancestries) that
you QC, filter, and analyze end-to-end in R + PLINK2 — no real genetic data
involved, and no local software to install.

## Quick start

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/malmarri/popgen-gwas-tutorial)

Click the badge, wait for the container to build (a couple of minutes the
first time), then open RStudio in the browser tab that appears — no password
required. The toy dataset is generated automatically on first build
(`R/simulate_toy_gwas.R` runs via `postCreateCommand`); if you ever need to
regenerate it by hand, run:

```r
source("R/simulate_toy_gwas.R")
```

Prefer to work locally? Any editor with the
[Dev Containers](https://containers.dev/) spec support (e.g. VS Code +
"Dev Containers" extension, with Docker installed) can open this same
`.devcontainer/` and get an identical environment.

## Agenda

Work through the notebooks in `notebooks/` in order — each one picks up
where the last left off:

| # | Notebook | Topic |
|---|---|---|
| 00 | `00_intro_and_setup.Rmd` | Orientation: what's in the dataset, sanity-check the simulation |
| 01 | `01_qc_missingness_heterozygosity.Rmd` | Sample QC: missingness & heterozygosity, choose your own thresholds |
| 02 | `02_pca_ancestry.Rmd` | PCA against a reference panel, identify and remove ancestry outliers |
| 03 | `03_maf_hwe_filtering.Rmd` | Variant QC: minor allele frequency & Hardy-Weinberg equilibrium filtering |
| 04 | `04_association_gwas.Rmd` | Running the association test with PLINK2, naive vs. ancestry-corrected |
| 05 | `05_interpretation_manhattan_qq.Rmd` | Manhattan/QQ plots, genomic inflation (λ), interpreting hits |

## What's deliberately hidden in the data

The dataset was built with specific, known problems and signals baked in, so
you can check your own work against ground truth once you've worked through
it honestly. The instructor's answer key is intentionally **not** in this
repository — it lives in a private repo instead, so students who peek at the
public source can't spoil it for themselves. (Instructors: see
`popgen-gwas-tutorial-solutions`.) The `instructor/` folder (truth tables)
is generated locally at runtime and gitignored for the same reason.

## Repository layout

```
.devcontainer/     RStudio + PLINK2 environment definition (Codespaces reads this)
R/                 simulation code (sim_utils.R, simulate_toy_gwas.R)
notebooks/         the six workshop notebooks (.Rmd)
data/              generated at runtime, gitignored — the PLINK filesets
results/           generated as you work through the notebooks, gitignored
instructor/        generated at runtime, gitignored — answer-key truth tables
```

Instructor answer key and build notes are kept in a separate private repo,
not published here.

## Requirements knowledge-wise

Basic command-line comfort and some R (reading a data frame, making a
ggplot) are assumed. No prior PLINK experience needed — every command is
explained inline.
