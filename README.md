# Toy GWAS Workshop

**Population genetics, statistical genetics, and GWAS — a hands-on tutorial**

A fully simulated case/control GWAS dataset (5 chromosomes, 25,000 SNPs,
~4,600 samples across a majority ancestry and minor ancestries) that
you QC, filter, and analyze end-to-end in R + PLINK2 — no real genetic data
involved, and no local software to install.

## Quick start

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/malmarri/popgen-gwas-tutorial)

Click the badge, wait for the container to build (a couple of minutes the
first time), then find RStudio: click the **PORTS** tab at the bottom of the
VS Code window, find port **8787**, and click its 🌐 globe icon to open it in
a new tab — no password required. (If it doesn't connect on the first try,
wait 30 seconds and click the globe icon again; the container needs a moment
to finish starting.) **Make sure you're actually in that RStudio tab** (URL
ending in `.app.github.dev`), not the Codespace's default VS Code/terminal
view — the notebooks need R, which only runs inside RStudio Server.

**Step 0, before anything else:** in RStudio's **Files** pane (bottom-right),
click **`popgen-gwas-tutorial.Rproj`** and confirm "Yes" to reopen as a
project. RStudio does not do this automatically, and every notebook needs
it — without it, R starts in your home folder instead of the project, and
`data/`, `results/`, and the notebooks' own file paths won't resolve. You'll
know it worked when the Console's top line shows the project name. If you
forget this step, the fix is one line in the Console:
`setwd("/workspaces/popgen-gwas-tutorial")`.

The toy dataset is generated automatically on first build
(`R/simulate_toy_gwas.R` runs via `postCreateCommand`); if you ever need to
regenerate it by hand, run this **in RStudio's Console** (bottom-left pane):

```r
source("R/simulate_toy_gwas.R")
```

If you're doing this from a plain terminal instead (e.g. the Codespace's
default shell, not RStudio's Console), use `Rscript` rather than `source`
— bash has its own unrelated `source` command that will not run R code:

```bash
Rscript R/simulate_toy_gwas.R
```

Either way, only run it from **one place at a time**. This script briefly
uses a few GB of memory; running it twice at once (e.g. once in a terminal
and once in RStudio's Console) roughly doubles that, and can get one of
the processes killed on a memory-constrained machine.

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
`popgen-gwas-tutorial-solutions`.) The `instructor/` truth tables aren't
generated at all unless an instructor explicitly asks for them
(`GWAS_WORKSHOP_TRUTH=1 Rscript R/simulate_toy_gwas.R`), so they won't be
sitting in your Codespace tempting you.

## Repository layout

```
.devcontainer/     RStudio + PLINK2 environment definition (Codespaces reads this)
R/                 simulation code (sim_utils.R, simulate_toy_gwas.R)
notebooks/         the six workshop notebooks (.Rmd)
data/              generated at runtime, gitignored — the PLINK filesets
results/           generated as you work through the notebooks, gitignored
instructor/        answer-key truth tables — only created on request, gitignored
```

Instructor answer key and build notes are kept in a separate private repo,
not published here.

## If something goes wrong (quick fixes)

**"Permission denied" writing to `data/` or `results/`** — open a Terminal
in RStudio (Tools > Terminal > New Terminal) and run:

```bash
chown -R rstudio:rstudio /workspaces/popgen-gwas-tutorial 2>/dev/null || sudo chown -R rstudio:rstudio /workspaces/popgen-gwas-tutorial
```

**`data/` is empty** — regenerate it. In RStudio's *Console* (bottom-left):

```r
source("R/simulate_toy_gwas.R")
```

Only run this in one place at a time — running it in a terminal and the
Console simultaneously doubles the memory it needs.

**"Can't find the project root"** — RStudio is pointed somewhere
unexpected. Either open `popgen-gwas-tutorial.Rproj`, or use
Session > Set Working Directory > To Project Directory, then re-run the
chunk.

**Notebook 05 says to run notebook 04 first** — it means it. The notebooks
build on each other's output files; run them 00 → 05 in order.

**You're in a terminal, not R** — if a prompt looks like
`root@codespaces-xxxx:/workspaces/...#` that's a *shell*, not R. There,
use `Rscript R/simulate_toy_gwas.R` (not R's `source()`, which is a
different, unrelated bash builtin). R code belongs in RStudio's Console.

## Requirements knowledge-wise

Basic command-line comfort and some R (reading a data frame, making a
ggplot) are assumed. No prior PLINK experience needed — every command is
explained inline.
