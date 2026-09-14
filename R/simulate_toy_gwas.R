# R/simulate_toy_gwas.R
#
# Generates the toy GWAS dataset used throughout the workshop notebooks.
# Everything is simulated -- there is no real genetic data involved, and
# the "chromosomes" are fake (5 of them, 5,000 SNPs each = 25,000 total).
#
# Run once per workshop:  Rscript R/simulate_toy_gwas.R
# (this also happens automatically when the devcontainer is built)
#
# Outputs:
#   data/study.{bed,bim,fam}                 raw study data, PRE-QC
#   data/reference.{bed,bim,fam}             reference panel, 4 populations
#
# Instructors only -- written ONLY when GWAS_WORKSHOP_TRUTH=1 is set, so that
# the answer key doesn't appear in every student's container (see below):
#   instructor/snp_truth.csv                 which SNPs are causal/decoy
#   instructor/sample_truth.csv              true ancestry + injected QC problems per sample

suppressPackageStartupMessages(library(genio))
source(file.path("R", "sim_utils.R"))

# Fixed seed -> every student gets an identical dataset. Chosen (out of
# ~20 candidates tested against real PLINK2 output after the 10->5
# chromosome change) specifically because it gives a clean, unambiguous
# single false positive at decoy_1 in the naive analysis with no other
# SNP spuriously crossing threshold, AND keeps causal_6 (the deliberately
# low-frequency/large-effect SNP) clearly detectable after correction --
# both properties took real searching to get simultaneously at 25,000
# SNPs (half the test count means less "room" for both to line up).
set.seed(7)

## ---- Study design constants ----------------------------------------------
N_CHR        <- 5
SNPS_PER_CHR <- 5000
N_SNPS       <- N_CHR * SNPS_PER_CHR  # 25,000

N_CASES_TARGET    <- 2000
N_CONTROLS_TARGET <- 2000

N_BAD_MAJORITY <- 130  # extra majority-ancestry samples that should FAIL QC
N_ANC_B        <- 250  # minority ancestry, moderate divergence
N_ANC_C        <- 250  # minority ancestry, larger divergence

FST_A <- 0.005
FST_B <- 0.01
FST_C <- 0.02

N_REF_PER_POP <- 80    # reference panel size per population
FST_REF_D     <- 0.20  # a 4th reference population NOT present in the study

dir.create("data", showWarnings = FALSE)

# The instructor truth tables name every causal SNP and every planted QC
# problem outright -- i.e. the entire answer key. They are NOT written by
# default, because this script runs automatically in each student's
# container and an "instructor/" folder sitting in the file browser is a
# spoiler one click away. Instructors: regenerate them with
#   GWAS_WORKSHOP_TRUTH=1 Rscript R/simulate_toy_gwas.R
WRITE_TRUTH <- Sys.getenv("GWAS_WORKSHOP_TRUTH") == "1"
if (WRITE_TRUTH) dir.create("instructor", showWarnings = FALSE)

## ---- 1. SNP map + arbitrary allele labels ---------------------------------
snp_map <- build_snp_map(N_CHR, SNPS_PER_CHR)
alleles <- assign_alleles(N_SNPS)
bim <- data.frame(
  chr  = snp_map$chr,
  id   = snp_map$id,
  posg = 0,
  pos  = snp_map$pos,
  alt  = alleles$alt,
  ref  = alleles$ref,
  stringsAsFactors = FALSE
)

## ---- 2. Pick the 6 causal SNPs + 1 decoy SNP -------------------------------
# With only 5 chromosomes, several special roles now share a chromosome
# (distinguished by position instead): chr1/chr2 hold three causal SNPs
# each, chr3 holds the decoy alone (it's the centerpiece of notebook 04,
# worth keeping easy to point at), chr4 holds both QC-artifact blocks,
# and chr5 is left completely clean as the null-chromosome contrast.
special_idx <- function(chrom, offset) which(bim$chr == chrom)[offset]

idx_causal <- c(
  causal_1 = special_idx(1, 1250), causal_2 = special_idx(1, 2500), causal_3 = special_idx(1, 3750),
  causal_4 = special_idx(2, 1250), causal_5 = special_idx(2, 2500), causal_6 = special_idx(2, 3750)
)
idx_decoy <- c(decoy_1 = special_idx(3, 2500))

causal_maf <- c(causal_1 = 0.30, causal_2 = 0.20, causal_3 = 0.15,
                causal_4 = 0.10, causal_5 = 0.05, causal_6 = 0.02)
causal_or  <- c(causal_1 = 1.4,  causal_2 = 1.4,  causal_3 = 1.35,
                causal_4 = 1.4,  causal_5 = 1.4,  causal_6 = 2.3)
log_or_causal <- log(causal_or)

decoy_freq_A <- 0.05  # ancestry A and B share this frequency (no confounding there)
decoy_freq_C <- 0.95  # ancestry C differs sharply -> drives the stratification example

# chr4, positions 500-800: block of SNPs with a genuine HWE violation
# (genotyping artefact). Same block size (300 SNPs) as before.
idx_hwe_break <- which(bim$chr == 4)[500:800]
# chr4, positions 3500-3700: block of SNPs with high missingness (failing
# assay). Same block size (200 SNPs) as before; well clear of the HWE
# block above so the two artefacts stay distinguishable by position.
idx_snp_miss  <- which(bim$chr == 4)[3500:3700]
# chr5 is left with no injected features at all -- a "clean null" chromosome

## ---- 3. Population allele frequencies -------------------------------------
p_anc <- runif(N_SNPS, 0.05, 0.50)  # ancestral/global frequency for ordinary SNPs

p_A <- bn_pop_freq(p_anc, FST_A)
p_B <- bn_pop_freq(p_anc, FST_B)
p_C <- bn_pop_freq(p_anc, FST_C)

# Causal SNPs: force IDENTICAL frequency across ancestries. Their signal
# should come only from the true genotype effect, not from stratification.
for (nm in names(idx_causal)) {
  i <- idx_causal[[nm]]
  p_A[i] <- causal_maf[[nm]]
  p_B[i] <- causal_maf[[nm]]
  p_C[i] <- causal_maf[[nm]]
}

# Decoy SNP: force a real frequency DIFFERENCE between A/B and C, no true effect.
p_A[[idx_decoy]] <- decoy_freq_A
p_B[[idx_decoy]] <- decoy_freq_A
p_C[[idx_decoy]] <- decoy_freq_C

## ---- 4. Majority ancestry (Ancestry A) case/control pool ------------------
# Retrospective case-control ascertainment: simulate genotypes + a liability
# score from the 6 causal SNPs, keep sampling batches until we have exactly
# N_CASES_TARGET cases and N_CONTROLS_TARGET controls.
# Pre-allocate the final-size matrices and fill columns in place as batches
# come in, rather than accumulating a list of chunks and concatenating
# afterward -- the list-then-cbind-then-subset approach briefly holds 2-3
# full-size (~400MB) copies of the same data at once, which is what was
# pushing this script's peak memory high enough to get OOM-killed on a
# memory-constrained Codespace machine. This produces byte-identical output
# to the old approach (same random draws, same order, just a different
# place to put them), because only the *storage* changed, not the RNG calls.
geno_cases    <- matrix(NA_integer_, nrow = N_SNPS, ncol = N_CASES_TARGET)
geno_controls <- matrix(NA_integer_, nrow = N_SNPS, ncol = N_CONTROLS_TARGET)
n_case_have <- 0; n_control_have <- 0
batch <- 1
while (n_case_have < N_CASES_TARGET || n_control_have < N_CONTROLS_TARGET) {
  batch_n <- 1500
  g <- sim_genotypes(batch_n, p_A)

  liability <- rep(0, batch_n)  # logit(0.5) baseline
  for (nm in names(idx_causal)) {
    liability <- liability + log_or_causal[[nm]] * g[idx_causal[[nm]], ]
  }
  is_case <- rbinom(batch_n, 1, plogis(liability)) == 1

  cases_here    <- which(is_case)
  controls_here <- which(!is_case)

  if (n_case_have < N_CASES_TARGET && length(cases_here) > 0) {
    take <- head(cases_here, N_CASES_TARGET - n_case_have)
    geno_cases[, (n_case_have + 1):(n_case_have + length(take))] <- g[, take, drop = FALSE]
    n_case_have <- n_case_have + length(take)
  }
  if (n_control_have < N_CONTROLS_TARGET && length(controls_here) > 0) {
    take <- head(controls_here, N_CONTROLS_TARGET - n_control_have)
    geno_controls[, (n_control_have + 1):(n_control_have + length(take))] <- g[, take, drop = FALSE]
    n_control_have <- n_control_have + length(take)
  }
  message(sprintf("  batch %d: cases %d/%d, controls %d/%d",
                   batch, n_case_have, N_CASES_TARGET, n_control_have, N_CONTROLS_TARGET))
  batch <- batch + 1
  if (batch > 50) stop("Could not reach target case/control counts -- check effect sizes.")
  # Without this, each batch's ~300MB `g` piles up unreclaimed rather than
  # being freed before the next batch allocates -- R doesn't trigger a
  # garbage collection on every reassignment, so across several batches
  # this was the single largest source of peak memory in the whole script.
  rm(g); invisible(gc(FALSE))
}

n_majority_clean     <- N_CASES_TARGET + N_CONTROLS_TARGET
pheno_majority_clean <- c(rep(2L, N_CASES_TARGET), rep(1L, N_CONTROLS_TARGET))  # 2=case, 1=control

## ---- 5. Extra majority-ancestry samples designed to FAIL QC ---------------
g_bad <- sim_genotypes(N_BAD_MAJORITY, p_A)
pheno_bad <- sample(c(1L, 2L), N_BAD_MAJORITY, replace = TRUE)

het_bad_idx  <- seq_len(N_BAD_MAJORITY %/% 2)
miss_bad_idx <- (N_BAD_MAJORITY %/% 2 + 1):N_BAD_MAJORITY
g_bad <- corrupt_heterozygosity(g_bad, het_bad_idx)
g_bad <- corrupt_sample_missingness(g_bad, miss_bad_idx)

## ---- 6. Minority ancestries B and C, with uneven case/control ascertainment
# This ascertainment imbalance, combined with the decoy SNP's real frequency
# gap between ancestries, is what creates a population-stratification false
# positive at decoy_1 in the naive (uncorrected) association test.
g_ancB <- sim_genotypes(N_ANC_B, p_B)
g_ancC <- sim_genotypes(N_ANC_C, p_C)
pheno_ancB <- sample(c(1L, 2L), N_ANC_B, replace = TRUE, prob = c(0.85, 0.15))  # mostly controls
pheno_ancC <- sample(c(1L, 2L), N_ANC_C, replace = TRUE, prob = c(0.15, 0.85))  # mostly cases

## ---- 7. Assemble the full pre-QC "study" fileset --------------------------
# Single cbind (cases, controls, bad, ancB, ancC) instead of stacking through
# an intermediate geno_majority_clean copy -- same column order and values,
# one fewer ~1GB full-matrix duplicate held in memory at once.
geno_study <- cbind(geno_cases, geno_controls, g_bad, g_ancB, g_ancC)
rm(geno_cases, geno_controls, g_bad, g_ancB, g_ancC); invisible(gc(FALSE))

# Give each causal SNP realistic-looking "LD shoulders" -- nearby SNPs on
# the same chromosome become partially correlated with it, correlation
# decaying with distance, so the Manhattan plot shows the skyscraper-with-
# shoulders shape of a real locus instead of a single isolated spike.
# Deliberately NOT applied to decoy_1: its whole lesson is that it looks
# associated with no real biology behind it, and giving it LD support
# would undermine that by making it look like a genuine locus.
for (nm in names(idx_causal)) {
  geno_study <- add_ld_block(geno_study, idx_causal[[nm]])
}

geno_study <- corrupt_hwe(geno_study, idx_hwe_break)
geno_study <- corrupt_snp_missingness(geno_study, idx_snp_miss)
invisible(gc(FALSE))

n_study <- ncol(geno_study)
fam_study <- data.frame(
  fam   = sprintf("F%05d", seq_len(n_study)),
  id    = sprintf("S%05d", seq_len(n_study)),
  pat   = 0, mat = 0,
  sex   = sample(1:2, n_study, replace = TRUE),
  pheno = c(pheno_majority_clean, pheno_bad, pheno_ancB, pheno_ancC),
  stringsAsFactors = FALSE
)

write_plink(file.path("data", "study"), X = geno_study, bim = bim, fam = fam_study)
rm(geno_study); invisible(gc(FALSE))

## ---- 8. Reference panel (separate individuals, identical SNP map) --------
p_refD <- bn_pop_freq(p_anc, FST_REF_D)
g_refA <- sim_genotypes(N_REF_PER_POP, p_A)
g_refB <- sim_genotypes(N_REF_PER_POP, p_B)
g_refC <- sim_genotypes(N_REF_PER_POP, p_C)
g_refD <- sim_genotypes(N_REF_PER_POP, p_refD)

geno_ref <- cbind(g_refA, g_refB, g_refC, g_refD)
rm(g_refA, g_refB, g_refC, g_refD); invisible(gc(FALSE))
n_ref <- ncol(geno_ref)
fam_ref <- data.frame(
  fam   = sprintf("R%05d", seq_len(n_ref)),
  id    = sprintf("R%05d", seq_len(n_ref)),
  pat   = 0, mat = 0,
  sex   = sample(1:2, n_ref, replace = TRUE),
  pheno = -9,
  stringsAsFactors = FALSE
)
write_plink(file.path("data", "reference"), X = geno_ref, bim = bim, fam = fam_ref)

## ---- 9. Instructor-only truth tables (answer key -- do not share) --------
if (WRITE_TRUTH) {
snp_truth <- data.frame(id = bim$id, role = "neutral", stringsAsFactors = FALSE)
snp_truth$role[idx_causal] <- names(idx_causal)
snp_truth$role[idx_decoy]  <- names(idx_decoy)
snp_truth$target_maf <- NA_real_
snp_truth$target_or  <- NA_real_
snp_truth$target_maf[idx_causal] <- causal_maf
snp_truth$target_or[idx_causal]  <- causal_or
snp_truth$target_or[idx_decoy]   <- 1.0
write.csv(snp_truth, file.path("instructor", "snp_truth.csv"), row.names = FALSE)

sample_truth <- data.frame(
  fam = fam_study$fam, id = fam_study$id,
  true_ancestry = c(
    rep("Ancestry_A", n_majority_clean + N_BAD_MAJORITY),
    rep("Ancestry_B", N_ANC_B),
    rep("Ancestry_C", N_ANC_C)
  ),
  qc_flag = c(
    rep("clean", n_majority_clean),
    rep("het_outlier", length(het_bad_idx)),
    rep("high_missing", length(miss_bad_idx)),
    rep("clean_but_ancestry_outlier", N_ANC_B + N_ANC_C)
  ),
  stringsAsFactors = FALSE
)
write.csv(sample_truth, file.path("instructor", "sample_truth.csv"), row.names = FALSE)
}

message("Done.")
message("  data/study.{bed,bim,fam}       -- ", n_study, " samples, ", N_SNPS, " SNPs (pre-QC)")
message("  data/reference.{bed,bim,fam}   -- ", n_ref, " samples (4 reference populations)")
if (WRITE_TRUTH) {
  message("  instructor/*_truth.csv         -- answer key, keep out of student view")
} else {
  message("  (instructor truth tables not written; set GWAS_WORKSHOP_TRUTH=1 for them)")
}
