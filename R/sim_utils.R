# R/sim_utils.R
#
# Helper functions used by simulate_toy_gwas.R to build the toy dataset.
# Nothing here touches real genetic data -- every "chromosome", position,
# and allele frequency is fabricated for teaching purposes.

#' Draw a population-specific allele frequency from an ancestral frequency
#' using the Balding-Nichols model (the standard way to fake realistic
#' population structure without a full coalescent simulator).
#'
#' @param p_anc vector of ancestral (global) allele frequencies
#' @param fst   population differentiation (F_ST) relative to the ancestral
#'              population -- bigger fst = more diverged population
bn_pop_freq <- function(p_anc, fst) {
  a <- p_anc * (1 - fst) / fst
  b <- (1 - p_anc) * (1 - fst) / fst
  p <- rbeta(length(p_anc), a, b)
  pmin(pmax(p, 0.001), 0.999)  # keep away from the 0/1 boundary
}

#' Simulate genotypes under Hardy-Weinberg equilibrium (random mating) for a
#' set of individuals, given one allele frequency per SNP.
#'
#' @param n_ind  number of individuals to simulate
#' @param p_vec  per-SNP allele frequency (length = number of SNPs)
#' @return integer matrix, SNPs x individuals (the orientation genio's
#'   write_plink() expects), values in {0, 1, 2}
sim_genotypes <- function(n_ind, p_vec) {
  # Looping per-SNP (instead of building one giant n_ind*n_snps probability
  # vector) keeps peak memory manageable at 50,000 SNPs x thousands of people.
  geno <- vapply(
    p_vec,
    function(p) rbinom(n_ind, size = 2, prob = p),
    FUN.VALUE = integer(n_ind)
  )
  t(geno)  # vapply gives n_ind x n_snps; transpose to n_snps x n_ind
}

#' Build a fake multi-chromosome SNP map (chromosome + position + ID only;
#' alleles are added separately by assign_alleles()).
build_snp_map <- function(n_chr, snps_per_chr, spacing_bp = 1000) {
  chr <- rep(seq_len(n_chr), each = snps_per_chr)
  pos <- rep(seq_len(snps_per_chr) * spacing_bp, times = n_chr)
  id  <- sprintf("chr%d_snp%05d", chr, rep(seq_len(snps_per_chr), times = n_chr))
  data.frame(chr = chr, id = id, pos = pos, stringsAsFactors = FALSE)
}

#' Assign an arbitrary biallelic REF/ALT letter pair to each SNP. These
#' labels carry no meaning of their own here -- ref vs. alt is not the same
#' thing as major vs. minor or risk vs. non-risk allele, and the notebooks
#' deliberately do not assume otherwise (see the MAF notebook for why that
#' distinction matters when reading real PLINK output too).
assign_alleles <- function(n_snps) {
  bases <- c("A", "C", "G", "T")
  ref <- sample(bases, n_snps, replace = TRUE)
  alt <- vapply(ref, function(r) sample(setdiff(bases, r), 1), character(1))
  data.frame(ref = ref, alt = alt, stringsAsFactors = FALSE, row.names = NULL)
}

#' Simulate sample contamination: for each targeted individual, flip a
#' fraction of their homozygous genotypes to heterozygous, inflating their
#' apparent heterozygosity rate.
corrupt_heterozygosity <- function(geno, ind_idx, flip_rate = 0.35) {
  for (j in ind_idx) {
    hom_idx <- which(geno[, j] %in% c(0L, 2L))
    n_flip <- round(length(hom_idx) * flip_rate)
    if (n_flip > 0) {
      flip_idx <- sample(hom_idx, n_flip)
      geno[flip_idx, j] <- 1L
    }
  }
  geno
}

#' Simulate poor-quality DNA/genotyping for a set of individuals by setting
#' a fraction of their genotype calls to missing (NA), genome-wide.
corrupt_sample_missingness <- function(geno, ind_idx, miss_rate = 0.15) {
  n_snps <- nrow(geno)
  for (j in ind_idx) {
    n_miss <- round(n_snps * miss_rate)
    miss_idx <- sample(n_snps, n_miss)
    geno[miss_idx, j] <- NA
  }
  geno
}

#' Simulate a failing assay/probe for a set of SNPs by setting a fraction of
#' calls to missing across all individuals.
corrupt_snp_missingness <- function(geno, snp_idx, miss_rate = 0.20) {
  n_ind <- ncol(geno)
  for (i in snp_idx) {
    n_miss <- round(n_ind * miss_rate)
    miss_idx <- sample(n_ind, n_miss)
    geno[i, miss_idx] <- NA
  }
  geno
}

#' Break Hardy-Weinberg equilibrium at a set of SNPs by depleting
#' heterozygotes (mimics a genotyping/batch artefact -- a real and common
#' cause of HWE failure that has nothing to do with population structure).
corrupt_hwe <- function(geno, snp_idx, dropout_rate = 0.5) {
  for (i in snp_idx) {
    het_idx <- which(geno[i, ] == 1L)
    n_drop <- round(length(het_idx) * dropout_rate)
    if (n_drop > 1) {
      drop_idx <- sample(het_idx, n_drop)
      to_zero <- sample(drop_idx, length(drop_idx) %/% 2)
      to_two  <- setdiff(drop_idx, to_zero)
      geno[i, to_zero] <- 0L
      geno[i, to_two]  <- 2L
    }
  }
  geno
}

#' Merge two PLINK filesets that were generated from the SAME underlying SNP
#' map (same variants, same ref/alt coding) by stacking their individuals,
#' restricted to whichever SNP IDs both filesets still have in common.
#'
#' This is a shortcut, not a general-purpose merge: real datasets from two
#' different sources almost never share a common marker set this cleanly,
#' which is why real pipelines use plink's own merge machinery (and have to
#' handle strand flips, allele mismatches, etc.). Here the study data and
#' the reference panel were both generated from the exact same simulated
#' SNP map, so after one side has been QC-filtered down to a subset of that
#' map, intersecting by SNP ID and stacking individuals in R is simpler and
#' has zero moving parts to debug live in a workshop.
merge_same_snps <- function(prefix1, prefix2, out_prefix) {
  d1 <- genio::read_plink(prefix1)
  d2 <- genio::read_plink(prefix2)

  common_ids <- intersect(d1$bim$id, d2$bim$id)
  if (length(common_ids) == 0) stop("No SNP IDs in common between the two filesets.")
  message(sprintf("Merging on %d SNPs common to both filesets (%d in '%s', %d in '%s').",
                   length(common_ids), nrow(d1$bim), prefix1, nrow(d2$bim), prefix2))

  i1 <- match(common_ids, d1$bim$id)
  i2 <- match(common_ids, d2$bim$id)
  stopifnot(identical(d1$bim$ref[i1], d2$bim$ref[i2]),
            identical(d1$bim$alt[i1], d2$bim$alt[i2]))

  geno <- cbind(d1$X[i1, , drop = FALSE], d2$X[i2, , drop = FALSE])
  fam  <- rbind(d1$fam, d2$fam)
  genio::write_plink(out_prefix, X = geno, bim = d1$bim[i1, ], fam = fam)
  invisible(NULL)
}
