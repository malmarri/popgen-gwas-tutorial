# R/plot_utils.R
#
# Shared helpers for reading PLINK2 --glm output and building Manhattan/QQ
# plots by hand in ggplot2, so the mechanics stay visible to students rather
# than hidden inside a single-purpose plotting package.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

#' Read a plink2 --glm logistic output file into a tidy data frame.
read_glm <- function(path) {
  df <- data.table::fread(path)
  # plink2's column names vary a little by version/build
  # (e.g. "#CHROM" vs "CHROM", "P" vs "LOG10_P") -- normalise defensively.
  names(df) <- sub("^#", "", names(df))
  if (!"P" %in% names(df) && "LOG10_P" %in% names(df)) {
    df$P <- 10^(-df$LOG10_P)
  }
  df <- as.data.frame(df)
  df[!is.na(df$P) & df$P > 0, ]
}

#' Genomic inflation factor, lambda_GC: the median chi-square test statistic
#' (back-calculated from the p-values) divided by the value expected under
#' the null with no true association anywhere (qchisq(0.5, df = 1)).
#' lambda_GC ~ 1 means "the bulk of the genome behaves as expected"; well
#' above 1 is a classic sign of uncorrected population stratification.
lambda_gc <- function(pvals) {
  chisq <- qchisq(1 - pvals, df = 1)
  median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
}

#' Manhattan plot, built directly so the mechanics (cumulative position
#' across chromosomes, alternating colour, -log10(p)) stay visible.
plot_manhattan <- function(df, title = "", sig_line = NULL) {
  chr_lengths <- df %>%
    group_by(CHROM) %>%
    summarise(chr_len = max(POS), .groups = "drop") %>%
    arrange(CHROM)
  chr_lengths$offset <- cumsum(c(0, head(chr_lengths$chr_len, -1)))

  df <- df %>%
    left_join(chr_lengths %>% select(CHROM, offset), by = "CHROM") %>%
    mutate(cum_pos = POS + offset)

  axis_df <- df %>%
    group_by(CHROM) %>%
    summarise(center = mean(cum_pos), .groups = "drop") %>%
    arrange(CHROM)

  p <- ggplot(df, aes(x = cum_pos, y = -log10(P), color = factor(CHROM %% 2))) +
    geom_point(size = 0.6, alpha = 0.7, show.legend = FALSE) +
    scale_color_manual(values = c("0" = "#4C72B0", "1" = "#8C8C8C")) +
    scale_x_continuous(labels = axis_df$CHROM, breaks = axis_df$center) +
    labs(title = title, x = "Chromosome", y = expression(-log[10](p))) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

  if (!is.null(sig_line)) {
    p <- p + geom_hline(yintercept = -log10(sig_line), linetype = "dashed", color = "firebrick")
  }
  p
}

#' QQ plot of observed vs. expected -log10(p) under the null.
plot_qq <- function(pvals, title = "") {
  pvals <- sort(pvals[pvals > 0])
  n <- length(pvals)
  expected <- -log10(ppoints(n))
  observed <- -log10(pvals)
  ggplot(data.frame(expected, observed), aes(expected, observed)) +
    geom_point(size = 0.6, alpha = 0.5, color = "#4C72B0") +
    geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = "dashed") +
    labs(title = title,
         x = expression(Expected ~ -log[10](p)),
         y = expression(Observed ~ -log[10](p))) +
    theme_minimal(base_size = 12)
}
