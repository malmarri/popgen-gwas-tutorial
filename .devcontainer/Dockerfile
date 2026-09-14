FROM rocker/rstudio:4.4.1

# --- System dependencies -----------------------------------------------
# unzip/wget: to fetch PLINK2
# lib*-dev: needed to compile a couple of the R packages below
RUN apt-get update && apt-get install -y --no-install-recommends \
    unzip \
    wget \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# --- PLINK 2.0 ------------------------------------------------------------
# PINNED to a specific dated build rather than the floating "_latest" link,
# so every student gets the identical binary no matter which day they create
# their Codespace, and so a new upstream release can't silently change flag
# behaviour or output formats mid-workshop. This exact build (v2.0.0-a.7.0LM,
# 25 Apr 2026) was verified to still accept every flag this workshop uses --
# notably `--glm ... allow-no-covars` -- and to still compute --hardy/--hwe
# over *founders* (which notebook 03 explains). To move to a newer build,
# get the dated link from https://www.cog-genomics.org/plink/2.0/ , update
# the URL below, and re-run the notebooks end to end before the workshop.
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20260425.zip -O /tmp/plink2.zip \
    && unzip -o /tmp/plink2.zip -d /usr/local/bin \
    && rm /tmp/plink2.zip \
    && chmod +x /usr/local/bin/plink2 \
    && plink2 --version

# --- R packages -------------------------------------------------------
RUN install2.r --error --skipinstalled \
    genio \
    dplyr \
    tidyr \
    ggplot2 \
    readr \
    data.table \
    rmarkdown \
    knitr
