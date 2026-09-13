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
# "latest" is cog-genomics's stable convenience link for the Linux x86_64
# build. If this ever breaks, grab the current link from
# https://www.cog-genomics.org/plink/2.0/ and update the URL below.
RUN wget -q https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip -O /tmp/plink2.zip \
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
