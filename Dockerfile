# Start from the official CRAN Debian test environment
FROM cran/debian:latest

# Install required system libraries for spatial and compiled packages
RUN apt-get update && apt-get install -y \
    cmake \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libsqlite3-dev \
    libtiff5-dev \
    libglu1-mesa-dev \
    mesa-common-dev \
    libssl-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libabsl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# Install essential helper R packages
RUN R -e "install.packages(c('remotes','rcmdcheck'), repos='https://cloud.r-project.org')"

# Set working directory
WORKDIR /pkg

# Install R dependencies defined in the package DESCRIPTION
# This happens *after* your package source is mounted into /pkg
CMD Rscript -e "\
  remotes::install_deps('/pkg', dependencies = TRUE, repos='https://cloud.r-project.org'); \
  rcmdcheck::rcmdcheck('/pkg', args='--no-manual', error_on='note')"
