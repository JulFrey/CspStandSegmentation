FROM rocker/geospatial:4.5

# # Install system dependencies
# RUN apt-get update && apt-get install -y \
# build-essential \
# libcurl4-openssl-dev \
# libssl-dev \
# libxml2-dev \
# libudunits2-dev \
# libgdal-dev \
# libgeos-dev \
# libproj-dev \
# libglpk-dev \
# libpng-dev \
# libtiff-dev \
# libgmp-dev \
# libfontconfig1-dev \
# git \
# curl \
# && apt-get clean

# Install CRAN packages and GitHub package using 'remotes'
RUN R -e "install.packages(c( \
    'remotes', 'Rcpp', 'dbscan', 'igraph', 'foreach', \
    'doParallel', 'magrittr', 'data.table', 'plumber', 'sf', 'stars', 'lidR' \
    ))"

# Install custom segmentation package from GitHub
RUN R -e "remotes::install_github('JulFrey/CspStandSegmentation')"

# Create working directory
WORKDIR /app

# Copy plumber script
COPY R/2025-06-03_plumber_service.R /app/api.R

# Expose port
EXPOSE 8000

# Start plumber
CMD ["R", "-e", "pr <- plumber::plumb('/app/api.R'); pr$run(host='0.0.0.0', port=8000)"]
