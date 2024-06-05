## CspStandSegmentation is an R-package for the segmentation of single trees from forest point clouds scanned with terrestrial, mobile or unmanned LiDAR systems <img src="https://github.com/JulFrey/CspStandSegmentation/blob/main/inst/figures/csp_logo.png" align="right" width = 300/>

Author: Julian Frey, University of Freiburg, Chair of Forest Growth and Dendroecology


## Installation

If you are working on Windows operating systems, you will need to install Rtools prior to installation: https://cran.r-project.org/bin/windows/Rtools/>. On Mac, Xcode is required. 

```R
install.packages(c('devtools', 'Rcpp', 'lidR', 'dbscan', 'igraph', 'foreach', 'parallel', 'doParallel','magrittr', 'data.table'))

devtools::install_github('https://github.com/JulFrey/CspStandSegmentation')

# Check if it is working
library(CspStandSegmentation)
example("csp_cost_segmentation")

```

## Usage
The package is firmly based on the `lidR` package and uses the las file structure. Smaller point clouds can be directly segmented using the ```csp_cost_segmentation``` function. This requires a set of tree positions (map) as starting points, which can be derived using the ```find_base_coordinates_raster``` function, which might require parameter optimization. Theoretically, tree positions might also come from field measurements or manual assignments.:

```R
# read example data
file = system.file("extdata", "pine_plot.laz", package="TreeLS")
tls = lidR::readTLSLAS(file)

# normalize height
tls <- TreeLS::tlsNormalize(tls)

# find tree positions as starting points for segmentation
map <- CspStandSegmentation::find_base_coordinates_raster(tls)

# segment trees
segmented <- tls |>
  CspStandSegmentation::add_geometry() |>
  CspStandSegmentation::csp_cost_segmentation(map, 1)

# show results
lidR::plot(segmented, color = "TreeID")
```

