## CspStandSegmentation ia a R-package for the segmenattion of single trees from forest point clouds scanned with terestrial, mobile or unmanned LiDAR systems

Author: Julian Frey, University of Freiburg, Chair of Forest Growth and Dendroecology

## installation
```R
install.packages('Rcpp', 'lidR', 'TreeLS', 'dbscan', 'igraph', 'foreach', 'parallel', 'doParallel','magrittr', 'data.table', 'future.apply')
devtools::install_github('https://github.com/JulFrey/CspStandSegmentation')
```

## Usage
The package is strongly based on the TreeLS and lidR packages and uses the las file structure. Smaller point clouds can be directly segmented using the ```csp_cost_segmentation``` function. This requires a set of tree positions (map) as starting points, which can be derived using ```TreeLS::treeMap``` function, which might require some parameter optimization. Theoretically tree positions might also come from field measurements or manual assignments.:

```R
# read example data
file = system.file("extdata", "pine_plot.laz", package="TreeLS")
tls = lidR::readTLSLAS(file)

# normalize height
tls <- TreeLS::tlsNormalize(tls)

# find tree positions as starting point for segmentation
map <- TreeLS::treeMap(tls)

# segment trees
segmented <- tls |>
  lidR::filter_poi(Classification != 2) |>
  add_geometry() |>
  csp_cost_segmentation(map, 1)

lidR::plot(segmented, color = "TreeID")
```

