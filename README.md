## CspStandSegmentation is an R-package for the segmentation of single trees from forest point clouds scanned with terrestrial, mobile or unmanned LiDAR systems <img src="https://github.com/JulFrey/CspStandSegmentation/blob/main/inst/figures/csp_logo.png" align="right" width = 300/>

Authors: Julian Frey and Zoe Schindler, University of Freiburg, Chair of Forest Growth and Dendroecology


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17294732.svg)](https://doi.org/10.5281/zenodo.17294732)  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![R-CMD-check](https://github.com/JulFrey/CspStandSegmentation/actions/workflows/r.yml/badge.svg)](https://github.com/JulFrey/CspStandSegmentation/actions/workflows/r.yml)




## Installation

If you are working on Windows operating systems, you will need to install Rtools prior to installation: https://cran.r-project.org/bin/windows/Rtools/>. On Mac, Xcode is required. 

```R
install.packages(c('devtools', 'Rcpp', 'lidR', 'dbscan', 'igraph', 'foreach', 'doParallel','magrittr', 'data.table'))

devtools::install_github('https://github.com/JulFrey/CspStandSegmentation')

# Check if it is working
library(CspStandSegmentation)
example("csp_cost_segmentation", run.dontrun = TRUE)

```

## Usage
The package is firmly based on the `lidR` package and uses the las file structure. Smaller point clouds can be directly segmented using the ```csp_cost_segmentation``` function. This requires a set of tree positions (map) as starting points, which can be derived using the ```find_base_coordinates_raster``` function, which might require parameter optimization. Theoretically, tree positions might also come from field measurements or manual assignments.:

```R
# read example data
file = system.file("extdata", "beech.las", package="CspStandSegmentation")
tls = lidR::readTLSLAS(file)

# find tree positions as starting points for segmentation
map <- CspStandSegmentation::find_base_coordinates_raster(tls)

# segment trees
segmented <- tls |>
  CspStandSegmentation::add_geometry(n_cores = parallel::detectCores()/2) |>
  CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = parallel::detectCores()/2)

# show results
lidR::plot(segmented, color = "TreeID")

# create inventory
inventory <- CspStandSegmentation::forest_inventory(segmented)
head(inventory)
lidR::plot(segmented, color = "TreeID") |> CspStandSegmentation::plot_inventory(inventory)
```

For large areas, the package can be used within the lidR LAScatalogue engine to cope with memory limitations. The following example shows how this can be done. The single tiles of segmented trees are saved in a folder in this example and merged afterwards. 

```R
# packages
library(lidR)
library(CspStandSegmentation)

# parameters
las_file <- "your_file.laz"
base_dir <- "~/your_project_folder/" # with trailing /
cores <- parallel::detectCores()/2 # number od cpu cores 
res <- 0.3 # voxel resolution for segmentation
chunk_size <- 50 # size of one tile in m excl. buffer
chunk_buffer <- 10 # buffer around tile in m

# main

# create dir for segmentation tiles
if(!dir.exists(paste0(base_dir,"segmentation_tiles/"))) {
  dir.create(paste0(base_dir,"segmentation_tiles/"))
}

uls = lidR::readTLSLAScatalog(paste0(base_dir,las_file), select = "XYZ0", chunk_size = chunk_size, chunk_buffer = chunk_buffer)
plot(uls, chunk_pattern = TRUE)
# plot(dtm,add = TRUE)
# sf::as_Spatial(sf::st_as_sf(map, coords = c("X", "Y"))) |> plot(add = TRUE)

opt_output_files(uls) <- paste0(base_dir,"segmentation_tiles/{ID}")
segmented <- catalog_apply(uls, function(cluster) {
  
  las <- suppressWarnings(readLAS(cluster)) # read files
  if (is.empty(las) ) return(NULL) # stop if empty
  message(str(cluster))
  # find tree positions as starting points for segmentation
  map <- CspStandSegmentation::find_base_coordinates_raster(las)
  
  # add the tile ID*100,000 to the TreeID to ensure unique IDs across all tiles
  map$TreeID <- map$TreeID + as.numeric(basename(cluster@save)) * 100000
  
  # only use seed positions within the tile+buffer and save the tile bbox to only return tree pos within the tile (excl. buffer)
  inv <- map
  invb <- map
  # the bbox includes the buffer, the bbbox excludes the buffer 
  bbox <- cluster@bbox
  bbbox <- cluster@bbbox
  inv <- inv[inv$X < bbox[1,2] & inv$X > bbox[1,1] & inv$Y < bbox[2,2] & inv$Y > bbox[2,1],]
  invb <- invb[invb$X < bbbox[1,2] & invb$X > bbbox[1,1] & invb$Y < bbbox[2,2] & invb$Y > bbbox[2,1],]
  if (nrow(inv) == 0) return(NULL) # stop if no tree pos in tile found
  if (is.empty(las) ) return(NULL) # stop if empty
  
  # Assign all points to trees
  las <- las |> add_geometry(n_cores = cores) |> csp_cost_segmentation(invb,res, N_cores = cores, V_w = 0.5)
  # las <- las |> csp_cost_segmentation(map,res, N_cores = cores, V_w = 0.5) # this is a faster version which does not make use of the geometric feature weights
  if (is.empty(las)) return(NULL)
  
  las <- las |>  filter_poi(TreeID %in% c(0,inv$TreeID)) # only return trees within the tile
  if (is.empty(las)) return(NULL) # stop if empty
  
  # remove unneccesary attributes for further processing 
  las <- las |> remove_lasattribute('x_vox') |> 
    remove_lasattribute('y_vox') |> 
    remove_lasattribute('z_vox') |> 
    remove_lasattribute('buffer') |>
    remove_lasattribute('Linearity') |>
    remove_lasattribute('Sphericity') |>
    remove_lasattribute('Verticality')
  
  # validate las
  las <- las  |>  las_quantize()  |> las_update()
  if (is.empty(las)) return(NULL)
  return(las)
}, .options = list(automerge = TRUE))

# merge segmented trees
segmented <- readTLSLAScatalog(paste0(base_dir,"segmentation_tiles/"), select = "xyz0", chunk_buffer = 0)
opt_merge(segmented) <- TRUE
opt_output_files(segmented) <- paste0("")
segmented <- catalog_apply(segmented, function(cluster) {
  las <- suppressWarnings(readLAS(cluster)) # read files
  if (is.empty(las) ) return(NULL) # stop if empty
  return(las)
}, .options = list(automerge = TRUE))

# write results in a single file
writeLAS(segmented, paste0(base_dir,"segmented.las"))
```

## Citation
If you publish work related to CspStandSegmentation please cite the following article:

Larysch E, Frey J, Schindler Z, Sprengel L, Hillenmeyer K, Kohnle U, Seifert T, Spiecker H (2025) Quantifying and mapping the ready-to-use veneer volume of European beech trees based on terrestrial laser scanning data. European Journal of Forest Research. https://doi.org/10.1007/s10342-025-01796-z

BibTex:
```
@article{larysch_2025,
	title = {Quantifying and mapping the ready-to-use veneer volume of European beech trees based on terrestrial laser scanning data},
	issn = {1612-4669},
	url = {https://link.springer.com/epdf/10.1007/s10342-025-01796-z},
	doi = {10.1007/s10342-025-01796-z},
	abstract = {Using 3D point clouds obtained with terrestrial laser scanning ({TLS}), we automatically and non-destructively quantified and mapped the estimated veneer wood volume of standing trees in differently structured beech stands. To mitigate climate change, we need to utilise wood for long-term carbon storage in products like construction wood and for substituting building materials based on fossil fuels. As the supply of wood from Norway spruce decreases, alternative species like beech must be considered for construction purposes. We present an approach to quantify and map the volume available for veneer production in beech forests. Our method is based on point clouds derived from {TLS}. We studied three forest plots, each with two different treatments (moderate vs. heavy thinning), resulting in varying stand basal areas ranging from 25 m2 to 36 m2 per hectare. We fitted different configurations of veneer rolls into point clouds of tree stems, choosing the configuration that yielded the highest volume of veneer wood. Our automatic optimisation algorithm ensured no misplaced veneer rolls. At the tree level, veneer wood volume was higher in intensely thinned stands. At the stand level, overall veneer volume was higher in moderately thinned stands, whereas the overall veneer share was higher in the heavily thinned stands. The veneer volume of a tree depended on diameter at breast height, crown base height, taper and curvature depth. Our approach detects all trees in a forest potentially ready for veneer production and shows the direct volumetric outcome under bark. This enables the planning of tree selection for harvest based on adaptable requirements for the veneer production.},
	journaltitle = {European Journal of Forest Research},
	author = {Larysch, Elena and Frey, Julian and Schindler, Zoe and Sprengel, Lars and Hillenmeyer, Katharina and Kohnle, Ulrich and Seifert, Thomas and Spiecker, Heinrich},
	urldate = {2025-06-25},
	date = {2025},
	langid = {english},
	file = {Full Text PDF:O\:\\Research\\Projects\\Confobi_IWW\\Literatur\\lit_database\\storage\\R7Q8BFU5\\Larysch et al. - 2025 - Quantifying and mapping the ready-to-use veneer volume of European beech trees based on terrestrial.pdf:application/pdf},
}
```

