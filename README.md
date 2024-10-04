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
plot(uls, chunk_pattern = T)
# plot(dtm,add = T)
# sf::as_Spatial(sf::st_as_sf(map, coords = c("X", "Y"))) |> plot(add = T)

opt_output_files(uls) <- paste0(base_dir,"segmentation_tiles/{ID}")
segmented <- catalog_apply(uls, function(cluster) {
  
  las <- suppressWarnings(readLAS(cluster)) # read files
  if (is.empty(las) ) return(NULL) # stop if empty
  print(str(cluster))
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
