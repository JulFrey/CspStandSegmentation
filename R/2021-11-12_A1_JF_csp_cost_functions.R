# load packages
invisible(lapply(c('lidR','dbscan', 'igraph', 'foreach'), require, character.only = TRUE))

# ------------------------------------------------------------------------------

# thx zoe https://github.com/zoeschindler/masterarbeit/blob/main/03_raster_calculation_functions.R

#' Add geometric features to a LAS object
#'
#' The function calls a fast cpp multi-core function to calculate eigenvalues
#' for the points in a point cloud based on the k nearest neighbors. Afterwards
#' it adds geometric features like Curvature, Linearity, Planarity, Sphericity,
#' Anisotrophy and Verticlity to the points itself.
#'
#' Details of the metrics can be found in: \ Hackel, T., Wegner, J.D. &
#' Schindler, K. (2016) Contour Detection in Unstructured 3D Point Clouds. In
#' 2016 IEEE Conference on Computer Vision and Pattern Recognition (CVPR).
#' Presented at the 2016 IEEE Conference on Computer Vision and Pattern
#' Recognition (CVPR), IEEE, Las Vegas, NV, USA, pp. 1610–1618.
#'
#' @param las A LAS object (see lidR::LAS)
#' @param k the k nearest neighbors to use for the eigenvalue calculation
#' @param n_cores The number of CPU cores to use
#' @return The function returns a single LAS object with the geometric features
#' attached to it in the LAS@data section.
#' @note %% ~~further notes~~
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @examples
#'
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- lidR::readLAS(LASfile, select = "xyz", filter = "-inside 481250 3812980 481300 3813030")
#'
#' las <- add_geometry(las, k = 5, n_cores = parallel::detectCores()-1)
#' summary(las@data)
#'
#'
#' @export add_geometry
add_geometry <- function(las, k = 10, n_cores = 1) {
  # necessary for raster_geometry
  # returns geometric features based on eigenvalues
  eigen <- eigen_decomposition(las, k, n_cores) # k neighbours, n cores
  las <- las |>
    add_lasattribute(eigen[,3] / (eigen[,1] + eigen[,2] + eigen[, 3]), 'Curvature', 'curvature') |>
    add_lasattribute((eigen[,1] - eigen[,2]) / eigen[,1], 'Linearity', 'linearity') |>
    add_lasattribute((eigen[,2] - eigen[,3]) / eigen[,1], 'Planarity', 'planarity') |>
    add_lasattribute(eigen[,3] / eigen[,1], 'Sphericity', 'sphericity') |>
    add_lasattribute((eigen[,1] - eigen[,3]) / eigen[,1], 'Anisotropy', 'anisotropy') |>
    add_lasattribute(1 - abs(eigen[,4]) ,'Verticality','verticality')
  return(las)
}

# ------------------------------------------------------------------------------

#' helper function to voxelize a las element
#'
#' Calculate voxel mean values for all numeric attributes in the las@data table
#' including the XYZ-coordinates.
#'
#' Returns a las element with XYZ-coordinates as the voxel center and
#' X_gr,Y_gr,Z_gr as the center of gravity (mean point coordinates) as well as
#' all other numeric columns voxel mean values with their original name.
#'
#' @param las a lidR::LAS element
#' @param res voxel resolution in meter
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @seealso \code{\link{voxelize_points}}
#' @examples
#'
#' # read example data
#' file = system.file("extdata", "beech.las", package="CspStandSegmentation")
#' tls = lidR::readTLSLAS(file)
#' tls |> voxelize_points_mean_attributes(1) |> lidR::plot(color = 'X_gr')
#'
#' @export voxelize_points_mean_attributes
voxelize_points_mean_attributes <- function(las, res) {

  # Checking resolution input validity
  if (length(res) == 1L) {
    res <- c(res, res)
  } else if (length(res) > 2L) {
    stop('Wrong resolution provided.')
  }

  # create voxel coordinates
  by <- lidR:::group_grid_3d(las@data$X, las@data$Y, las@data$Z, res, c(0, 0, 0.5*res[2]))

  # add mean attributes
  voxels <- las@data[,lapply(.SD, mean), by = by]
  if (length(names(las@data)) > 3) {
    data.table::setnames(voxels, c('X', 'Y', 'Z', 'X_gr', 'Y_gr', 'Z_gr', names(las@data)[4:length(names(las@data))]))
  } else {
    data.table::setnames(voxels, c('X', 'Y', 'Z', 'X_gr', 'Y_gr', 'Z_gr'))
  }

  # convert voxels to LAS object
  output <- LAS(voxels, header = las@header, crs = st_crs(las), check = FALSE, index = las@index)
  return(output)
}

# ------------------------------------------------------------------------------

#' Add voxel coordinates to a las file
#'
#' Adds the collums x_vox, y_vox and z_vox in the given ressolution to the las
#' element. This is convenient if information has been derived in voxel space
#' and these should be attached to the original points.
#'
#' Voxel coordinates derived with this function are identical to those derived
#' by lidR::voxelize.
#'
#' @param las an element of lidR::LAS class
#' @param res voxel ressolution in [m]
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @examples
#'
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- lidR::readLAS(LASfile, select = "xyz", filter = "-inside 481250 3812980 481300 3813030")
#'
#' las <- add_voxel_coordinates(las,res = 1)
#'
#' lidR::plot(las, color = 'z_vox')
#'
#' @export add_voxel_coordinates
add_voxel_coordinates <- function(las, res) {

  # create voxel coordinates
  vox <- lidR:::group_grid_3d(las@data$X, las@data$Y, las@data$Z, c(res, res), c(0, 0, 0.5*res))

  # add voxel coordinates to LAS
  las <- las |>
    add_lasattribute(vox[[1]], 'x_vox', 'x_vox') |>
    add_lasattribute(vox[[2]], 'y_vox', 'y_vox') |>
    add_lasattribute(vox[[3]], 'z_vox', 'z_vox')
  return(las)
}

# ------------------------------------------------------------------------------

#' Add all las_attributes from las@data to the header of a las element
#'
#' The helper function adds all headings from las@data which are not part of
#' lidR:::LASATTRIBUTES to the las header using lidR::add_lasattribute. Only
#' attributes that are included in the header got saved when using
#' lidR::writeLAS, this is a convenient way to add them.
#'
#' @param las an element of lidR::LAS class
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @examples
#'
#' LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
#' las <- lidR::readLAS(LASfile, select = "xyz", filter = "-inside 481250 3812980 481300 3813030")
#'
#' las@data$noise <- runif(nrow(las@data))
#' las@data$noiseZ <- las@data$var1 * las@data$Z
#'
#' las <- add_las_attributes(las)
#'
#' @export add_las_attributes
add_las_attributes <- function(las) {

  # Add attributes from data.table permanently to attributes
  names <- names(las@data)
  names <- names[!(names %in% lidR:::LASATTRIBUTES)]
  for (name in names) {
    if (!with(las@data, is.numeric(get(name)))) {
      next
    }
    las <- las |>
      add_lasattribute(name = name, desc = name)
  }
  return(las)
}

# ------------------------------------------------------------------------------

# V_w, L_W, S_w are the weights for 1-verticality, sphericity, linearity

#' helper function for csp_cost_segemntation
#'
#' The function performs a Dijkstra algorithm on a 3D voxel file to assign
#' every voxel to the closest seed point using the igraph package.
#'
#' @param vox a LAS S4 element with XYZ voxel coordinates in the @data slot.
#' @param adjacency_df a data.frame with voxel ids (row numbers) in the first
#' column and a neighboring voxel ID in the second column and the weight
#' (distance) in the third column. Might be generated using the dbscan::frNN
#' function (which requires reshaping the data).
#' @param seeds seed points for tree positions.
#' @param v_w,l_w,s_w weights for verticality, linearity spericity see
#' \code{\link{csp_cost_segmentation}}
#' @param N_cores Number of CPU cores for multi-threading
#' @param Voxel_size Edge length used to create the voxels. This is only
#' important to gain comparable distance weights on different voxel sizes.
#' Should be greater than 0.
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @seealso \code{\link{csp_cost_segmentation}}
#' @examples
#'
#' @export comparative_shortest_path
comparative_shortest_path <- function(vox = vox, adjacency_df = adjacency_df, seeds, v_w = 0, l_w = 0, s_w = 0, N_cores = parallel::detectCores() - 1, Voxel_size) {

  # update weights
  adjacency_df$weight <- with(vox@data[adjacency_df$adjacency_list], adjacency_df$weight^2 + ((1 - Verticality) * v_w + Sphericity * s_w + Linearity * l_w) * Voxel_size)
  adjacency_df$weight[adjacency_df$weight < 0] <- 0.001 * Voxel_size # catch negative weights

  # set distances to seeds 0
  adjacency_df$weight[adjacency_df$adjacency_list_id %in% seeds$SeedID] <- 0

  #-----------------------
  # compute dijkstra matrix for each seed (trunk)
  # and weigh matrix by DBH^2/3 (Tao et al 2015.)
  #-----------------------

  # build graph
  vox_graph <- adjacency_df |>
    igraph::graph_from_data_frame(directed = F) |>
    igraph::simplify()

  # calculate a distance (weight) graph per seed using Dijkstra
  doParallel::registerDoParallel(cores = N_cores)
  dists_list <- foreach::foreach(
    t = 1:nrow(seeds),
    .noexport = c('las', 'map', 'vox', 'tree_seeds', 'ground', 'dtm', 'adjacency_df', 'inv'),
    .errorhandling = c('pass')) %dopar% {
      return(igraph::distances(vox_graph, as.character(seeds$SeedID[t]), algorithm = 'dijkstra'))
    }
  doParallel::stopImplicitCluster()

  unreachable <- which(sapply(dists_list,function(x) is.character(x[[1]])))
  if(length(unreachable) > 0) {
    warning('Not all base positions could be reached by the graph. Try a lower resolution or a different approach to find tree base positions. Error messages for the unreachable TreeIDs:')
    warning(paste0(unreachable , paste(":",paste(dists_list[unreachable]), collapse = " "), "\n"), call. = F)
    seeds <- seeds[-unreachable,]
    dists_list <- dists_list[-unreachable]
    }

  # Combine to matrix
  dist_matrix <- simplify2array(dists_list)[1,,]

  # get seed with minimum distance
  min_matrix <- apply(dist_matrix, 1, which.min)
  min_dist_matrix <- suppressWarnings(apply(dist_matrix, 1, min, na.rm = T))
  min_matrix <- data.table::data.table(PointID = as.integer(igraph::V(vox_graph)$name), TreeID = seeds$TreeID[as.integer(min_matrix)], dist = min_dist_matrix)
  min_matrix$TreeID[min_dist_matrix == Inf] <- 0 # set SeedIDs 0 for voxels which any seed can't reach

  # assign voxels to seeds (minimum cost/distance to trunk)
  vox <- vox |>
    remove_lasattribute('TreeID') |>
    add_attribute(as.integer(rownames(vox@data)), 'PointID')
  vox@data <- merge(vox@data, min_matrix, by = 'PointID')
  return(vox)
}

# ------------------------------------------------------------------------------

# This is the main function
# It requires a las point cloud of a forest patch and
# a forest inventory as it can be calculated by CspStandSegmentation::find_base_coordinates_raster
# Preferable geometric features should be computed for the point cloud prior
# to the use of this function using 'add_geometry()'.

#' Comparative Shortest Path with cost weighting tree segmentation
#'
#' Segments single trees from forest point clouds based on tree positions
#' (xyz-coordinates) provided in the map argument.
#'
#' The whole point cloud is voxelized in the given resolution and the center of
#' gravity for the points inside is calculated as voxel coordinate. A graph is
#' build, which connects the voxel coordinates based on a db-scan algorithm. The
#' distances between the voxel coordinates is weighted based on geometric
#' features computed for the points in the voxels. Distances along planar
#' and/or vertical faces like stems are weighted shorter than distances through
#' voxels with high sphericity like leaves and clusters of twigs. This
#' avoids small individuals spreading into the upper canopy.
#' For every voxel center, the weighted distance in the network is calculated to
#' all tree locations from the map argument. The TreeID of the map argument
#' with the shortest distance is assigned to the voxel. All points in the point
#' cloud receive the TreeID from their parent voxel.
#'
#' @param las A lidR LAS S4 object.
#' @param map A data.frame, including the columns
#' X, Y, Z, TreeID, with X and Y depicting the location of the trees. Can be generated using CspStandSegmentation::find_base_coordinates_raster
#' @param Voxel_size The voxel size (3D resolution) for the routing graph to
#' determine the nearest map location for every point in the point cloud.
#' @param V_w verticality weight. Since trunks are vertical structures, routing
#' through voxels with high verticality can be rated 'cheaper.' should be a
#' number between 0 and 1 with 0 meaning no benefit for more vertical
#' structures.
#' @param L_w Linearity weight. Similar to V_w but for linearity, higher
#' values indicate a malus for linear shapes (usually branches).
#' @param S_w Spericity weight. Similar to V_w but for sphericity, higher
#' values indicate a malus for spherical shapes (usually small branches and
#' leaves).
#' @param N_cores number of CPU cores used for parallel routing using the
#' foreach package.
#' @return Returns a copy of the las point cloud with an additional field for
#' the TreeID.
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @seealso \code{\link{comparative_shortest_path}}
#'
#' @value The input LAS object with the TreeID field added to the @data slot.
#'
#' @examples
#'
#' # read example data
#' file = system.file("extdata", "beech.las", package="CspStandSegmentation")
#' tls = lidR::readTLSLAS(file)
#'
#' # Find tree positions as starting points for segmentation
#' map <- CspStandSegmentation::find_base_coordinates_raster(tls)
#' # segment trees
#' segmented <- tls |>
#' CspStandSegmentation::add_geometry() |>
#'   CspStandSegmentation::csp_cost_segmentation(map, 1)
#'
#' lidR::plot(segmented, color = "TreeID")
#'
#' @export csp_cost_segmentation
csp_cost_segmentation <- function(las, map, Voxel_size = 0.3, V_w = 0, L_w = 0, S_w = 0, N_cores = 1) {
  # check if TreeID already exists
  if ('TreeID' %in% names(las@data)) {
    warning("'las' already contains TreeIDs', which will be overwritten.")
    las <- las |>
      lidR::remove_lasattribute('TreeID')
  }
  # Check if geometric features exist in las and compute dummies if not
  if(!all(c('Sphericity', 'Linearity', 'Verticality') %in% names(las@data))) {
    warning("no geometric features in las.  V_w, L_w and/or S_w weights will be ignored.Use 'las <- las |> add_gemetry()' prior to calling this function.")
    las <- las |>
      lidR::add_lasattribute(0, 'Sphericity', 'Sphericity') |>
      lidR::add_lasattribute(0, 'Linearity', 'Linearity') |>
      lidR::add_lasattribute(0, 'Verticality', 'Verticality')
  }

  # Voxelize las with mean attributes
  vox <- voxelize_points_mean_attributes(las, res = Voxel_size)

  if (typeof(map) == 'S4') {
    inv <- map@data[map@data$TreePosition,]
    if (nrow(inv) == 0) {
      inv <- aggregate(map@data[map@data$Z > 1 & map@data$Z < 1.5,], by = list(map@data$TreeID[map@data$Z > 1 & map@data$Z < 1.5]), median)
    }
  } else {
    inv <- map
  }

  # Add seeds
  vox <- vox |>
    lidR::add_lasattribute(0, 'TreeID', 'TreeID')
  vox@data <- vox@data[,c('X', 'Y', 'Z', 'X_gr', 'Y_gr', 'Z_gr', 'Sphericity', 'Linearity', 'Verticality', 'TreeID')]

  # Lift the starting points if map doesn't have Z values
  if (sum(inv$Z) == 0) {
    inv$Z <- 0.5
  }
  if(las@header$`Min Z` > max(inv$Z)) {stop('The minimum Z value of the point cloud is higher than the tree positions. Is the inventory without Z values and the point clound not normalized?')}

  inv <- inv |>
    cbind(X_gr = inv$X) |>
    cbind(Y_gr = inv$Y) |>
    cbind(Z_gr = inv$Z) |>
    cbind(Sphericity = 0) |>
    cbind(Linearity = 0) |>
    cbind(Verticality = 0)
  vox@data <- rbind(vox@data, inv[,c('X', 'Y', 'Z', 'X_gr', 'Y_gr', 'Z_gr', 'Sphericity', 'Linearity', 'Verticality', 'TreeID')])

  # Seed positions
  seed_range <- (nrow(vox@data) - nrow(inv) + 1):nrow(vox@data)
  tree_seeds <- data.frame(SeedID = seed_range, TreeID = vox@data$TreeID[seed_range])
  rm(seed_range)

  # Use dbscan to calculate a matrix of neighboring points
  neighborhood_list <- dbscan::frNN(vox@data[,c('X_gr', 'Y_gr', 'Z_gr')], Voxel_size * 2, bucketSize = 22)
  # dbscan::frNN(vox@data[,c('X_gr', 'Y_gr', 'Z_gr')], Voxel_size * 2, bucketSize = 22)  # voxel size * 1.42 (sqrt(1^2 + 1^2)) 1.73

  # The result has to be disentangled we get the adjacent voxel IDs first
  adjacency_list <- unlist(neighborhood_list$id)

  # Then we grab the origin voxel using cpp
  adjacency_list_id <- fast_unlist(neighborhood_list$id, length(adjacency_list)) + 1 # +1 because of cpp counting

  # We do the same with the distances
  dists_vec <- fast_unlist_dist(neighborhood_list$dist, length(adjacency_list))

  # Compile to a data frame
  adjacency_df <- data.frame(adjacency_list_id,adjacency_list, weight = dists_vec) #, TreeID = vox@data$TreeID[adjacency_list_id]
  rm(adjacency_list, adjacency_list_id, dists_vec, neighborhood_list)

  # Calculate CSP, including the weights
  vox2 <- comparative_shortest_path(vox = vox, adjacency_df = adjacency_df, v_w = V_w, l_w = L_w, s_w = S_w, Voxel_size = Voxel_size, N_cores = N_cores, seeds = tree_seeds)

  las <- las |>
    add_voxel_coordinates(Voxel_size)
  las@data <- merge(las@data, vox2@data[,c('X', 'Y', 'Z', 'TreeID')], by.x = c('x_vox', 'y_vox', 'z_vox'), by.y = c('X', 'Y', 'Z'))
  las <- add_las_attributes(las)
  return(las)
}

# ------------------------------------------------------------------------------

# Function to calculate tree start points based on a raster density approach
#'
#' Find stem base position using a density raster approach
#' @param las an element of lidR::LAS class
#' @param zmin lower search boundary
#' @param zmax upper search boundary
#' @param q quantile of raster density to assign a tree region
#' @param eps search radius to merge base points
#' @param res raster resolution
#' @param quantile raster density quantile to assign a tree region
#' @param merge_radius search radius to merge base points
#' @param normalized boolean flag if the point cloud is normalized
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @examples
#' @export find_base_coordinates_raster
find_base_coordinates_raster <- function(las, res = 0.1, zmin = 0.5, zmax = 2, q = 0.975, eps = 0.2, normalized = F){
  if (!normalized) {
    las <- las |>
      lidR::classify_ground(lidR::csf(class_threshold = 0.05, cloth_resolution = 0.05), last_returns = F)
    dtm <- lidR::rasterize_terrain(las, 0.5, lidR::tin())
    las <- las |> lidR::normalize_height(lidR::tin(), dtm = dtm)
  }
  slice <- las |>  lidR::filter_poi(Z > zmin & Z < zmax)
  density <- lidR::grid_metrics(slice, length(Z), res = res)
  height <- lidR::grid_metrics(slice, mean(Z), res = res)
  seed_rast <- terra::as.points(terra::rast(density > quantile(terra::values(density),probs = q, na.rm = T)))
  seed_rast <- terra::subset(seed_rast, seed_rast$layer == 1) |> as.data.frame(geom = 'XY')
  seed_rast <- seed_rast |>  cbind(data.frame(cluster = dbscan::dbscan(seed_rast[,c("x","y")], eps = eps, minPts = 1)$cluster) )
  seed_rast <- aggregate(seed_rast, by = list(seed_rast$cluster), mean)[,3:5]
  if(normalized){
    z_vals <- terra::extract(height, seed_rast[,1:2])
    if(any(is.na(z_vals))) z_vals[is.na(z_vals)] <- mean(c(zmin, zmax)) # catch NA's
  } else {
    z_vals <- terra::extract(dtm, seed_rast[,1:2])[,2]
    if(any(is.na(z_vals))) z_vals[is.na(z_vals)] <- mean(terra::values(dtm)) # catch NA's
  }

  seed_rast <- cbind(seed_rast, z_vals)[,c(1,2,4,3)]
  colnames(seed_rast) <- c('X','Y','Z','TreeID')
  return(seed_rast)
}

# ------------------------------------------------------------------------------

# own function to calculate tree start points

#' Find stem base position using a geometric feature filtering and clustering
#' approach
#' @param las an element of lidR::LAS class
#' @param zmin lower search boundary
#' @param zmax upper search boundary
#' @param res cluster search radius
#' @param min_verticality minimum verticality >0 & <1 for a point to be
#' considered a stem point
#' @param min_planarity minimum planarity >0 & <1 for a point to be considered
#' a stem point
#' @param min_cluster_size minimum number of points in the cluster to be considered
#' a tree, if NULL median cluster size is chosen
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @author Julian Frey <julian.frey@@wwd.uni-freiburg.de>
#' @examples
#' @export find_base_coordinates_geom
find_base_coordinates_geom <- function(las, zmin = 0.5, zmax = 2, res = 0.5, min_verticality = 0.9, min_planarity = 0.5, min_cluster_size = NULL) {

  Zref <- T # flag if a normalized point cloud was given
  if (!('Zref' %in% names(las@data))) {
    las <- las |>
      classify_ground(csf(class_threshold = 0.05, cloth_resolution = 0.05), last_returns = F)
    las <- las |>
      normalize_height(grid_terrain(las, res = 0.25, algorithm = knnidw(), full_raster = TRUE))
    Zref <- F
  }

  slice <- las |>
    filter_poi(Classification != 2 & Z > zmin & Z < zmax)
  if (lidR::is.empty(slice)) {
    stop('No points found in the specified zmin/xmax range. Is your point cloud normalized?')
  }

  slice <- slice |>
    add_geometry() |>
    filter_poi(Planarity > min_planarity & Verticality > min_verticality)
  if (lidR::is.empty(slice)) {
    stop('No points found in the specified planarity/verticality range. Try lower parameters (> 0 & < 1)')
  }

  if (is.null(min_cluster_size)) {
    cluster <- slice@data[,1:3] |>
      dbscan::dbscan(res, 1)
    slice@data$Cluster <- cluster$cluster
    slice <- slice |>
      filter_poi(Cluster %in% unique(cluster$cluster)[table(cluster$cluster) > median(table(cluster$cluster))])
  }
  map <- aggregate(slice@data[,1:2], by = list(slice@data$Cluster), mean)
  if (Zref) {
    Z <- aggregate(slice@data[,3], by = list(slice@data$Cluster), min)
  } else {
    Z <- aggregate(slice@data[,'Zref'], by = list(slice@data$Cluster), min)
  }
  map <- data.frame(map[,2:3], Z = Z[,2], TreeID = Z[,1])
}
