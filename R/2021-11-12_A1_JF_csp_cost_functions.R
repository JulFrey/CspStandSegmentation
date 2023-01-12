# load packages
invisible(lapply(c('lidR','TreeLS', 'dbscan', 'igraph', 'foreach'), require, character.only = TRUE))

# ------------------------------------------------------------------------------

# thx zoe https://github.com/zoeschindler/masterarbeit/blob/main/03_raster_calculation_functions.R
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

voxelize_points_mean_attributes <- function(las, res) {

  # checking resolution input validity
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

add_las_attributes <- function(las) {

  # add attributes from data table permanently to attributes
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
comparative_shortest_path <- function(vox = vox, adjacency_df = adjacency_df, seeds, v_w = 0, l_w = 0, s_w = 0, N_cores = parallel::detectCores() - 1, Voxel_size) {

  # update weights
  adjacency_df$weight <- with(vox@data[adjacency_df$adjacency_list], adjacency_df$weight^2 + ((1 - Verticality) * v_w + Sphericity * s_w + Linearity * l_w) * Voxel_size)
  adjacency_df$weight[adjacency_df$weight < 0] <- 0.01 * Voxel_size # catch negative weights

  #-----------------------
  # compute dijkstra matrix for each seed (trunk)
  # and weigh matrix by DBH^2/3 (Tao et al 2015.)
  #-----------------------

  # build graph
  vox_graph <- adjacency_df |>
    igraph::graph_from_data_frame(directed = F) |>
    igraph::simplify()

  # calculate a distance (weight) graph per seed using dijkstra
  doParallel::registerDoParallel(cores = N_cores)
  dists_list <- foreach::foreach(
    t = 1:nrow(seeds),
    .noexport = c('las', 'map', 'vox', 'tree_seeds', 'ground', 'dtm', 'adjacency_df', 'inv'),
    .errorhandling = c('remove')) %dopar% {
      return(igraph::distances(vox_graph, as.character(seeds$SeedID[t]), algorithm = 'dijkstra'))
    }
  doParallel::stopImplicitCluster()

  # combine to matrix
  dist_matrix <- simplify2array(dists_list)[1,,]

  # get seed with minimum distance
  min_matrix <- apply(dist_matrix, 1, which.min)
  min_dist_matrix <- suppressWarnings(apply(dist_matrix, 1, min, na.rm = T))
  min_matrix <- data.table::data.table(PointID = as.integer(igraph::V(vox_graph)$name), TreeID = seeds$TreeID[as.integer(min_matrix)], dist = min_dist_matrix)
  min_matrix$TreeID[min_dist_matrix == Inf] <- 0 # set SeedIDs 0 for voxels which can't be reached by any seed

  # assign voxels to seeds (minimum cost/distance to trunk)
  vox <- vox |>
    remove_lasattribute('TreeID') |>
    add_attribute(as.integer(rownames(vox@data)), 'PointID')
  vox@data <- merge(vox@data, min_matrix, bz = 'PointID')
  return(vox)
}

# ------------------------------------------------------------------------------

# this is the main function
# it requires a normalized las point cloud of a forest patch with already calculated geometric features
# using the add geometry function,
# a forest inventory as it can be calculated by TreeLS::tlsInventory

csp_cost_segmentation <- function(las, map, Voxel_size = 0.3, V_w = 0, L_w = 0, S_w = 0, N_cores = 1) {

  if ('TreeID' %in% names(las@data)) {
    las <- las |>
      remove_lasattribute('TreeID')
  }

  vox <- voxelize_points_mean_attributes(las, res = Voxel_size)

  if (typeof(map) == 'S4') {
    inv <- map@data[map@data$TreePosition,]
    if (nrow(inv) == 0) {
      inv <- aggregate(map@data[map@data$Z > 1 & map@data$Z < 1.5,], by = list(map@data$TreeID[map@data$Z > 1 & map@data$Z < 1.5]), median)
    }
  } else {
    inv <- map
  }

  # add seeds
  vox <- vox |>
    add_lasattribute(0, 'TreeID', 'TreeID')
  vox@data <- vox@data[,c('X', 'Y', 'Z', 'X_gr', 'Y_gr', 'Z_gr', 'Sphericity', 'Linearity', 'Verticality', 'TreeID')]

  # lift the starting points if LAS is normalized or map doesn't have Z values
  if (sum(inv$Z) == 0 | 'Zref' %in% names(las)) {
    inv$Z <- 1.3
  } else {
    inv$Z <- inv$Z + 1
  }
  inv <- inv |>
    cbind(X_gr = inv$X) |>
    cbind(Y_gr = inv$Y) |>
    cbind(Z_gr = inv$Z) |>
    cbind(Sphericity = 0) |>
    cbind(Linearity = 0) |>
    cbind(Verticality = 0)
  vox@data <- rbind(vox@data, inv[,c('X', 'Y', 'Z', 'X_gr', 'Y_gr', 'Z_gr', 'Sphericity', 'Linearity', 'Verticality', 'TreeID')])

  # possible seeds
  seed_range <- (nrow(vox@data) - nrow(inv) + 1):nrow(vox@data)
  tree_seeds <- data.frame(SeedID = seed_range, TreeID = vox@data$TreeID[seed_range])
  rm(seed_range)

  # use dbscan to calculate a matrix of neighboring points
  neighborhood_list <- dbscan::frNN(vox@data[,4:6], Voxel_size * 2, bucketSize = 22) # voxel size * 1.42 (sqrt(1^2 + 1^2)) 1.73

  # the result has to be disentangled we get the adjacent voxel ids first
  adjacency_list <- unlist(neighborhood_list$id)

  # then we grab the origin voxel using cpp
  adjacency_list_id <- fast_unlist(neighborhood_list$id, length(adjacency_list)) + 1 # +1 because of cpp counting

  # we do the same with the distances
  dists_vec <- fast_unlist_dist(neighborhood_list$dist, length(adjacency_list))

  # compile to a data frame
  adjacency_df <- data.frame(adjacency_list_id,adjacency_list, weight = dists_vec) #, TreeID = vox@data$TreeID[adjacency_list_id]

  # clean
  rm(adjacency_list, adjacency_list_id, dists_vec, neighborhood_list)

  # calculate csp including the weights
  vox2 <- comparative_shortest_path(vox = vox, adjacency_df = adjacency_df, v_w = V_w, l_w = L_w, s_w = S_w, Voxel_size = Voxel_size, N_cores = N_cores, seeds = tree_seeds)

  las <- las |>
    add_voxel_coordinates(Voxel_size) #|> remove_lasattribute('Radius')
  las@data <- merge(las@data, vox2@data[,c('X', 'Y', 'Z', 'TreeID')], by.x = c('x_vox', 'y_vox', 'z_vox'), by.y = c('X', 'Y', 'Z'))
  las <- add_las_attributes(las)
  return(las)
}

# ------------------------------------------------------------------------------

# own function to calculate tree start points to get rid of TreeLS (since its not on cran)
#' Find seeds for tree segmentation based on a density raster and cluster analyses
#'
#' @param las A LAS point cloud object.
#' @param res Resolution of the density raster.
#' @param eps Search radius to combine positive raster cells should be >res.
#' @param q Threshold quantile for a stem in the density raster.
#' @param zmin,zmax Lower and Upper boundary for the density raster.
#' @return A data.frame with seed coordinates.
#' @examples
#' file = system.file("extdata", "pine_plot.laz", package = "TreeLS")
#' tls = readTLS(file) |> classify_ground(csf(), last_returns = F) |> normalize_height(tin())
#' find_seeds(las = tls, res = 0.2, eps = 0.4, q = 0.97, 0.5, 2)
find_base_coordinates_raster <- function(las, res = 0.1, zmin = 0.5, zmax = 2, q = 0.975, eps = 0.2){
  slice <- las |>  filter_poi(Z > zmin & Z < zmax)
  density <- grid_metrics(slice, length(Z), res = res)
  height <- grid_metrics(slice, mean(Z), res = res)
  seed_rast <- terra::as.points(terra::rast(density > quantile(values(density),probs = q, na.rm = T)))
  seed_rast <- terra::subset(seed_rast, seed_rast$layer == 1) |> as.data.frame(geom = 'XY')
  seed_rast <- seed_rast |>  cbind(data.frame(cluster = dbscan::dbscan(seed_rast[,c("x","y")], eps = eps, minPts = 1)$cluster) )
  seed_rast <- aggregate(seed_rast, by = list(seed_rast$cluster), mean)[,3:5]
  z_vals <- extract(height, seed_rast[,1:2])
  seed_rast <- cbind(seed_rast, z_vals)[,c(1,2,4,3)]
  names(seed_rast) <- c('X','Y','Z','TreeID')
  return(seed_rast)
}

# ------------------------------------------------------------------------------

# own function to calculate tree start points to get rid of TreeLS (since its not on cran)
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
