#packages
invisible(lapply( c("lidR","TreeLS", "dbscan", "igraph", "foreach"), require, character.only = TRUE))


# thx zoe https://github.com/zoeschindler/masterarbeit/blob/main/03_raster_calculation_functions.R
add_geometry <- function(las, k = 10,  n_cores = 1) {
  # necessary for raster_geometry
  # returns geometric features based on eigenvalues
  eigen <- eigen_decomposition(las, k, n_cores) # 20 neighbours, n cores
  las <- add_lasattribute(las, eigen[, 3] / (eigen[, 1] + eigen[, 2] + eigen[, 3]), "Curvature", "curvature")
  las <- add_lasattribute(las, (eigen[, 1] - eigen[, 2]) / eigen[, 1], "Linearity", "linearity")
  las <- add_lasattribute(las, (eigen[, 2] - eigen[, 3]) / eigen[, 1], "Planarity", "planarity")
  las <- add_lasattribute(las, eigen[, 3] / eigen[, 1], "Sphericity", "sphericity")
  las <- add_lasattribute(las, (eigen[, 1] - eigen[, 3]) / eigen[, 1], "Anisotropy", "anisotropy")
  las <- add_lasattribute(las, 1 - abs(eigen[,4]) ,"Verticality","verticality")
  return(las)
}


voxelize_points_mean_attributes = function(las, res){
  Intensity <- NULL

  if (length(res) == 1L)  {
    res <- c(res,res)
  }else if (length(res) > 2L){
    stop("Wrong resolution provided.")
  }


  by <- lidR:::group_grid_3d(las@data$X, las@data$Y, las@data$Z, res, c(0,0,0.5*res[2]))


  voxels <- las@data[, lapply(.SD, mean), by = by]
  if(length(names(las@data)) > 3){
    data.table::setnames(voxels, c("X","Y","Z","X_gr","Y_gr","Z_gr",names(las@data)[4:length(names(las@data))]))
  } else {
    data.table::setnames(voxels, c("X","Y","Z","X_gr","Y_gr","Z_gr"))
  }
  output <- LAS(voxels, header = las@header, crs = st_crs(las), check = FALSE, index = las@index)
  return(output)
}


add_voxel_coordinates <- function(las, res){
  vox <- lidR:::group_grid_3d(las@data$X, las@data$Y, las@data$Z, c(res,res), c(0,0,0.5*res))
  las <- add_lasattribute(las, vox[[1]], "x_vox","x_vox") |>
    add_lasattribute(vox[[2]], "y_vox", "y_vox") |>
    add_lasattribute(vox[[3]], "z_vox", "z_vox")
  return(las)

}

add_las_attributes <- function(las){
  names <- names(las@data)
  names <- names[!(names %in% lidR:::LASATTRIBUTES)]
  for(a in names){
    if(!with(las@data, is.numeric(get(a)))){ next
    }
    las <- las |> add_lasattribute(name = a, desc = a)
  }
  return(las)
}





# V_w, L_W, S_w are the weights for 1-verticality, sphericity, linearity
comparative_shortest_path <- function(vox = vox, adjacency_df = adjacency_df, seeds, v_w = 0,l_w = 0,s_w = 0, N_cores = parallel::detectCores()-1, Voxel_size){

  #update weights
  adjacency_df$weight <- with(vox@data[adjacency_df$adjacency_list], adjacency_df$weight^2 + ((1-Verticality) * v_w + Sphericity * s_w + Linearity * l_w) * Voxel_size)
  adjacency_df$weight[adjacency_df$weight < 0] <- 0.01 * Voxel_size # catch negative weights

  #-----------------------
  # compute dijkstra matrix for each seed (trunk)
  # and weigh matrix by DBH^2/3 (Tao et al 2015.)
  #-----------------------


  # build graph
  vox_graph <- igraph::graph_from_data_frame(adjacency_df, directed = F) |> igraph::simplify()

  # calculate a distance (weight) graph per seed using dijkstra
  doParallel::registerDoParallel(cores = N_cores)
  dists_list <- foreach::foreach(t = 1:nrow(seeds), .noexport = c("las", "map","vox","tree_seeds","ground","dtm",'adjacency_df', 'inv'), .errorhandling=c('remove')) %dopar% {
    return(igraph::distances(vox_graph, as.character(seeds$SeedID[t]), algorithm  = "dijkstra") )
  }
  doParallel::stopImplicitCluster()


  #combine to matrix
  dist_matrix <- simplify2array(dists_list)[1,,]

  #get seed with minimum distance
  min_matrix <- apply(dist_matrix,1, which.min)
  min_dist_matrix <- suppressWarnings(apply(dist_matrix,1, min, na.rm = T))

  min_matrix <- data.table::data.table(PointID = as.integer(igraph::V(vox_graph)$name), TreeID = seeds$TreeID[as.integer(min_matrix)], dist = min_dist_matrix)

  min_matrix$TreeID[min_dist_matrix == Inf] <- 0 # set SeedIDs 0 for voxels which can't be reached by any seed

  # assign voxels to seeds (minimum cost/distance to trunk)


  vox <- vox |>remove_lasattribute("TreeID") |> add_attribute(as.integer(rownames(vox@data)), "PointID")
  vox@data <- merge(vox@data, min_matrix, bz = "PointID")
  return(vox)
}


# this is the main function
# it requires a normalized las point cloud of a forest patch with already calculated geometric features
# using the add geometry function,
# a forest inventory as it can be calculated by TreeLS::tlsInventory

csp_cost_segmentation <- function(las, map, Voxel_size = 0.3, V_w = 0,L_w = 0,S_w = 0, N_cores = 1){

  if("TreeID" %in% names(las@data)) las <- remove_lasattribute(las, "TreeID")

  vox <- voxelize_points_mean_attributes(las, res = Voxel_size)

  if(typeof(map) == "S4"){
    inv <- map@data[map@data$TreePosition,]
    if(nrow(inv) == 0){
      inv <- aggregate(map@data[map@data$Z > 1 & map@data$Z < 1.5,], by = list(map@data$TreeID[map@data$Z > 1 & map@data$Z < 1.5]), median)
    }
  } else {
    inv <- map
  }




  # add Seeds
  vox <- add_lasattribute(vox,0,"TreeID","TreeID")
  vox@data <- vox@data[,c("X", "Y","Z","X_gr", "Y_gr","Z_gr","Sphericity","Linearity","Verticality", "TreeID")]
  if(sum(inv$Z) == 0 | 'Zref' %in% names(las)) {inv$Z <- 1.3} else inv$Z <- inv$Z + 1 #lift the starting points if las is normalized or map doesn't have Z values
  inv <- inv |> cbind(X_gr = inv$X) |> cbind(Y_gr = inv$Y) |> cbind(Z_gr = inv$Z) |> cbind(Sphericity = 0) |> cbind(Linearity = 0) |> cbind(Verticality = 0)
  vox@data <- rbind(vox@data, inv[,c("X", "Y","Z","X_gr", "Y_gr","Z_gr","Sphericity","Linearity","Verticality", "TreeID")])

  #possible seeds
  seed_range <- (nrow(vox@data)-nrow(inv)):nrow(vox@data)
  tree_seeds <- data.frame(SeedID = seed_range, TreeID = vox@data$TreeID[seed_range])
  rm(seed_range)


  # use dbscan to calculate a matrix of neighboring points
  neighborhood_list <- dbscan::frNN(vox@data[,4:6], Voxel_size * 2, bucketSize = 22) # voxel size * 1.42 (sqrt(1^2 + 1^2)) 1.73

  # the result has to be disentangled we get the adjacent voxel ids first
  adjacency_list <- unlist(neighborhood_list$id)
  # then we grab the origin voxel using cpp
  adjacency_list_id <- fast_unlist(neighborhood_list$id, length(adjacency_list))+1 #+1 because of cpp counting
  # we do the same with the distances
  dists_vec <- fast_unlist_dist(neighborhood_list$dist, length(adjacency_list))

  #compile to a data frame
  adjacency_df <- data.frame(adjacency_list_id,adjacency_list, weight = dists_vec) #, TreeID = vox@data$TreeID[adjacency_list_id]

  #clean
  rm(adjacency_list, adjacency_list_id, dists_vec, neighborhood_list)

  # calculate csp including the weights
  vox2 <- comparative_shortest_path(vox = vox, adjacency_df = adjacency_df, v_w = V_w, l_w = L_w, s_w = S_w, Voxel_size = Voxel_size, N_cores = N_cores, seeds = tree_seeds)

  las <- add_voxel_coordinates(las, Voxel_size) #|> remove_lasattribute("Radius")
  las@data <- merge(las@data, vox2@data[,c("X","Y","Z","TreeID")], by.x = c("x_vox","y_vox","z_vox"), by.y = c("X","Y","Z"))
  las <- add_las_attributes(las)
  return(las)
}

# own function to calculate tree start points to get rid of TreeLS (since its not on cran)
find_base_coordinates_raster <- function(las, zmin = 0.5, zmax = 2, res = 0.1, quantile = 0.8, merge_radius = 0.5){
  slice <- las |> filter_poi(Classification != 2 &Z > zmin & Z < zmax)
  if(lidR::is.empty(las)){
    stop('no points found in the specified zmin/xmax range. Is your point cloud normalized?')
  }
  density <- slice |> grid_metrics(~length(Z), res = res)

  density <- density > quantile(density, quantile)
  poly <- density |> terra::rast() |>  terra::as.polygons(dissolve = T) |> sf::st_as_sf()
  centroids <- poly[[2]][2] |> st_cast("POLYGON") |> sf::st_centroid()
  cluster <- sf::st_coordinates(centroids) |> dbscan::dbscan(merge_radius,1)
  centroids <- sf::st_as_sf(centroids)
  centroids['cluster'] <- cluster$cluster
  centroids <- centroids[[1]][which(!duplicated(cluster$cluster))]
  map <- data.frame(st_coordinates(centroids), Z = 0, TreeID = 1:length(centroids))
  return(map)
}


# own function to calculate tree start points to get rid of TreeLS (since its not on cran)
find_base_coordinates_geom <- function(las, zmin = 0.5, zmax = 2, res = 0.5, min_verticality = 0.9, min_planarity = 0.5, min_cluster_size = NULL){

  Zref = T # flag if a normalized point cloud was given
  if(!('Zref' %in% names(las@data))){
    las <- las |> classify_ground(csf(class_threshold = 0.05, cloth_resolution = 0.05), last_returns = F)
    las <- las |> normalize_height(grid_terrain(las, res = 0.25, algorithm = knnidw(), full_raster = TRUE))
    Zref = F
  }


  slice <- las |> filter_poi(Classification != 2 &Z > zmin & Z < zmax)
  if(lidR::is.empty(las)){
    stop('no points found in the specified zmin/xmax range. Is your point cloud normalized?')
  }
  slice <- slice |> add_geometry() |> filter_poi(Planarity > min_planarity & Verticality > min_verticality)

  if(lidR::is.empty(las)){
    stop('no points found in the specified planarity/verticality range. Try lower parameters (>0 & < 1)')
  }

  if(is.null(min_cluster_size)){
    cluster <- slice@data[,1:3] |> dbscan::dbscan(res,1)
    slice@data$Cluster <- cluster$cluster
    slice <- slice |> filter_poi(Cluster %in% unique(cluster$cluster)[table(cluster$cluster)>median(table(cluster$cluster))])
  }
  map <- aggregate(slice@data[,1:2], by = list(slice@data$Cluster), mean)
  if(Zref) {
    Z <- aggregate(slice@data[,3], by = list(slice@data$Cluster), min)} else {
    Z <- aggregate(slice@data[,'Zref'], by = list(slice@data$Cluster), min)
    }
  map <- data.frame(map[, 2:3], Z = Z[,2], TreeID = Z[,1])
}
