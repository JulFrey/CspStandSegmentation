###########
# functions
###########

#' Helper function to compute distances from a point to the circle
#' @param point numeric vector of length 2 c(X,Y)
#' @param circle numeric vector of length 3 c(center_X, center_Y, radius)
#' @return numeric distance from the point to the circle
#' @export point_circle_distance
point_circle_distance <- function(point, circle) {
  return(abs(sqrt(sum((point - circle[1:2])^2)) - circle[3]))
}

#' Returns the angle between the center of the circle and a point in degrees
#' @param point numeric vector of length 2 c(X,Y)
#' @param circle numeric vector of length 3 c(center_X, center_Y, radius)
#' @return numeric angle in degrees
#' @export point_center_angle
point_center_angle <- function(point, circle){
  ang <- atan2(point[2] - circle[2], point[1] - circle[1])
  return(ang * 180 / pi)
}

#' Suppress only the cat() output
#' @param f function to be called
#' @param ... parameters to the function
#' @return the return value of the function
#' @export suppress_cat
suppress_cat <- function(f, ...) {
  null_device <- if (.Platform$OS.type == "windows") "nul" else "/dev/null"
  con <- file(null_device, "w") # Open connection to null device
  sink(con)                     # Redirect output to null device
  on.exit({                     # Ensure cleanup
    sink()
    close(con)
  })
  f(...)                        # Call the function and capture its return value
}

#' RANSAC circle fitting algorithm specially adapted for tree DBH estimation
#'
#' This function fits a circle to a set of points using the RANSAC algorithm it maximizes the points that are in the circle and the number of filled 36 degree angle segments
#' Therefore, this function searches for the most complete circle with the highest number of points represented.
#'
#' @param data numeric matrix with 2 columns (X, Y) representing the point cloud
#' @param n_iterations integer maximum number of iterations
#' @param distance_threshold numeric maximum distance from a point to the circle to be considered an inlier
#' @param min_inliers integer minimum number of inliers to consider the circle as valid
#' @return a list with the following elements:
#' circle: the center coordinates and radius of the circle
#' inliers: number of points within the circles dist threshold
#' angle_segs: number of populated 10deg angular segments of the circle using the distance_threshold
#' n_iter: number of iterations run
#' @export ransac_circle_fit
ransac_circle_fit <- function(data,n_iterations = 1000L,distance_threshold = 0.01,min_inliers = 3L) {
  # Ensure matrix
  data <- as.matrix(data)
  n <- nrow(data)

  # catch if less than 3 points are given
  if (n < 3L) {
    return(list(circle      = c(mean(data[, 1]), mean(data[, 2]), NA_real_),
                inliers     = NA_integer_,
                angle_segs  = NA_integer_,
                n_iter      = 0L))
  }

  # Pre-extract coordinates once
  x <- data[, 1]
  y <- data[, 2]

  # Initialize best circle
  best_circle     <- NULL
  best_inliers    <- 0L
  best_angle_segs <- 0L

  # how many points to sample per iteration
  n_pts_sample <- min(5L, n)

  # point densities for weighted sampling
  p_densities <- dbscan::pointdensity(data, eps = 0.05)

  for (i in seq_len(n_iterations)) {
    # Randomly sample points (weighted)
    idx <- sample.int(n, n_pts_sample, prob = p_densities + 1)
    sample_points <- data[idx, , drop = FALSE]

    # Fit circle; keep tryCatch very tight and avoid pipe
    circle <- tryCatch(
      {
        CspStandSegmentation::suppress_cat(conicfit::CircleFitByPratt, sample_points)
      },
      warning = function(w) NULL,
      error   = function(e) NULL
    )

    # If fit failed, skip iteration
    if (is.null(circle))
      next

    cx <- circle[1L]
    cy <- circle[2L]
    r  <- circle[3L]

    # Vectorised distances: distance from each point to circle
    # dist_to_center = sqrt((x - cx)^2 + (y - cy)^2)
    # distance_to_circle = |dist_to_center - r|
    dx <- x - cx
    dy <- y - cy
    dist_to_center <- sqrt(dx*dx + dy*dy)
    distances <- abs(dist_to_center - r)

    if(any(is.na(c(cx, cy, r, distances)))){ next }

    # Count inliers
    inlier_mask <- distances < distance_threshold
    inliers     <- sum(inlier_mask)

    if (inliers < min_inliers){ next }

    # Vectorised angles only for inliers
    # atan2 in degrees in [0, 360)
    ang <- atan2(dy[inlier_mask], dx[inlier_mask]) * 180 / pi
    ang[ang < 0] <- ang[ang < 0] + 360

    # 10-degree bins
    angle_segs <- length(unique(floor(ang / 10)))

    # Update best circle if it has more angle segments and/or inliers
    if (angle_segs > best_angle_segs ||
        (angle_segs == best_angle_segs && inliers >= best_inliers)) {
      best_circle     <- list(circle     = circle,
                              inliers    = inliers,
                              angle_segs = angle_segs,
                              n_iter     = i)
      best_angle_segs <- angle_segs
      best_inliers    <- inliers
    }

    # optionally early stop: 36 segments is the max for 10° bins
    if (best_angle_segs >= 30L) {
      break
    }
  }

  # if no valid circle found, fall back to mean center, NA radius
  if (is.null(best_circle)) {
    best_circle <- list(circle     = c(mean(x), mean(y), NA_real_),
                        inliers    = 0L,
                        angle_segs = 0L,
                        n_iter     = 0L)
  }

  best_circle
}

#' Function to perform a forest inventory based on a segmented las object (needs to contain TreeID)
#'
#' This function estimates a taper curve for evry tree and returns the DBH at 1.3m, its position in XY coordinates, the tree height and the trees 2D projection area.
#'
#' @param las lidR las object with the segmented trees
#' @param slice_min the minimum height of a slice for stems to estimate the taper curve
#' @param slice_max the maximum height of a slice for stems to estimate the taper curve
#' @param increment the increment of the slices
#' @param width the width of the slices
#' @param max_dbh the maximum DBH allowed
#' @param n_cores number of cores to use
#' @param tree_id_col Column name for the instance segmantation ID (TreeIDs)
#' @param non_tree_id tree_id_col value for non tree elements (can be a vector of IDs)
#' @param use_stem_segmentation logical whether to use only points classified as stem for DBH estimation
#' @param semantic_colname character name of the semantic segmentation column (only needed if use_stem_segmentation is TRUE)
#' @param stem_semantic_label integer semantic label value for stem points (only needed if use_stem_segmentation is TRUE)
#'
#' @returns a data.frame with the TreeID, X, Y, DBH, quality_flag, Height and ConvexHullArea
#'
#' @examples
#' \donttest{
#' # read example data
#' file = system.file("extdata", "beech.las", package="CspStandSegmentation")
#' las = lidR::readTLSLAS(file)
#'
#' # find tree positions as starting points for segmentation
#' map <- CspStandSegmentation::find_base_coordinates_raster(las)
#'
#' # segment trees
#' segmented <- las |>
#'   CspStandSegmentation::add_geometry(n_cores = 2) |>
#'   CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = 2)
#'
#' # perform inventory
#' inventory <- CspStandSegmentation::forest_inventory(segmented, n_cores = 2)
#' }
#' @export forest_inventory
#' @import data.table
forest_inventory <- function(las,
                             slice_min = 0.3, slice_max = 4, increment = 0.2,
                             width = 0.1, max_dbh = 1, n_cores = 1, tree_id_col = "TreeID",
                             non_tree_id = 0,
                             use_stem_segmentation = FALSE,
                             semantic_colname = NULL,
                             stem_semantic_label = NULL) {

  # Temporarily disable data.table progress
  old_opt <- options(datatable.showProgress = FALSE)
  on.exit(options(old_opt), add = TRUE)

  # check valid inputs
  if(width > (increment/2)) {
    stop("Overlapping slices are not allowed. The maximum width is increment*0.5.")
  }
  if (!tree_id_col %in% names(las)) {
    stop("The las object does not contain the tree_id_col attribute")
  }

  # remove no Tree elements
  las@data <- las@data[!is.na(get(tree_id_col)) & !(get(tree_id_col) %in% non_tree_id)]

  # normalize heights if not already done and filter points in the dbh slice
  if ("Zref" %in% names(las)) {
    dbh_slice <- lidR::filter_poi(las, Z > slice_min & Z < slice_max)
    dbh_slice@data$Znorm <- dbh_slice@data$Z
    las@data$Znorm <- las@data$Z
  } else if ("Znorm" %in% names(las)) {
    dbh_slice <- lidR::filter_poi(las, Znorm > slice_min & Znorm < slice_max)
    dbh_slice@data$Zref <- dbh_slice@data$Z
    dbh_slice@data$Z <- dbh_slice@data$Znorm
  } else {
    las_norm <- lidR::normalize_height(lidR::classify_ground(las, lidR::csf()), lidR::tin())
    las <- lidR::add_lasattribute(las, las_norm$Z, "Znorm", "Z normalized")
    dbh_slice <- lidR::filter_poi(las, Znorm > slice_min & Znorm < slice_max)
    dbh_slice@data$Zref <- dbh_slice@data$Z
    dbh_slice@data$Z <- dbh_slice@data$Znorm
  }

  # if specified, filter only stem points
  if(use_stem_segmentation == TRUE){
    if(is.null(semantic_colname) | is.null(stem_semantic_label)){
      stop("If use_stem_segmentation is TRUE, semantic_colname and stem_semantic_label must be provided.")
    }
    if(!semantic_colname %in% names(las)){
      stop("The las object does not contain the specified semantic_colname attribute.")
    }
    dbh_slice <- lidR::filter_poi(dbh_slice, get(semantic_colname) == stem_semantic_label)
  }


  # get BBox for later quality control
  bbox <- c(min_x = las@header$`Min X`, max_x = las@header$`Max X`, min_y = las@header$`Min Y`, max_y = las@header$`Max Y`)

  starts <- seq(slice_min, slice_max - width, by = increment)
  slices <- data.table::data.table(
    slice_id = seq_along(starts),
    z_lo     = starts,
    z_hi     = starts + width
  )

  if(!all(c("Planarity", "Linearity") %in% names(dbh_slice))){
    dbh_slice <- dbh_slice |> CspStandSegmentation::add_geometry(n_cores = n_cores)
  }

  .na2true <- function(x) ifelse(is.na(x) | is.infinite(x),TRUE, x)

  #points_per_stem <- aggregate(dbh_slice$TreeID, by = list(dbh_slice$TreeID),FUN = length)
  points_per_stem <- dbh_slice@data[,.N, by = tree_id_col]

  t1 <- Sys.time()
  IDs <- points_per_stem[N > 3, ..tree_id_col]

  message("Fit a DBH value to every tree.")

  .fit_circle <- function(slice) {
    # 'slice' is a data.table with columns X, Y, Planarity, Verticality

    if (nrow(slice) < 3) {
      return(rep(NA_real_, 3))
    } else if (nrow(slice) < 100) {
      planes <- slice
    } else {
      q <- 1 - sqrt(100 / nrow(slice)) + 0.05
      planes <- slice[
        Planarity  > quantile(Planarity,  q) &
          Verticality > quantile(Verticality, q)
      ]
    }

    ransac_circle_fit(
      planes[, .(X, Y)],
      n_iterations = 500
    )$circle
  }

  .fit_circles <- function(tree) {
    tree[slices,
         on = .(Z >= z_lo, Z < z_hi),{
           pars <- .fit_circle(.SD)
           .(
             r = pars[3L],
             X = pars[1L],
             Y   = pars[2L]
           )
         },
         by = .EACHI,
         .SDcols = c("X", "Y", "Planarity", "Verticality")]
  }

  .spline_predict <- function(tree){

    # return position but no dbh if less than 3 points
    if (nrow(tree) < 3) {
      return(c(X = mean(tree$X), Y = mean(tree$Y), Z = mean(tree$Zref), DBH = NA, quality_flag = 1))
    }

    # fit circles allong trunk
    t_seq <- .fit_circles(tree)
    names(t_seq)[1:2] <- c("Zmin","Zmax")

    # calculate mean dbh height
    Z <- mean(tree$Zref - tree$Znorm) + 1.3

    # quality control: remove unrealistic large circles
    t_seq$r[abs(t_seq$r - median(t_seq$r, na.rm = TRUE)) > 2 * sd(t_seq$r, na.rm = TRUE)] <- NA

    # quality control: check if circles do not get way larger allong the trunk
    for (s in 2:nrow(t_seq)) {
      if (.na2true(t_seq$r[s] > suppressWarnings(min(t_seq$r[1:(s -1)], na.rm = TRUE) * 1.3)) | .na2true(t_seq$r[s] < 0.04) | all(is.na(t_seq$r[1:s]))) {
        t_seq[s, c("X", "Y", "r")] <- NA
      }
    }

    # if there are too little circles to fit a spline return the mean circle diameter
    if ((sum(!is.na(t_seq$r)) <= 3) | (sum(!is.na(t_seq$X)) <= 3) | (sum(!is.na(t_seq$Y)) <= 3)) {
      dbh <- mean(t_seq$r, na.rm = TRUE) * 2
      if (is.na(dbh) | dbh > max_dbh | dbh < 0) {
        return(c(X = mean(tree$X), Y = mean(tree$Y), Z = Z, DBH = NA, quality_flag = 2))
      }
      else {
        return(c(X = mean(tree$X), Y = mean(tree$Y), Z = Z, DBH = dbh, quality_flag = 2))
      }
    }

    # fit a spline for X,Y and DBH each along the Z direction
    spline_r <- suppressWarnings(with(t_seq[!is.na(t_seq$r),
    ], smooth.spline(Zmin, r, df = 3)))
    spline_x <- suppressWarnings(with(t_seq[!is.na(t_seq$X),
    ], smooth.spline(Zmin, X, df = 20)))
    spline_y <- suppressWarnings(with(t_seq[!is.na(t_seq$Y),
    ], smooth.spline(Zmin, Y, df = 20)))
    # predict DBH X and Y at DBH height
    r <- as.numeric(predict(spline_r, 1.3)$y)
    x <- as.numeric(predict(spline_x, 1.3)$y)
    y <- as.numeric(predict(spline_y, 1.3)$y)
    dbh <- r * 2
    if (is.na(dbh) | dbh > max_dbh | dbh < 0) {
      return(c(X = mean(tree$X), Y = mean(tree$Y), Z = Z, DBH = NA, quality_flag = 3))
    } else if (x > (mean(tree$X) + 5) | x < (mean(tree$X) - 5) | y > (mean(tree$Y) + 5) | y < (mean(tree$Y) - 5)) {
      return(c(X = mean(tree$X), Y = mean(tree$Y), Z = Z, quality_flag = 3))
    } else {
      return(c(X = x, Y = y, Z = Z, DBH = dbh, quality_flag = 4))
    }
  }


  dbh_results <- dbh_slice@data[,{
    pars <- .spline_predict(.SD)
    .(
      X = pars[1L],
      Y = pars[2L],
      Z = pars[3L],
      DBH   = pars[4L],
      quality_flag = pars[5L]
    )
  }, by = tree_id_col, .SDcols = c("X", "Y","Z", "Zref","Znorm", "Planarity", "Verticality")]

  message("DBH estimated. Calculating tree heights.")

  .Zdif <- function(x) max(x) - min(x)
  heights <- las@data[,.(Height = .Zdif(Z)), by = tree_id_col]
  message("Tree heights calculated. Calculating convex hull areas.")

  .convex_hull_area <- function(x, y) {
    if (length(x) < 3L) return(NA_real_)
    idx <- chull(x, y)
    xx  <- x[idx]
    yy  <- y[idx]
    0.5 * abs(sum(xx * c(yy[-1], yy[1]) - yy * c(xx[-1], xx[1])))
  }

  dt <- subset(las@data, Znorm > 0.5)
  cpa <- dt[, .(ConvexHullArea = .convex_hull_area(X, Y)), by = tree_id_col]

  dbh_results <- merge(dbh_results, heights, by = tree_id_col)
  dbh_results <- merge(dbh_results, cpa, by = tree_id_col)
  return(dbh_results)
}

#' Function to perform a forest inventory based on a segmented las object (needs to contain TreeID)
#' This version is a faster but more simplistic approach than forest_inventory() for the DBH estimates
#'
#' @param las lidR las object with the segmented trees
#' @param slice_min the minimum height of a DBH slice
#' @param slice_max the maximum height of a DBH slice
#' @param max_dbh the maximum DBH allowed
#' @param n_cores number of cores to use
#'
#' @return a data.frame with the TreeID, X, Y, DBH, quality_flag, Height and ConvexHullArea
#'
#' @export forest_inventory_simple
#' @import data.table
forest_inventory_simple <- function(las, slice_min = 1.2, slice_max = 1.4, max_dbh = 1, n_cores = max(c(1, parallel::detectCores()/2 - 1)))
{
  if (!"TreeID" %in% names(las)) {
    stop("The las object does not contain a TreeID attribute")
  }
  las <- lidR::filter_poi(las, !is.na(TreeID) & TreeID > 0)
  if ("Zref" %in% names(las)) {
    dbh_slice <- lidR::filter_poi(las, Z > slice_min & Z <
                                    slice_max)
  } else if ("Znorm" %in% names(las)) {
    dbh_slice <- lidR::filter_poi(las, Znorm > slice_min &
                                    Znorm < slice_max)
  }  else {
    las_norm <- lidR::normalize_height(lidR::classify_ground(las,
                                                             lidR::csf()), lidR::tin())
    las <- lidR::add_lasattribute(las, las_norm$Z, "Znorm",
                                  "Z normalized")

  }

  na2true <- function(x) ifelse(is.na(x) | is.infinite(x), TRUE, x)

  t1 <- Sys.time()
  message("Fast version.")
  message("Fit a DBH value to every tree:")
  dbh_results <- dbh_slice@data[
    , {
      if (.N < 3L) {
        # too few points for this tree
        .(
          X            = mean(X),
          Y            = mean(Y),
          DBH          = NA_real_,
          quality_flag = 1L
        )
      } else {
        circle <- CspStandSegmentation::ransac_circle_fit(
          cbind(X, Y),
          n_iterations = 500
        )
        .(
          X            = circle$circle[1],
          Y            = circle$circle[2],
          DBH          = min(circle$circle[3] * 2, 2),  # cap at 2
          quality_flag = 4L
        )
      }
    },
    by = TreeID
  ]

  # Add Zdiff column (no need for lidR::add_attribute if you just need it in @data)
  dt <- dbh_slice@data
  data.table::setDT(dt)
  dt[, Zdiff := Z - Znorm]
  dbh_slice@data <- dt

  # Compute mean Zdiff per TreeID
  tree_pos_height <- dbh_slice@data[
    , .(Z = mean(Zdiff, na.rm = TRUE)),  # use na.rm = TRUE if needed
    by = TreeID
  ]

  .Zdif <- function(x) max(x) - min(x)
  heights <- las@data[,.(Height = .Zdif(Z)), by = TreeID]
  message("Tree heights calculated. Calculating convex hull areas.")

  .convex_hull_area <- function(x, y) {
    if (length(x) < 3L) return(NA_real_)
    idx <- chull(x, y)
    xx  <- x[idx]
    yy  <- y[idx]
    0.5 * abs(sum(xx * c(yy[-1], yy[1]) - yy * c(xx[-1], xx[1])))
  }

  dt <- subset(las@data, Znorm > 0.5)
  cpa <- dt[, .(ConvexHullArea = .convex_hull_area(X, Y)), by = TreeID]

  dbh_results2 <- merge(dbh_results, heights, by = "TreeID")
  dbh_results2 <- merge(dbh_results2, cpa, by = "TreeID")
  dbh_results2 <- merge(dbh_results2, tree_pos_height, by = "TreeID")
  return(dbh_results2)
}


#' Function to plot the inventory results into a lidR 3d plot of the point cloud
#' @param plot lidR 3d plot
#' @param inventory data.frame with the inventory results
#' @param cex numeric size of the labels
#' @param label_col character color of the labels
#' @param col color of the spheres
#' @return the plot with the inventory results
#' @examples
#' \donttest{
#' # read example data
#' file = system.file("extdata", "beech.las", package="CspStandSegmentation")
#' las = lidR::readTLSLAS(file)
#'
#' # find tree positions as starting points for segmentation
#' map <- CspStandSegmentation::find_base_coordinates_raster(las)
#'
#' # segment trees
#' segmented <- las |>
#'   CspStandSegmentation::add_geometry(n_cores = 2) |>
#'   CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = 2)
#'
#' # perform inventory
#' inventory <- CspStandSegmentation::forest_inventory(segmented, n_cores = 2)
#'
#' # plot the results
#' \dontrun{
#' x <- lidR::plot(segmented, color = "TreeID")
#' plot_inventory(x, inventory)
#' }
#' }
#'
#' @export plot_inventory
plot_inventory <- function(plot, inventory,col = NA,cex = 1.5, label_col = "white"){
  if(is.na(col)){
    col <- rainbow(max(inventory$TreeID) - min(inventory$TreeID))
  }
  # generate a circle for evry dbh estimation
  rgl::spheres3d(inventory$X - plot[1], inventory$Y - plot[2], inventory$Z+1.3, radius = inventory$DBH/2,col = col)
  for(i in 1:nrow(inventory)){
    rgl::texts3d(inventory$X[i] - plot[1] + inventory$DBH[i], inventory$Y[i] - plot[2] + inventory$DBH[i], inventory$Z[i] + 1.3, text = inventory$TreeID[i], adj = c(0.5,0.5), col = label_col, cex = cex, family = "mono", font = 2)
    rgl::lines3d(c(inventory$X[i] - plot[1], inventory$X[i] - plot[1]), c(inventory$Y[i] - plot[2], inventory$Y[i] - plot[2]), c(inventory$Z[i], inventory$Height[i]), col = ifelse(length(col) >= i, col[i], col), lwd = 2)
  }
}

