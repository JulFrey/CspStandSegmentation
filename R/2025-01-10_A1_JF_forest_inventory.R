###########
# functions
###########

#' Helper function to compute distances from a point to the circle
#' @param point numeric vector of length 2 c(X,Y)
#' @param circle numeric vector of length 3 c(center_X, center_Y, radius)
#' @return numeric distance from the point to the circle
point_circle_distance <- function(point, circle) {
  return(abs(sqrt(sum((point - circle[1:2])^2)) - circle[3]))
}

#' Returns the angle between the center of the circle and a point in degrees
#' @param point numeric vector of length 2 c(X,Y)
#' @param circle numeric vector of length 3 c(center_X, center_Y, radius)
#' @return numeric angle in degrees
point_center_angle <- function(point, circle){
  ang <- atan2(point[2] - circle[2], point[1] - circle[1])
  return(ang * 180 / pi)
}

#' Suppress only the cat() output
#' @param f function to be called
#' @return the return value of the function
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
ransac_circle_fit <- function(data, n_iterations = 1000, distance_threshold = 0.01, min_inliers = 3) {
  # catch if less than 3 points are given
  if (nrow(data) < 3) {
    return(list(circle = c(mean(data[,1]), mean(data[,2]), NA), inliers = NA, angle_segs = NA, n_iter = 0))
  }

  # Initialize best circle
  best_circle <- NULL
  best_inliers <- 0
  best_angle_segs <- 0
  n_pts_sample <- min(c(5, nrow(data)))

  # calculate point densities in neighbourhood
  p_densities <- dbscan::pointdensity(data, eps = 0.05)


  for (i in 1:n_iterations) {
    # Randomly sample 3 points
    sample_points <- data[sample(1:nrow(data), n_pts_sample, prob = p_densities + 1), ] |> as.matrix()

    # Compute circle parameters
    tryCatch({
      circle <- suppress_cat(conicfit::CircleFitByPratt, sample_points)

      # Count inliers
      distances <- apply(data, 1, point_circle_distance, circle = circle)
      inliers <- sum(distances < distance_threshold)
      angle_segs <- length(unique(floor(apply(data, 1, point_center_angle, circle = circle)/10)))

      # # Update best circle if more inliers are found
      # if (inliers > best_inliers && inliers >= min_inliers) {
      #   best_circle <- list(circle = circle, inliers = inliers)
      #   best_inliers <- inliers
      # }

      # Update best circle if it has more angle segments and/or inliers
      if (angle_segs >= best_angle_segs && inliers >= best_inliers) {
        best_circle <- list(circle = circle, inliers = inliers, angle_segs = angle_segs, n_iter = i)
        best_angle_segs <- angle_segs
        best_inliers <- inliers
      }

      if((best_angle_segs == nrow(data) | best_angle_segs >= 30) ){ #& inliers == nrow(data)
        break
      }

    },
    warnings = function(w){},
    error = function(e) {}) |> suppressWarnings()
  }
  best_circle$n_iter <- i
  return(best_circle)
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
#'
#' @returns a data.frame with the TreeID, X, Y, DBH, quality_flag, Height and ConvexHullArea
#'
#' @examples
#' # read example data
#' file = system.file("extdata", "beech.las", package="CspStandSegmentation")
#' tls = lidR::readTLSLAS(file)
#'
#' # find tree positions as starting points for segmentation
#' map <- CspStandSegmentation::find_base_coordinates_raster(tls)
#'
#' # segment trees
#' segmented <- tls |>
#'   CspStandSegmentation::add_geometry(n_cores = parallel::detectCores()/2) |>
#'   CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = parallel::detectCores()/2)
#'
#' # show results
#' \dontrun{lidR::plot(segmented, color = "TreeID")}
#'
#' # perform inventory
#' inventory <- CspStandSegmentation::forest_inventory(segmented)
forest_inventory <- function(tls, slice_min = 0.3, slice_max = 4, increment = 0.2, width = 0.1, max_dbh = 1, n_cores = max(c(1,parallel::detectCores()/2 - 1))){
  # check if the TreeID attribute is present

  if(!"TreeID" %in% names(tls)){
    stop("The las object does not contain a TreeID attribute")
  }
  tls <- lidR::filter_poi(tls, !is.na(TreeID) & TreeID > 0)
  if("Zref" %in% names(tls)){
    dbh_slice <- tls |> lidR::filter_poi(Z > slice_min & Z < slice_max)
  } else if("Znorm" %in% names(tls)){
    dbh_slice <- tls |> lidR::filter_poi(Znorm > slice_min & Znorm < slice_max)
  } else {
    tls_norm <- tls |> lidR::classify_ground(lidR::csf()) |> lidR::normalize_height(lidR::tin())
    tls <- tls |> lidR::add_lasattribute(tls_norm$Z, "Znorm", "Z normalized")
    dbh_slice <- tls |> lidR::filter_poi(Znorm > slice_min & Znorm < slice_max)
  }

  # create sequence of heights
  seq <- seq(slice_min, slice_max, by = increment)
  seq <- data.frame(id = 1:length(seq), Zmin = seq- width *0.5, Zmax = seq + width*0.5)

  # DBH estimation
  #dbh_slice <- dbh_slice |> CspStandSegmentation::add_geometry(n_cores = ceiling(parallel::detectCores()/2 - 1))
  na2true <- function(x) ifelse(is.na(x) | is.infinite(x), TRUE, x)

  points_per_stem <- aggregate(dbh_slice$TreeID, by = list(dbh_slice$TreeID), FUN = length)

  t1 <- Sys.time()
  #dbh_results <- data.frame(TreeID = numeric(), X = numeric(), Y = numeric(), DBH = numeric(), quality_flag = numeric())

  IDs <- points_per_stem$Group.1[points_per_stem$x > 3]

  las_list <- lapply(IDs, function(id) tls |> lidR::filter_poi(TreeID == id))
  # initialize progress bar
  print("Fit a DBH value to every tree:")
  #pb = txtProgressBar(min = 0, max = length(IDs), initial = 0, style = 3)
  #i <- 1

  require(foreach)
  cl=parallel::makeCluster(n_cores) #I have 8 core 16 threads but use 14 here. Please edit this accordingly.
  doParallel::registerDoParallel(cl)
  dbh_results <- foreach::foreach(tree = las_list, .combine=rbind, .export = c("ransac_circle_fit", "seq", "max_dbh")) %do% { #, .errorhandling = "remove"
    # get the tree points
    t <- tree@data$TreeID[1]
    if(nrow(tree) < 3){
      dbh <- NA
      return(data.frame(TreeID = t,X = mean(tree$X), Y = mean(tree$Y), DBH = dbh, quality_flag = 1))
      warning(paste("to little points in tree",t))
    }

    # create sequence of heights
    seq <- seq(slice_min, slice_max, by = increment)
    seq <- data.frame(id = 1:length(seq), Zmin = seq- width *0.5, Zmax = seq + width*0.5)
    t_seq <- cbind(seq, X = NA,Y = NA, r = NA)

    rs <- numeric()
    xs <- numeric()
    ys <- numeric()
    for(s in seq$id){
      # get the slice
      slice <- tree |> lidR::filter_poi(Znorm > t_seq$Zmin[s] & Znorm < t_seq$Zmax[s])

      if(nrow(slice) < 3){
        rs <- c(rs, NA)
        xs <- c(xs, NA)
        ys <- c(ys, NA)
        next
      } else if(nrow(slice) < 100){
        planes <- slice
      } else {
        # estimate the planarity and verticality of the points
        slice <- slice |> CspStandSegmentation::add_geometry()
        # use the most planar and vertical points to estimate the DBH
        planes <- slice |> lidR::filter_poi(Planarity > quantile(Planarity, 0.95) & Verticality > quantile(Verticality, 0.95))
        q <- 0.95
        while(nrow(planes) < 100){
          q <- q - 0.05
          planes <- slice |> lidR::filter_poi(Planarity > quantile(Planarity, q) & Verticality > quantile(Verticality, q))
          if(q < 0.2){
            planes <- slice
            break
          }
        }
      }
      planes <- slice
      circle <- ransac_circle_fit(planes@data[, c("X", "Y")], n_iterations = 500)
      # circle
      # plot(Y ~ X, data = tree@data, col = "black", pch = ".", asp = 1)
      # points(planes@data$X, planes@data$Y, col = "cornflowerblue", pch = 20)
      # if(!is.null(circle)){
      #   circ <- conicfit::calculateCircle(circle$circle[1],circle$circle[2],circle$circle[3])
      #   points(circ[,1], circ[,2], pch = 20, col = "red")
      # }
      # t_seq$r[s] <- circle$circle[3]
      # t_seq$X[s] <- circle$circle[1]
      # t_seq$Y[s] <- circle$circle[2]
      rs <- c(rs, circle$circle[3])
      xs <- c(xs, circle$circle[1])
      ys <- c(ys, circle$circle[2])
    }

    t_seq$r <- rs
    t_seq$X <- xs
    t_seq$Y <- ys
    # remove outliers of t_seq$r based on the median and the 2*SD
    t_seq$r[abs(t_seq$r - median(t_seq$r, na.rm = TRUE)) > 2 * sd(t_seq$r, na.rm = TRUE)] <- NA

    # remove rows from the sequence were the circle was getting bigger > 1.3 times the smallest previous circle

    for(s in 2:nrow(t_seq)){
      if(na2true(t_seq$r[s] > suppressWarnings(min(t_seq$r[1:(s-1)], na.rm = TRUE) * 1.3 )) | na2true(t_seq$r[s] < 0.04) | all(is.na(t_seq$r[1:s]))){
        t_seq[s, c("X","Y","r")] <- NA
      }
    }
    #rm(na2true)

    # if there are less than 3 points in the sequence, the DBH is not estimated
    if((sum(!is.na(t_seq$r)) <= 3) | (sum(!is.na(t_seq$X)) <= 3) | (sum(!is.na(t_seq$Y)) <= 3)){
      dbh <- mean(t_seq$r, na.rm = TRUE) * 2
      # check if dbh is within limits
      if(is.na(dbh) | dbh > max_dbh | dbh < 0){
        return(data.frame(TreeID = t, X = mean(tree$X), Y = mean(tree$Y), DBH = NA, quality_flag = 2))
      } else {
        return(data.frame(TreeID = t, X = mean(tree$X), Y = mean(tree$Y), DBH = dbh, quality_flag = 2))
      }
    }

    # fit 3 splines to the radius and displacement
    spline_r <- with(t_seq[!is.na(t_seq$r),], smooth.spline(Zmin, r, df = 3)) |> suppressWarnings()
    spline_x <- with(t_seq[!is.na(t_seq$X),], smooth.spline(Zmin, X, df = 20)) |> suppressWarnings()
    spline_y <- with(t_seq[!is.na(t_seq$Y),], smooth.spline(Zmin, Y, df = 20)) |> suppressWarnings()

    # predict the radius and displacement at 1.3m
    r <- predict(spline_r, 1.3)$y |> as.numeric()
    x <- predict(spline_x, 1.3)$y |> as.numeric()
    y <- predict(spline_y, 1.3)$y |> as.numeric()

    dbh <- r * 2
    # if the DBH is bigger than the maximum allowed or negative, the DBH is not estimated
    if(is.na(dbh) | dbh > max_dbh | dbh < 0){
      return(data.frame(TreeID = t, X = mean(tree$X), Y = mean(tree$Y),  DBH = NA, quality_flag = 3))
    } else {
      return(data.frame(TreeID = t, X = x, Y = y,  DBH = dbh, quality_flag = 4))
    }
    # update progress bar
    #setTxtProgressBar(pb, i)
    #i <- i + 1
  }

  parallel::stopCluster(cl)
  # close progress bar
  #close(pb)

  # Calculate the mean tree height above ground for dbh visualization
  dbh_slice <- lidR::add_attribute(dbh_slice, dbh_slice@data$Z - dbh_slice@data$Znorm, "Zdiff")
  tree_pos_height <- aggregate(dbh_slice$Zdiff, by = list(dbh_slice$TreeID), FUN = mean)
  names(tree_pos_height) <- c("TreeID", "Z")

  # calculate the tree heights by aggregating the original las file
  heights <- aggregate(tls@data$Z, by = list(tls@data$TreeID), FUN = function(x) max(x) - min(x))
  names(heights) <- c("TreeID", "Height")

  print("Tree heights calculated.")

  # add convex hull area to the results
  convhull_area <- function(xy){
    xy <- xy |> as.data.frame()
    if(nrow(xy) < 3){
      return(NA)
    }
    ch <- chull(xy)
    return(abs(0.5 * sum(xy[ch,1] * c(tail(xy[ch,2], -1), head(xy[ch,2], 1)) - c(tail(xy[ch,1], -1), head(xy[ch,1], 1)) * xy[ch,2])))
  }

  # initialize progress bar
  print("Calculating convex hull areas.")
  pb = txtProgressBar(min = 0, max = length(IDs), initial = 0, style = 3)
  i <- 1

  cpa <- data.frame(TreeID = numeric(), ConvexHullArea = numeric())
  for(t in IDs){
    tree <- tls |> lidR::filter_poi(TreeID == t & Znorm > 0.5)
    if(nrow(tree) < 3){
      cpa <- rbind(cpa, data.frame(TreeID = t, ConvexHullArea = NA))
      next
    }
    cpa <- rbind(cpa, data.frame(TreeID = t, ConvexHullArea = convhull_area(tree@data[, c("X", "Y")])))
    # update progress bar
    setTxtProgressBar(pb, i)
    i <- i + 1
  }

  # merge the results with the heights
  dbh_results <- merge(dbh_results, heights, by = "TreeID")
  dbh_results <- merge(dbh_results, cpa, by = "TreeID")
  dbh_results <- merge(dbh_results, tree_pos_height, by = "TreeID")
  return(dbh_results)
}

#' Function to plot the inventory results into a lidR 3d plot of the point cloud
#' @param plot lidR 3d plot
#' @param inventory data.frame with the inventory results
#' @param cex numeric size of the labels
#' @param label_col character color of the labels
#' @param col color of the spheres
#' @return the plot with the inventory results
#' @examples
#' \dontrun{
#' x <- lidR::plot(segmented, color = "TreeID")
#' plot_inventory(x, inventory)
#' }
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

