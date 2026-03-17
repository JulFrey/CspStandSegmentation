#' Point distance function
#'
#' calculates euclidean distances for n dimensions
#'
#' @param p1 point 1
#' @param p2 point 2
#'
#' @return the distance between the two points
#' @export p_dist
#'
#' @examples
#' p_dist(c(0,0), c(3,4))
p_dist <- function(p1, p2){
  if(length(p1) != length(p2)){
    stop("p1 and p2 must have the same length")
  }
  sqrt(sum((p1 - p2)^2))
}

#' Point distance function
#'
#' calculates euclidean distances for n dimensions between a matrix of points and a single point
#'
#' @param mat matrix with points as rows
#' @param p point to calculate distances
#'
#' @return the distances between every row of mat and p
#' @export p_mat_dist
#'
#' @examples
#' p_mat_dist(as.matrix(cbind(runif(100),runif(100))), c(3,4))
p_mat_dist <- function(mat, p){
  mat2 <- mat
  for(c in 1:ncol(mat)){
    mat2[,c] <- (mat[,c] - p[c])^2
  }
  return(sqrt(.rowSums(mat2, nrow(mat2), ncol(mat2))))
}

#' Farthest Distance Sampling (Farthest Point Sampling)
#'
#' This function selects n points from a matrix of points such that the minimum
#' distance between any two points is maximized. Either a fixed number of points
#' can be selected (using the `n` parameter) or points can be selected until a
#' specified minimum spacing is violated (using the `spacing` parameter).
#' This version is memory efficient and can handle large matrices.
#'
#' @param mat a matrix of points with one row for each point and one column for
#' each dimension, can also be a las object then only XYZ will be used
#' @param n the number of points to select, or if <1 the proportion of points to
#' select
#' @param spacing If supplied, sampling continues until the minimum nearest-
#' neighbour distance among selected points drops below `spacing`.
#' @param ret the type of output to return. Options are "idx" (default) to return
#' the indices of the selected points, "mat" to return the selected points.
#' @param scale logical. If TRUE, the dimensions are scaled to have a mean of 0
#' and a standard deviation of 1 before calculating distances.
#'
#' @return a vector of indices or a matrix of points
#' @export fds
#'
#' @examples
#' mat <- matrix(rnorm(1000), ncol = 10)
#' sample <- fds(mat, n = 50, ret = "mat")
#' str(sample)
#'
fds <- function(mat, n = NULL, spacing = NULL, ret = "idx", scale = FALSE){
  # check the inputs
  was_las <- FALSE
  if(!is.matrix(mat)){
    if(is(mat, "LAS")){
      was_las <- TRUE
      las <- mat
      mat <- as.matrix(las@data[,c("X", "Y", "Z")])
    } else {
      stop("mat must be a matrix or a LAS object")
    }
  }

  if(ret != "idx" & ret != "mat"){
    stop("ret must be 'idx' or 'mat'")
  }


  # n iteration handling
  # spacing-restricted sampling
  if(!is.null(spacing)){
    if(!is.numeric(spacing) || length(spacing) != 1L || spacing <= 0){
      stop("'spacing' must be a single positive number")
    }
    max_n <- nrow(mat)
  } else { # n-restricted sampling
    if(n >= nrow(mat)){
      warning("n is greater or equal than the number of points in mat. Returning all points.")
      if(ret == "mat"){
        if(was_las){
          return(las)
        }
        return(mat)
      }
      return(1:nrow(mat))
    }
    if(n == 1){
      if(was_las){
        return(las[sample(1:nrow(mat), 1),])
      }
      if(ret == "mat"){
        return(mat[sample(1:nrow(mat), 1),])
      } else {
        return(sample(1:nrow(mat), 1))
      }
    }
    if(n < 0){
      stop("n must be greater than 0.")
    }
    if(n < 1){
      n <- round(nrow(mat) * n)
    }
  }

  # scale dimensions if requested
  if(scale){
    for(c in 1:ncol(mat)){
      mat[,c] <- scale(mat[,c])
    }
  }

  # select the first point randomly
  idx <- which.max(mat[,1])

  # calculate a vector of distances from the first point
  dists <- p_mat_dist(mat, mat[idx,])

  # determine the number of iterations (either n or max_n)
  n_iter <- if(is.null(spacing)) n else max_n

  # select all further points in a loop
  for(i in 2:n_iter){
    # select the next point
    idx <- c(idx, which.max(dists))

    # calculate the distances from the new point
    dists2 <- p_mat_dist(mat, mat[idx[i],])
    dists <- pmin(dists, dists2)

    # spacing-based stopping rule (minimum nearest neighbour distance)
    if(!is.null(spacing)){
      if(length(idx) > 1){
        # compute nearest neighbour distances among selected points
        nn <- FNN::get.knn(mat[idx, , drop = FALSE], k = 1)
        # minimum nearest neighbor distance among selected points
        min_nn <- min(nn$nn.dist, na.rm = TRUE)

        # if the minimum nearest neighbour distance drops below spacing, stop
        # sampling
        if(min_nn < spacing){
          # remove last added point (it violated spacing)
          idx <- idx[-length(idx)]
          break
        }
      }
    }
  }

  # return the selected points
  if(ret == "mat"){
    if(was_las){
      las@data <- las@data[idx,]
      return(las)
    } else {
      return(mat[idx,])
    }
  } else {
    return(idx)
  }
}
