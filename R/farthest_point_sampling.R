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
#' This function selects n points from a matrix of points such that the minimum distance between any two points is maximized.
#' This version is memory efficient and can handle large matrices.
#'
#' @param mat a matrix of points with one row for each point and one column for each dimension, can also be a las object then only XYZ will be used
#' @param n the number of points to select, or if <1 the proportion of points to select
#' @param ret the type of output to return. Options are "idx" (default) to return the indices of the selected points, "mat" to return the selected points.
#' @param scale logical. If TRUE, the dimensions are scaled to have a mean of 0 and a standard deviation of 1 before calculating distances.
#'
#' @return a vector of indices or a matrix of points
#' @export fds
#'
#' @examples
#' mat <- matrix(rnorm(1000), ncol = 10)
#' sample <- fds(mat, 50, ret = "mat")
#' str(sample)
#'
#'  #plot(mat, col = "black", pch = 19)
#'  #points(sample, col = "red", pch = 19)
#'
fds <- function(mat, n, ret = "idx", scale = FALSE){
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

  # select all further points in a loop
  for(i in 2:n){
    # select the next point
    idx <- c(idx, which.max(dists))
    # calculate the distances from the new point
    dists2 <- p_mat_dist(mat, mat[idx[i],])
    dists <- pmin(dists, dists2)
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
