#' Point distance function
#'
#' calculates eucleadian distances for n dimensions
#'
#' @param p1 point 1
#' @param p2 point 2
#'
#' @return the distance between the two points
#' @export p_dist
#'
#' @examples
#' dist(c(0,0), c(3,4))
p_dist <- function(p1, p2){
  if(length(p1) != length(p2)){
    stop("p1 and p2 must have the same length")
  }
  sqrt(sum((p1 - p2)^2))
}

#' Farthest Distance Sampling (Farthest Point Sampling)
#'
#' This function selects n points from a matrix of points such that the minimum distance between any two points is maximized.
#' This version is memory efficient and can handle large matrices.
#'
#' @param mat a matrix of points with one row for each point and one collumn for each dimension
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
#' \dontrun{
#'   plot(mat, col = "black", pch = 19)
#'   points(sample, col = "red", pch = 19)
#' }
fds <- function(mat, n, ret = "idx", scale = F){
  # check the inputs
  if(ret != "idx" & ret != "mat"){
    stop("ret must be 'idx' or 'mat'")
  }
  if(n >= nrow(mat)){
    warning("n is greater or equal than the number of points in mat. Returning all points.")
    if(ret == "mat"){
      return(mat)
    }
    return(1:nrow(mat))
  }
  if(n == 1){
    return(sample(1:nrow(mat), 1))
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
  idx <- sample(1:nrow(mat), 1)
  # calculate a vector of distances from the first point
  dists <- apply(mat, 1, function(x) p_dist(mat[idx,], x))

  # select all further points in a loop
  for(i in 2:n){
    # select the next point
    idx <- c(idx, which.max(dists))
    # calculate the distances from the new point
    dists <- pmin(dists, apply(mat, 1, function(x) p_dist(mat[idx[i],], x)))
  }

  # return the selected points
  if(ret == "mat"){
    return(mat[idx,])
  } else {
    return(idx)
  }
}

