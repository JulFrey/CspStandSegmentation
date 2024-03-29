# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' helper function to unlist IDs generated by dbscan::frNN
#'
#' @param list a list element created by dbscan::frNN
#' @param l the expected length of the result
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#'   \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#'     'comp2'} %% ...
#'       @note %% ~~further notes~~
#'       @author Dr. Julian Frey <julian.frey@@iww.uni-freiburg.de>
#'       @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#'       @references %% ~put references to the literature/web site here ~
#'       @examples
#'
#' %%
#'
#'   @export fast_unlist
fast_unlist <- function(list, l) {
    .Call(`_CspStandSegmentation_fast_unlist`, list, l)
}

#' helper function to unlist distances computed by dbscan::frNN %% ~~function
#' to do ... ~~
#'
#'
#' %% ~~ If necessary, more details than the description above ~~
#'
#' @param list a list element created by dbscan::frNN
#' @param l the expected length of the result
#' @return %% ~Describe the value returned %% If it is a LIST, use %%
#' \item{comp1 }{Description of 'comp1'} %% \item{comp2 }{Description of
#' 'comp2'} %% ...
#' @note %% ~~further notes~~
#' @author Dr. Julian Frey <julian.frey@@iww.uni-freiburg.de>
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @examples
#'
#'
#'
#' @export fast_unlist_dist
fast_unlist_dist <- function(list, l) {
    .Call(`_CspStandSegmentation_fast_unlist_dist`, list, l)
}

#' Fast Eigenvalues decomposition for k nearest neighbors using a C++ function
#'
#' C++ helper function to compute eigenvalues for geometric feature
#' calculation.
#'
#' %% ~~ If necessary, more details than the description above ~~
#'
#' @param las LAS element
#' @param k k nearest neighbors
#' @param ncpu number of cpu cores to use
#' @return The function returns for every point the 3 eigenvalues and the
#' third element of the third eigenvector. These values are needed to compute
#' planarity, linerity, verticality etc. in the add_geometry function
#' @note %% ~~further notes~~
#' @author Dr. Julian Frey <julian.frey@@iww.uni-freiburg.de>
#' @seealso \link{add_geometry}
#' @references %% ~put references to the literature/web site here ~
#' @examples
#'
#'
#'
#' @export eigen_decomposition
eigen_decomposition <- function(las, k, ncpu = 1L) {
    .Call(`_CspStandSegmentation_eigen_decomposition`, las, k, ncpu)
}

