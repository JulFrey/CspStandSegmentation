\name{csp_cost_segmentation}
\alias{csp_cost_segmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparative Shortest Path with cost weighting tree segementation
%%  ~~function to do ... ~~
}
\description{
Segments single trees from forest point clouds based on tree positions (xy-coordinates) provided in the map-argument.
}
\usage{
csp_cost_segmentation(las, map, Voxel_size = 0.3, V_w = 0,L_w = 0,S_w = 0, N_cores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{A lidR LAS S4 object.}
  \item{map}{A TreeLS map object, or a data.frame including the columns X,Y,Z,TreeID, with X and Y depicting the location of the trees.}
  \item{Voxel_size}{The voxel size (3D ressolution) for the routing graph to determine nearest map location for every point in the pointcloud.}
  \item{V_w}{verticality weight. Since trunks are vertical structures routing through voxels with high verticality can be rated 'cheaper'. should be a number between 0 and 1 with 0 meaning no benefit vor more vertical structures.}
  \item{L_w}{Linearity weight. Simmilar to V_w but for linearity, higher values indicate a malus for linear shapes  (usually branches).}
  \item{S_w}{Spericity weight. Simmilar to V_w but for sphericity, higher values indicate a malus for spherical shapes (usually small branches and leaves).}
  \item{N_cores}{number of CPU cores used for parrallel routing using the foreach package.}
}
\details{
The whole point cloud is voxelized in the given resolution and the centre of gravity for the points inside is calculated as voxel coordinate. A graph is build which connects the voxel-coordinates based on db-scan algorithm. The distances between the voxel-coordinates is weighted based on geometric features computed for the points in the voxels. Distances along planar and/or vertical faces like stems are weighted shorter than distances through voxels with a high sphericity like leaves and clusters of twigs. This avoids, that small trees and regrowth takes over points from mature trees. For every voxel-centre the weighted distance in the network is calculated to all tree-locations from the map-argument. The TreeID of the map argument with the shortest distance is assigned to the voxel. All points in the point cloud get the TreeID from their parent voxel.}
\value{
Returns a copy of the las point cloud with an additional field for the TreeID.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Julian Frey <julian.frey@iww.uni-freiburg.de>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{comparative_shortest_path}}
}
\examples{
# read example data
file = system.file("extdata", "pine_plot.laz", package="TreeLS")
tls = lidR::readTLSLAS(file)

# normalize height
tls <- TreeLS::tlsNormalize(tls)

# find tree positions as starting point for segmentation
map <- TreeLS::treeMap(tls)

# segment trees
segmented <- tls |>
  lidR::filter_poi(Classification != 2) |>
  add_geometry() |>
  csp_cost_segmentation(map, 1)

lidR::plot(segmented, color = "TreeID")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
