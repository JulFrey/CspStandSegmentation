\name{csp_cost_segmentation}
\alias{csp_cost_segmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparative Shortest Path with cost weighting tree segementation function
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
  \item{S_w}{Linearity weight. Simmilar to V_w but for sphericity, higher values indicate a malus for spherical shapes (usually small branches and leaves).}
  \item{N_Cores}{number of CPU cores used for parrallel routing using the foreach package.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a copy of the las point cloud with an additional field for the TreeID.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dr. Julian Frey
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
segmented <- tls %>%
  lidR::filter_poi(Classification != 2) %>%
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
