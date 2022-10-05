\name{comparative_shortest_path}
\alias{comparative_shortest_path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{helper function for csp_cost_segemntation
%%  ~~function to do ... ~~
}
\description{
The function performs a Dijkstra algorithem on a 3D voxel file to assign every voxel to the closest seed point using the igraph package.
}
\usage{
comparative_shortest_path(vox, adjacency_df, seeds, v_w = 0,l_w = 0,s_w = 0,
N_cores = parallel::detectCores()-1, Voxel_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vox}{a LAS S4 element with XYZ voxel coordinates in the @data slot.}
  \item{adjacency_df}{a data.frame with voxel ids (row numbers) in the first collumn and a neighbouring voxel id in the second collumn and the weight (distance) in the third collumn. Might be generated using the dbscan::frNN function (which requires reshaping the data).}
    \item{seeds}{seed points for tree positions.}
  \item{v_w,l_w,s_w}{wheights for verticality, linerity spericity see \code{\link{csp_cost_segmentation}} }
  \item{N_cores}{Number of cpu cores for multithreading}
  \item{Voxel_size}{Edge length used to create the voxels. This is only important to gain compareable distance weights on different voxel sizes. Should be greater than 0.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
\code{\link{csp_cost_segmentation}}
}
\examples{
%% ~~no examples~~
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
