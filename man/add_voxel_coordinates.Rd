\name{add_voxel_coordinates}
\alias{add_voxel_coordinates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add voxel coordinates to a las file
%%  ~~function to do ... ~~
}
\description{
Adds the collums x_vox, y_vox and z_vox in the given ressolution to the las element. This is convinient if informations have been derived in voxel space and these should be attached to the original points.
}
\usage{
add_voxel_coordinates(las, res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{an element of lidR::LAS class}
  \item{res}{voxel ressolution in [m]}
}
\details{
Voxel coordinates derived with this function are identical as those derived by lidR::voxelize.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
las <- lidR::readLAS(LASfile, select = "xyz", filter = "-inside 481250 3812980 481300 3813030")

las <- add_voxel_coordinates(las,res = 1)

lidR::plot(las, color = 'z_vox')
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
