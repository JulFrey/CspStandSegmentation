\name{voxelize_points_mean_attributes}
\alias{voxelize_points_mean_attributes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{helper function to voxelize a las element
%%  ~~function to do ... ~~
}
\description{
Calculate voxel mean values for all numeric attributes in the las@data table including the XYZ-coordinates.
}
\usage{
voxelize_points_mean_attributes(las, res)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{a lidR::LAS element}
  \item{res}{voxel resolution in meter}
}
\details{
Returns a las element with XYZ-coordinates as the voxel center and X_gr,Y_gr,Z_gr as the center of gravity (mean point coordinates) as well as all other numeric collumns voxel mean  values with their original name.
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
\code{\link{voxelize_points}}
}
\examples{
# read example data
file = system.file("extdata", "pine_plot.laz", package="TreeLS")
tls = lidR::readTLSLAS(file)
tls |> voxelize_points_mean_attributes(1) |> lidR::plot(color = 'X_gr')
}
