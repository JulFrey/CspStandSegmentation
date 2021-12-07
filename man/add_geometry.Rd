\name{add_geometry}
\alias{add_geometry}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add geometric features to a LAS object
}
\description{
The function calls a fast cpp multicore function to calculate eigenvalues for the points in a point cloud based on the k neerest neighbors. Afterwards it adds geometric features like Curvature, Linearity, Planarity, Sphericity, Anisotrophy and Verticlity to the points itself.}
\usage{
add_geometry(las, k , n_cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{A LAS object (see lidR::LAS)}
  \item{k}{the k neerest neighbors to use for the eigenvalue calculation}
  \item{n_cores}{The number of CPU cores to use}
%%     ~~Describe \code{x} here~~
}

\details{
Details to the metrics can be found in:
\\
Hackel, T., Wegner, J.D. & Schindler, K. (2016) Contour Detection in Unstructured 3D Point Clouds. In 2016 IEEE Conference on Computer Vision and Pattern Recognition (CVPR). Presented at the 2016 IEEE Conference on Computer Vision and Pattern Recognition (CVPR), IEEE, Las Vegas, NV, USA, pp. 1610–1618.

}
\value{
The function returns a single LAS object with the geometric features attached to it in the LAS@data section. 
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

las <- add_geometry(las, k = 5, n_cores = parallel::detectCores()-1)
summary(las@data)

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
