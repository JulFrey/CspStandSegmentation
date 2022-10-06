\name{find_base_coordinates_raster}
\alias{find_base_coordinates_raster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find stem base position using a density raster approach
%%  ~~function to do ... ~~
}
\description{

}
\usage{
find_base_coordinates_raster(las, zmin = 0.5, zmax = 2, res = 0.1, quantile = 0.8, merge_radius = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{an element of lidR::LAS class}
  \item{zmin}{lower search boundary}
  \item{zmax}{upper search boundary}
  \item{res}{raster ressolution}
  \item{quantile}{raster density quantile to assign a tree region}
  \item{merge_radius}{search radius to merge base points}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

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
