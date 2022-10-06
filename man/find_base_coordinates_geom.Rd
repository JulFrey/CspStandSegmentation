\name{find_base_coordinates_geom}
\alias{find_base_coordinates_geom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find stem base position using a geometric feature filtering and clustering approach
%%  ~~function to do ... ~~
}
\description{

}
\usage{
find_base_coordinates_geom(las, zmin = 0.5, zmax = 2, res = 0.5, min_verticality = 0.9, min_planarity = 0.5, min_cluster_size = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{an element of lidR::LAS class}
  \item{zmin}{lower search boundary}
  \item{zmax}{upper search boundary}
  \item{res}{cluster search radius}
  \item{min_verticality}{minimum verticality >0 & <1 for a point to be considered a stem point}
  \item{min_planarity}{minimum planarity >0 & <1 for a point to be considered a stem point}
  \item{min_cluster_size}{minimum number of points in cluster to be considered a tree, if NULL median cluster size is choosen}
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
