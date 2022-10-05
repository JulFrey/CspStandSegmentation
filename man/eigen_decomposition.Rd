\name{eigen_decomposition}
\alias{eigen_decomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast Eigenvalues decomposition for k nearest neighbours using a C++ function
%%  ~~function to do ... ~~
}
\description{
C++ helper function to compute eigenvalues for geometric feature calculation.
}
\usage{
eigen_decomposition(las, k, ncpu)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{LAS element}
  \item{k}{k nearest neighbours}
  \item{ncpu}{number of cpu cores to use}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function resturns for every point the 3 eigenvalues and the third element of the third eigenvector. These values are needed to compute planarity, linerity, verticality etc. in the add_geometry function
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
\link{add_geometry}}
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
