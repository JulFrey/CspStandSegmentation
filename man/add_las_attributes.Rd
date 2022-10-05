\name{add_las_attributes}
\alias{add_las_attributes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Add all las_attributes from las@data to the header of a las element
%%  ~~function to do ... ~~
}
\description{
The helperfunction adds all headings from las@data which are nor part of lidR:::LASATTRIBUTES to the las header using lidR::add_lasattribute. Only attributes that are included in the header got saved when using lidR::writeLAS, this is a convinient way to add them.
}
\usage{
add_las_attributes(las)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{las}{an element of lidR::LAS class}
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
LASfile <- system.file("extdata", "MixedConifer.laz", package="lidR")
las <- lidR::readLAS(LASfile, select = "xyz", filter = "-inside 481250 3812980 481300 3813030")

las@data$noise <- runif(nrow(las@data))
las@data$noiseZ <- las@data$var1 * las@data$Z

las <- add_las_attributes(las)
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
