% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2021-11-12_A1_JF_csp_cost_functions.R
\name{add_voxel_coordinates}
\alias{add_voxel_coordinates}
\title{Add voxel coordinates to a las file}
\usage{
add_voxel_coordinates(las, res)
}
\arguments{
\item{las}{an element of lidR::LAS class}

\item{res}{voxel ressolution in [m]}
}
\value{
las file with additional voxel coordinates
}
\description{
Adds the collums x_vox, y_vox and z_vox in the given ressolution to the las
element. This is convenient if information has been derived in voxel space
and these should be attached to the original points.
}
\details{
Voxel coordinates derived with this function are identical to those derived
by lidR::voxelize.
}
\examples{

file = system.file("extdata", "beech.las", package="CspStandSegmentation")
las = lidR::readTLSLAS(file)

las <- add_voxel_coordinates(las,res = 1)

lidR::plot(las, color = 'z_vox')

}
\author{
Julian Frey <julian.frey@wwd.uni-freiburg.de>
}
