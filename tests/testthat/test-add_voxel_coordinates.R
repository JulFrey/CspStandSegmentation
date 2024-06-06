las <- suppressMessages(suppressWarnings(lidR::LAS(data.frame(X=runif(10),Y=runif(10),Z=runif(10)))))

testthat::test_that("errors", {
  testthat::expect_error(
    add_voxel_coordinates(1,1),
    "las has to be a LAS object."
  )

  testthat::expect_error(
    add_voxel_coordinates(las, -1),
    "res has to be numeric and positive."
  )
})

testthat::test_that("output type", {
  testthat::expect_s4_class(
    add_voxel_coordinates(las, 1),
    "LAS"
  )
  testthat::expect_equal(
    ncol(add_voxel_coordinates(las, 1)@data),
    6
  )
})
