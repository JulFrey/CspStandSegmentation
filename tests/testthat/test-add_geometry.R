las <- suppressMessages(suppressWarnings(lidR::LAS(data.frame(X=runif(10),Y=runif(10),Z=runif(10)))))

testthat::test_that("errors", {
  testthat::expect_error(
    add_geometry(1,1),
    "las has to be a LAS object."
  )

  testthat::expect_error(
    add_geometry(las, -1),
    "k has to be one positive integer."
  )
})

testthat::test_that("output type", {
  testthat::expect_s4_class(
    add_geometry(las, 1),
    "LAS"
  )
})
