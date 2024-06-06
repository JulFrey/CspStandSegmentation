las <- las <- suppressMessages(suppressWarnings(lidR::LAS(data.frame(X=runif(10),Y=runif(10),Z=runif(10)))))
las@data$test <- 1

testthat::test_that("errors", {
  testthat::expect_error(
    add_las_attributes(1),
    "las has to be a LAS object."
  )
})

testthat::test_that("output type", {
  testthat::expect_s4_class(
    add_las_attributes(las),
    "LAS"
  )
  testthat::expect_equal(
    length(add_las_attributes(las)@header@VLR),
    1
  )
})
