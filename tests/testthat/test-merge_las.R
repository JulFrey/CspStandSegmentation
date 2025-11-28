test_that("las_merge handles various input scenarios", {

  .LAS_silent <- function(x){
    las <- suppressMessages(suppressWarnings(lidR::LAS(x)))
  }

  # Load required package
  skip_if_not_installed("lidR")

  n <- 100

  # Create test data
  las1 <- .LAS_silent(data.frame(X = runif(n),
                               Y = runif(n),
                               Z = runif(n),
                               Intensity = rep(1L, n)))

  las2 <- .LAS_silent(data.frame(X = runif(n) + 2,
                               Y = runif(n),
                               Z = runif(n),
                               Intensity = rep(2L, n),
                               Classification = rep(1L, n)))

  las3 <- .LAS_silent(data.frame(X = runif(n) + 4,
                               Y = runif(n),
                               Z = runif(n),
                               Intensity = rep(3L, n),
                               ReturnNumber = rep(2L, n)))

  lidR::st_crs(las1) <- 4326
  lidR::st_crs(las2) <- 4326
  lidR::st_crs(las3) <- 4326

  lasList <- list(las1, las2, las3)

  # Test 1: Basic merging without OCI and fill=FALSE
  test_result <- las_merge(las1, las2, las3, oci = FALSE, fill = FALSE)
  expect_true(lidR::is(test_result, "LAS"))
  expect_equal(lidR::npoints(test_result), n * 3)
  expect_equal(colnames(test_result@data), c("X", "Y", "Z", "Intensity"))

  # Test 2: Basic merging with OCI
  test_result_oci <- las_merge(las1, las2, las3, oci = TRUE, fill = FALSE)
  expect_true("oci" %in% colnames(test_result_oci@data))
  expect_true(all(test_result_oci@data$oci[1:n] == 1))
  expect_true(all(test_result_oci@data$oci[(n + 1):(2 * n)] == 2))

  # Test 3: merging with list input
  test_result_list <- las_merge(lasList, oci = FALSE, fill = FALSE)
  expect_equal(lidR::npoints(test_result_list), n * 3)

  # Test 4: fill=TRUE with different columns
  test_result_fill <- las_merge(las1, las2, las3, oci = TRUE, fill = TRUE)
  expected_cols <- c("X", "Y", "Z", "Intensity", "Classification", "ReturnNumber", "oci")
  expect_true(all(expected_cols %in% colnames(test_result_fill@data)))
  expect_true(all(is.na(test_result_fill@data$Classification[1:n])))

  # Test 5: Error handling for non-LAS input
  expect_error(las_merge(las1, data.frame(X=1:10), las3),
               "All inputs must be LAS objects")

  # Test 6: Error handling for different CRS
  las4 <- .LAS_silent(data.frame(X = runif(n) + 4,
                                 Y = runif(n),
                                 Z = runif(n),
                                 Intensity = rep(3L, n),
                                 ReturnNumber = rep(2L, n)))
  expect_warning(las_merge(las4, las4),
                 "Some inputs do not have a CRS")

  # Temporarily modify CRS of one LAS
  original_crs <- lidR::st_crs(las3)
  lidR::st_crs(las3) <- 32632

  expect_error(las_merge(las1, las3),
               "All inputs must have the same CRS")


})

test_that("las_merge maintains data integrity", {

  skip_if_not_installed("lidR")

  .LAS_silent <- function(x){
    las <- suppressMessages(suppressWarnings(lidR::LAS(x)))
  }

  n <- 150
  # Create LAS with unique data in each column
  las1 <- .LAS_silent(data.frame(
    X = runif(n),
    Y = runif(n),
    Z = runif(n),
    Color = rep(1, n),
    Classification = rep(2L, n)
  ))

  las2 <- .LAS_silent(data.frame(
    X = runif(n, min = 1, max = 2),
    Y = runif(n),
    Z = runif(n),
    Intensity = rep(3L, n)
  ))

  lidR::st_crs(las1) <- 4326
  lidR::st_crs(las2) <- 4326

  # Test data integrity with fill=FALSE
  merged <- las_merge(las1, las2, fill = FALSE)
  expect_true(all(merged@data$X[1:n] %in% las1@data$X))
  expect_true(all(merged@data$Y[1:n] %in% las1@data$Y))

  # Test data integrity with fill=TRUE
  merged_full <- las_merge(las1, las2, fill = TRUE)
  expect_true(all(merged_full@data$Classification[1:n] == 2))
  expect_true(all(merged_full@data$Intensity[(n+1):(2*n)] == 3))
  expect_true(all(is.na(merged_full@data$Classification[(n+1):(2*n)])))
})

test_that("las_merge edge case", {
  skip_if_not_installed("lidR")

  .LAS_silent <- function(x){
    las <- suppressMessages(suppressWarnings(lidR::LAS(x)))
  }

  # Test with zero points
  las_empty <- .LAS_silent(data.frame(X = numeric(), Y = numeric(),Z = numeric() ))
  las1 <- .LAS_silent(data.frame(X = 1:5, Y = 1:5, Z = 1:5) * 1.1)
  lidR::st_crs(las_empty) <- 4326
  lidR::st_crs(las1) <- 4326

  result <- las_merge(las_empty, las1, fill = FALSE)
  expect_equal(lidR::npoints(result), 5)
})
