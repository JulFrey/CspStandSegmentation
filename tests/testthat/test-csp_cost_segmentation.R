las <- suppressMessages(suppressWarnings(lidR::LAS(data.frame(X=runif(10),Y=runif(10),Z=runif(10)))))
map <- data.frame(X=0:1,Y=0:1,Z=0:1,TreeID = 1:2)

testthat::test_that("errors", {
  testthat::expect_error(
    csp_cost_segmentation(1,1),
    "las has to be a LAS object."
  )

  testthat::expect_error(
    csp_cost_segmentation(las, data.frame(X=1,Y=1,Z=1,Tree = 1)),
    "map has to be a data.frame with collumn names X,Y,Z,TreeID."
  )
  testthat::expect_error(
    csp_cost_segmentation(las, map, "a"),
    "Voxel_size, V_w, L_w, S_w and N_cores have to be numeric."
  )
})
