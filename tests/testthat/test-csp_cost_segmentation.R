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

testthat::test_that("multi-source routing matches repeated Dijkstra", {
  edges <- data.frame(
    from = c(1L, 2L, 3L, 4L, 1L),
    to = c(2L, 3L, 4L, 5L, 5L),
    weight = c(1, 1, 1, 1, 10)
  )
  graph <- igraph::graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = as.character(1:6)
  )
  seeds <- c(1L, 5L)
  expected_distances <- igraph::distances(
    graph,
    v = as.character(seeds),
    algorithm = "dijkstra"
  )
  expected_seed <- max.col(-t(expected_distances), ties.method = "first")
  reachable <- is.finite(expected_distances[cbind(expected_seed, seq_len(igraph::vcount(graph)))])
  endpoints <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  actual <- multi_source_dijkstra(
    endpoints[, 1],
    endpoints[, 2],
    igraph::E(graph)$weight,
    igraph::vcount(graph),
    seeds
  )

  testthat::expect_identical(actual$seed_index[reachable], expected_seed[reachable])
  testthat::expect_equal(
    actual$distance,
    expected_distances[cbind(expected_seed, seq_len(igraph::vcount(graph)))]
  )
  testthat::expect_true(is.infinite(actual$distance[6]))
  testthat::expect_true(is.na(actual$seed_index[6]))
})

testthat::test_that("explicit seed links are free and survive edge simplification", {
  vox <- suppressMessages(suppressWarnings(lidR::LAS(data.frame(
    X = c(0, 1, 0.5), Y = 0, Z = 0,
    Verticality = 0, Sphericity = 0, Linearity = 0, TreeID = 0
  ))))
  adjacency <- data.frame(
    adjacency_list_id = c(1L, 2L, 3L),
    adjacency_list = c(2L, 1L, 1L),
    weight = c(1, 1, 0),
    seed_link = c(FALSE, FALSE, TRUE)
  )
  output <- comparative_shortest_path(
    vox = vox,
    adjacency_df = adjacency,
    seeds = data.frame(SeedID = 3L, TreeID = 99L),
    Voxel_size = 1,
    N_trees = 1
  )

  testthat::expect_equal(output@data$TreeID[1:2], c(99L, 99L))
  testthat::expect_equal(output@data$dist[1:2], c(0, 1))
})

testthat::test_that("negative routing costs are clamped to zero", {
  vox <- suppressMessages(suppressWarnings(lidR::LAS(data.frame(
    X = c(0, 1, 0.5), Y = 0, Z = 0,
    Verticality = 1, Sphericity = 0, Linearity = 0, TreeID = 0
  ))))
  adjacency <- data.frame(
    adjacency_list_id = c(1L, 2L, 3L),
    adjacency_list = c(2L, 1L, 1L),
    weight = c(0.1, 0.1, 0),
    seed_link = c(FALSE, FALSE, TRUE)
  )
  output <- comparative_shortest_path(
    vox = vox,
    adjacency_df = adjacency,
    seeds = data.frame(SeedID = 3L, TreeID = 99L),
    v_w = -1,
    Voxel_size = 1,
    N_trees = 1
  )

  testthat::expect_true(all(output@data$dist >= 0))
})
