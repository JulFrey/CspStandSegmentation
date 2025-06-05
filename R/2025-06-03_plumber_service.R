# generate a plumber service to run segmentation and/or forest inventory from docker
# based on the example:
# tls = lidR::readTLSLAS(file)
#
# # find tree positions as starting points for segmentation
# map <- CspStandSegmentation::find_base_coordinates_raster(tls)
#
# # segment trees
# segmented <- tls |>
#   CspStandSegmentation::add_geometry(n_cores = parallel::detectCores()/2) |>
#   CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = parallel::detectCores()/2)
#
# # show results
# lidR::plot(segmented, color = "TreeID")
#
# # create inventory
# inventory <- CspStandSegmentation::forest_inventory(segmented)

library(plumber)

#* @apiTitle CspStandSegmentation
#* @post /segmentation_inv
#* @param file path to the las file
#* @param n_cores number of cores to use for processing
#* @post /process
#*
function(file, n_cores = parallel::detectCores()/2) {
  library(lidR)
  library(CspStandSegmentation)
  # print cores

  print(paste("Using", n_cores, "cores for processing."))

  # Convert n_cores safely
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() / 2
  } else {
    n_cores <- as.numeric(n_cores)
    if (is.na(n_cores) || n_cores <= 0) {
      stop("Invalid value for n_cores.")
    }
  }

  # Read the LAS file
  print(paste("Reading file:", file))
  tls <- lidR::readTLSLAS(file)

  # Find tree positions as starting points for segmentation
  map <- CspStandSegmentation::find_base_coordinates_raster(tls)
  if (is.null(map) || nrow(map) == 0) {
    stop("No valid tree positions found in the input file.")
  }
  # Segment trees
  segmented <- tls |>
    CspStandSegmentation::add_geometry(n_cores = n_cores) |>
    CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = n_cores)

  # Create inventory
  inventory <- CspStandSegmentation::forest_inventory(segmented)

  # ensure results folder exists
  if (!dir.exists("results")) {
    dir.create("results")
  }
  # Return the segmented data and inventory
  write.csv(inventory,"results/inventory.csv", row.names = F)
  lidR::writeLAS(segmented, "results/segmented.las")
  return(list(
    segmented = "results/segmented.las",
    inventory = "results/inventory.csv"
  ))
}

#* @apiTitle CspStandSegmentation
#* @post /segmentation
#* @param file path to the las file
#* @param n_cores number of cores to use for processing
#* @post /process
#*
function(file, n_cores = parallel::detectCores()/2) {
  library(lidR)
  library(CspStandSegmentation)

  # Convert n_cores safely
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() / 2
  } else {
    n_cores <- as.numeric(n_cores)
    if (is.na(n_cores) || n_cores <= 0) {
      stop("Invalid value for n_cores.")
    }
  }

  # Read the LAS file
  print(paste("Reading file:", file))
  tls <- lidR::readTLSLAS(file)

  # Find tree positions as starting points for segmentation
  map <- CspStandSegmentation::find_base_coordinates_raster(tls)
  if (is.null(map) || nrow(map) == 0) {
    stop("No valid tree positions found in the input file.")
  }
  # Segment trees
  segmented <- tls |>
    CspStandSegmentation::add_geometry(n_cores = n_cores) |>
    CspStandSegmentation::csp_cost_segmentation(map, 1, N_cores = n_cores)

  # ensure results folder exists
  if (!dir.exists("results")) {
    dir.create("results")
  }

  lidR::writeLAS(segmented, "results/segmented.las")
  return(list(
    segmented = "results/segmented.las"
  ))
}

#* @apiTitle CspStandSegmentation
#* @post /inventory
#* @param file path to the las file
#* @param id_name ExtraByte name for the TreeID default TreeID
#* @post /process
#*
function(file, id_name = "TreeID") {
  library(lidR)
  library(CspStandSegmentation)

  # check if TreeID in las
  header <- rlas::read.lasheader(file)
  if (!id_name %in% names(header$`Variable Length Records`$Extra_Bytes$`Extra Bytes Description`)) {
    stop("The LAS file does not contain the id_name attribute.")
  }

  # Read the LAS file
  print(paste("Reading file:", file))
  tls <- lidR::readTLSLAS(file)

  # rename the ids to TreeID
  if (id_name != "TreeID") {
    tls <- lidR::add_lasattribute(tls, tls[id_name], "TreeID", "single instance tree ids")
  }

  # Create inventory
  inventory <- CspStandSegmentation::forest_inventory(tls)

  # ensure results folder exists
  if (!dir.exists("results")) {
    dir.create("results")
  }

  write.csv(inventory, "results/inventory.csv", row.names = FALSE)
  return(list(
    inventory = "results/inventory.csv"
  ))
}

