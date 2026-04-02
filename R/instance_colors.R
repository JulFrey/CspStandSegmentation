################################################################################
# PRETTY COLORING BASED ON GRAPH COLORING ALGORITHM
################################################################################

# HELPERS ======================================================================

#' Return a palette-generating function
#'
#' get_pal returns a function that produces a character vector of hex colours
#' given a required length. It exposes a set of named palettes used by
#' color_ids().
#'
#' @param name Character scalar. Name of the palette to use. Supported names:
#'   "sky", "sea", "cozy", "fairy", "winter", "rainbow", "pastel", "candy",
#'   "boring". Default: "rainbow".
#'
#' @return A function with signature \code{function(n)} that returns \code{n}
#'   hex colour strings.
#'
#' @examples
#' get_pal("sky")
#'
#' @export
get_pal <- function(name = "rainbow") {
  # modified from coolors.co
  palettes <- list(
    "sky" = c("#f2c13a","#ef7239","#f23184","#955ae8","#5c96f5"),
    "sea" = c("#d9ed92","#99d98c","#52b69a","#168aad","#1e6091","#154366","#0C304C"),
    "cozy" = c("#41764c","#a7c957","#eadbb3","#d68b8b","#bc4749"),
    "fairy" = c("#cd92ef","#ffadcb","#fed2e2","#b1e2fc","#9bc1ff"),
    "winter" = c("#007DA3","#00afb9","#fdfcdc","#fdca9b","#ee6258"),
    "rainbow" = c("#f52e2e","#f5982e","#f5d42e","#dcf636","#acf636","#36f6a6","#36e9f6","#3d90ee","#7336f6","#c236f6"),
    "pastel" = c("#fcf7b7","#fed9bb","#ffbfc3","#f7aae7","#c4a6f5","#90bbf8","#7cd9f8","#7beff9","#85f9e0","#a4fdac"),
    "candy" = c("#9137ff","#ff5ce4","#ff4545","#fee440","#00bbf9","#00f5bc"),
    "boring" = c("#DEE2E6", "#191C1F")
  )
  return(colorRampPalette(palettes[[name]]))
}


# FUNCTIONS ====================================================================
#' Assign colours to instances in a LAS object using a graph-colouring heuristic
#'
#' color_ids computes a per-instance colour assignment and merges RGB (R,G,B)
#' and a numeric \code{color_ID} back into \code{las@data}. Colour assignment is
#' driven by spatial adjacency (k-nearest neighbours on instance centroids)
#' and maximizing colour LAB distance from neighbouring instances.
#'
#' @param las A \code{LAS} object (from the \pkg{lidR} package).
#' @param col Either a character string naming a palette supported by
#'   \code{get_pal()} or a character vector of hex colours. Default: "sky".
#' @param n_col Integer. Number of discrete colours to generate from the palette.
#'   Default: 10.
#' @param n_neighbors Integer. Number of nearest neighbours to treat as
#'   adjacency when assigning colours. Default: 10.
#' @param instance_id Character. Column name in \code{las@data} that identifies
#'   instances (for example tree or stand IDs). Default: "ID".
#' @param ground_id Value used to identify ground or background instance;
#'   assigned \code{ground_color}. Default: 0.
#' @param ground_color Character. Hex colour used for the ground instance.
#'   Default: "#ffffff".
#' @param overwrite_rgb Logical. If \code{TRUE} existing R/G/B columns in
#'   \code{las@data} will be overwritten. Default: \code{TRUE}.
#'
#' @return A \code{LAS} object with R, G, B (0-255 integers) and \code{color_ID}
#'   merged into \code{las@data}. If \code{overwrite_rgb} is \code{FALSE} and
#'   the LAS already had RGB columns, those channels are preserved.
#'
#' @details
#' The function computes instance centroids, builds a k-nearest-neighbour
#' graph, sorts instances by local neighbour density and assigns discrete
#' colours greedily so neighbouring instances receive maximally different
#' colours in LAB space.
#'
#' @examples
#' \dontrun{
#' # read example data
#' \donttest{
#' file = system.file("extdata", "beech.las", package="CspStandSegmentation")
#' las = lidR::readTLSLAS(file)
#'
#' # Find tree positions as starting points for segmentation
#' map <- CspStandSegmentation::find_base_coordinates_raster(las)
#'
#' # segment trees
#' segmented <- las |>
#' CspStandSegmentation::csp_cost_segmentation(map, 1, S_w = 0.5)
#' las_col <- color_ids(las, col = "sky", instance_id = "PredInstance")
#' lidR::plot(las_col, color = "RGB")
#' }
#' }
#'
#' @importFrom colorspace hex2RGB
#' @importFrom RANN nn2
#' @importFrom data.table data.table
#' @export
color_ids <- function(
    las,
    col = "sky",
    n_col = 10,
    n_neighbors = 10,
    instance_id = "ID",
    ground_id = 0,
    ground_color = "#ffffff",
    overwrite_rgb = TRUE) {
  
  # check inputs
  if(!lidR::is(las, "LAS")) stop("las must be a LAS object")
  
  had_rgb <- all(c("R", "G", "B") %in% names(las@data))
  # check if las has rgb values and warn if they will be overwritten
  if(had_rgb & overwrite_rgb) {
    warning("LAS already has RGB values, they will be overwritten")
  }

  df <- las@data

  # strip off any existing RGB columns
  if(overwrite_rgb) {
    df <- df[, !names(df) %in% c("R", "G", "B")]
  }

  # assign color palette
  pal_fun <- ifelse(length(col) > 1, colorRampPalette(col), get_pal(col))
  pal <- pal_fun(n_col)
  
  # rename columns
  names(df)[names(df) == instance_id] <- "ID"
  
  # calculate ID centers
  centers <- df[ID != ground_id, .(X = mean(X), Y = mean(Y), Z = min(Z)), by = ID]
  
  # create color lookup
  color_RGB <- colorspace::hex2RGB(pal)@coords * 255
  color_LAB <- as(colorspace::hex2RGB(pal), "LAB")
  color_lookup <- data.table::data.table(
    "color_ID" = 1:n_col,
    "color_HEX" = pal,
    "R" = as.integer(color_RGB[,1]),
    "G" = as.integer(color_RGB[,2]),
    "B" = as.integer(color_RGB[,3])
  )
  
  # calculate LAB color distances (matrix)
  color_distances <- as.matrix(dist(color_LAB@coords))
  
  # calculate spatial distances (graphs)
  kdtree <- RANN::nn2(centers[, .(X, Y, Z)], k = n_neighbors)
  
  # get mean distance of n nearest neighbours per instance
  mean_distances <- rowMeans(kdtree$nn.dists)
  
  # sort IDs depending on average distance to neighbours
  # smallest distance first -> most close neighbours first
  sorted_center_ids <- centers[order(mean_distances), ID]
  
  # prepare output data
  out <- data.table::data.table(
    "ID" = centers$ID,
    "color_ID" = 0
  )
  
  # loop through all IDs
  for (curr_center_id in sorted_center_ids) {
    
    # extract nearest neighbours (or: within a certain distance?)
    kdd_idx <- which(centers$ID == curr_center_id)
    curr_nn_ids <- kdtree$nn.idx[kdd_idx,][-1]
    curr_nn_cols <- out$color_ID[out$ID %in% centers$ID[curr_nn_ids]]
    
    # check colors of nearest neighbours
    if (sum(curr_nn_cols) == 0) {
      
      # assign random color
      out[out$ID == curr_center_id,"color_ID"] <- sample(color_lookup[,color_ID], 1)
    } else {
      
      # get unused colors
      colors_unused <- setdiff(color_lookup$color_ID, unique(curr_nn_cols))
      colors_used <- setdiff(color_lookup$color_ID, colors_unused)
      
      # check if all colors have been used
      if (length(colors_unused) == 1) {
        
        # use left over color
        col_new <- colors_unused
        
      } else {
        
        # check if no colors left
        if (length(colors_unused) == 0) {
          
          # throw warning
          warning("not enough colors or too many neighbours")
          
          # set all colors to unused
          colors_unused <- color_lookup$color_ID
        }
        
        # get total color distances of unused to used colors
        if (length(colors_used) == 1) {
          col_dist_unused <- color_distances[colors_unused,colors_used]
        } else {
          col_dist_unused <- rowSums(color_distances[colors_unused,colors_used])
        }
        
        # get most different color
        col_new <- as.numeric(names(col_dist_unused[which.max(col_dist_unused)]))
      }
      
      # assign new color
      out[out$ID == curr_center_id,"color_ID"] <- col_new
    }
  }
  
  # assign RGB color to points
  out <- merge(out, color_lookup)
  
  # add ground color
  ground_rgb <- as.integer(colorspace::hex2RGB(ground_color)@coords * 255)
  out <- rbind(
    out,
    data.table::data.table(
      "color_ID" = 0,
      "ID" = ground_id,
      "color_HEX" = ground_color,
      "R" = ground_rgb[1],
      "G" = ground_rgb[2],
      "B" = ground_rgb[3]
      ))
  
  # rename column
  names(out)[names(out) == "ID"] <- instance_id
  
  # remove color_HEX for merging
  out_merge <- out[,.(get(instance_id), R, G, B, color_ID)]
  names(out_merge)[1] <- instance_id

  # remove RGB if overwrite_rgb is FALSE
  if(!overwrite_rgb & had_rgb) {
    out_merge <- out_merge[, !names(out_merge) %in% c("R", "G", "B")]
  }

  # merge into las data
  las@data <- merge(las@data, out_merge, by = instance_id, all.x = TRUE)
  return(las)
}


# read data
las <- lidR::readLAS("D:/github/trees.laz", select = "0")

# add color values
las <- color_ids(las, col = get_pal(), instance_id = "PredInstance", ground_color = "#ff0000")

# show results
lidR::plot(las, color = "RGB")

# ==============================================================================