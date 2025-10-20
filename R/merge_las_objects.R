#' Makes one las-object from multiple las-objects
#'
#' This function merges multiple las objects into one las object. The function checks if all inputs are las-objects and if they have the same CRS. The function will also add a column oci with the original cloud index to each las-object. The function will then rbind all data by the minimum set of columns. If the fill argument is set to False, columns which do not exist in all las objects will be removed. If the fill argument is set to True, missing columns will be filled with NA.
#'
#'
#' @param ... any number of las objects
#' @param oci add a column with the original cloud index
#' @param fill fill missing columns with NA if it is set to False columns which do not exist in all las-objects will be removed
#'
#' @returns A single las-object with only the overlapping column name
#'
#' @examples
#' las1 <- lidR::LAS(data.frame(X = runif(100), Y = runif(100), Z = runif(100)))
#' las2 <- lidR::LAS(data.frame(X = runif(100) + 2, Y = runif(100), Z = runif(100)))
#' las3 <- lidR::LAS(data.frame(X = runif(100) + 4, Y = runif(100), Z = runif(100)))
#' merged <- las_merge(las1, las2, las3)
#' str(merged)
#' \donttest{
#' merged |> lidR::plot(color = "oci")
#' }
#'
#' @export las_merge
las_merge <- function(..., oci = TRUE, fill = FALSE){
  is_las <- function(x) lidR::is(x, "LAS")

  # check if only las objects are passed
  if(!all(sapply(list(...), is_las))) stop("All inputs must be LAS objects")

  # check if the crs is the same
  crs <- lapply(list(...), function(x) lidR::st_crs(x))
  wkt <- sapply(crs, function(x) x$wkt)
  if(any(is.na(wkt))) warning("Some inputs do not have a CRS")
  # check if all crs are the same
  if(length(unique(wkt)) > 1) stop("All inputs must have the same CRS")

  # put all objects in a list
  las_list <- list(...)
  # add original cloud index to each object
  if(oci){
    for(i in 1:length(las_list)){
      las_list[[i]] <- las_list[[i]] |> lidR::add_lasattribute(i, "oci", "original cloud index")
    }
  }

  # rbind all data by the minimum set of columns
  las_1 <- las_list[[1]]
  if(fill){
    las_1@data <- data.table::rbindlist(lapply(las_list, function(x) x@data), fill = TRUE)
  } else {
    # build a minimum set of columns
    cols <- lapply(las_list, function(x) colnames(x@data))
    # find the intersection of columns
    cols <- Reduce(intersect, cols)
    las_1@data <- do.call(rbind, lapply(las_list, function(x) x@data[,..cols]))
  }

  # quantize the data
  las_1 <- las_1 |> lidR::las_quantize() |> lidR::las_update()

  # remove extra byte arguments from header which do not exist any longer
  exbytes_header <- names(las_1@header@VLR$Extra_Bytes$`Extra Bytes Description`)
  exbytes_data <- names(las_1@data)
  exbytes <- setdiff(exbytes_header, exbytes_data)
  if(length(exbytes) > 0){
    for(i in exbytes){
      las_1@header@VLR$Extra_Bytes$`Extra Bytes Description`[[i]] <- NULL
    }
  }

  return(las_1)
}
