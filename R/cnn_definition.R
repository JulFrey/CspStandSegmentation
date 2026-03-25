# ============================================================
# cnn_definition.R
# CircleCNN definitions + checkpoint loader
#
# This file is intentionally "definition-only":
# - defines checkpoint-compatible CircleCNN architectures
# - defines circlecnn_load_checkpoint()
# - does NOT run training or load data
#
# Supported checkpoints:
# - EfficientNetV2-S based CircleCNN
# - MobileNetV3-Large based CircleCNN
#
# Loader logic:
# - inspects ckpt$model_state
# - reconstructs the matching architecture automatically
# ============================================================

if (!requireNamespace("torch", quietly = TRUE)) {
  stop("Package 'torch' is required. Install it first.")
}
if (!requireNamespace("torchvision", quietly = TRUE)) {
  stop("Package 'torchvision' is required. Install it first.")
}

# ------------------------------------------------------------
# Helper: append CoordConv grids if input arrives as [B,2,H,W]
# Input contract used by inventory inference:
# - C=2: (density, zrel/verticality)
# - C=4: (density, zrel/verticality, xgrid, ygrid)
# ------------------------------------------------------------
.add_coordconv_if_needed <- function(x) {
  x <- x$to(dtype = torch::torch_float())

  sx <- as.integer(x$size())
  if (length(sx) != 4) {
    stop("CircleCNN forward expected a 4D tensor [B,C,H,W].")
  }

  B <- sx[1]
  C <- sx[2]
  H <- sx[3]
  W <- sx[4]

  if (C == 2L) {
    dev <- x$device
    xlin <- torch::torch_linspace(-1, 1, steps = W, device = dev)$
      view(c(1, 1, 1, W))$
      expand(c(B, 1, H, W))
    ylin <- torch::torch_linspace(-1, 1, steps = H, device = dev)$
      view(c(1, 1, H, 1))$
      expand(c(B, 1, H, W))
    x <- torch::torch_cat(list(x, xlin, ylin), dim = 2)
    C <- 4L
  }

  if (C != 4L) {
    stop(paste0("CircleCNN forward expected C=2 or C=4 channels, got C=", C))
  }

  x
}

# ------------------------------------------------------------
# EfficientNetV2-S CircleCNN
# - 4 -> 3 adapter conv
# - learnable gains on first 2 channels (chan_loggain)
# - EfficientNetV2-S backbone
# ------------------------------------------------------------
CircleCNN_EfficientNetV2S <- torch::nn_module(
  initialize = function(in_channels = 4,
                        dropout_p = 0.10,
                        use_pretrained = TRUE) {

    self$adapter <- torch::nn_conv2d(in_channels, 3, kernel_size = 1, bias = FALSE)

    # Learnable gains for first 2 channels (density, zrel/verticality)
    self$chan_loggain <- torch::nn_parameter(
      torch::torch_tensor(c(0, log(2)), dtype = torch::torch_float())
    )

    self$backbone <- torchvision::model_efficientnet_v2_s(
      pretrained = use_pretrained,
      progress = TRUE
    )

    self$drop <- torch::nn_dropout(dropout_p)

    feat_dim <- 1280
    self$head_stem <- torch::nn_linear(feat_dim, 1)
    self$head_reg  <- torch::nn_linear(feat_dim, 3)
  },

  forward = function(x) {
    x <- .add_coordconv_if_needed(x)

    # Apply learnable gains only to channels 1:2
    gain <- torch::torch_exp(self$chan_loggain)$view(c(1, 2, 1, 1))
    x_main <- x$narrow(2, 1, 2) * gain
    x_rest <- x$narrow(2, 3, 2)
    x <- torch::torch_cat(list(x_main, x_rest), dim = 2)

    # Adapter 4 -> 3
    x <- self$adapter(x)

    # EfficientNet features
    feat_map <- self$backbone$features(x)
    feat_map <- self$backbone$avgpool(feat_map)

    # Flatten [B,C,1,1] -> [B,C]
    feat <- feat_map$view(c(feat_map$size(1), -1))
    feat <- self$drop(feat)

    logit_stem <- self$head_stem(feat)
    reg_raw <- self$head_reg(feat)

    cx <- reg_raw$select(2, 1)$unsqueeze(2)
    cy <- reg_raw$select(2, 2)$unsqueeze(2)
    r  <- torch::nnf_softplus(reg_raw$select(2, 3))$unsqueeze(2)
    reg <- torch::torch_cat(list(cx, cy, r), dim = 2)

    list(logit_stem = logit_stem, reg = reg)
  }
)

# ------------------------------------------------------------
# MobileNetV3-Large CircleCNN
# - 4 -> 3 adapter conv
# - no chan_loggain parameter
# - MobileNetV3-Large backbone
# ------------------------------------------------------------
CircleCNN_MobileNetV3Large <- torch::nn_module(
  initialize = function(in_channels = 4,
                        dropout_p = 0.10,
                        use_pretrained = TRUE) {

    self$adapter <- torch::nn_conv2d(in_channels, 3, kernel_size = 1, bias = FALSE)

    self$backbone <- torchvision::model_mobilenet_v3_large(
      pretrained = use_pretrained,
      progress = TRUE
    )

    self$drop <- torch::nn_dropout(dropout_p)

    # MobileNetV3-Large feature dimension after avgpool/flatten
    feat_dim <- 960
    self$head_stem <- torch::nn_linear(feat_dim, 1)
    self$head_reg  <- torch::nn_linear(feat_dim, 3)
  },

  forward = function(x) {
    x <- .add_coordconv_if_needed(x)

    # Adapter 4 -> 3
    x <- self$adapter(x)

    # MobileNet features
    feat_map <- self$backbone$features(x)

    # Keep compatibility across torchvision versions
    if (!is.null(self$backbone$avgpool)) {
      feat_map <- self$backbone$avgpool(feat_map)
    } else {
      feat_map <- torch::nnf_adaptive_avg_pool2d(feat_map, output_size = c(1, 1))
    }

    # Flatten [B,C,1,1] -> [B,C]
    feat <- feat_map$view(c(feat_map$size(1), -1))
    feat <- self$drop(feat)

    logit_stem <- self$head_stem(feat)
    reg_raw <- self$head_reg(feat)

    cx <- reg_raw$select(2, 1)$unsqueeze(2)
    cy <- reg_raw$select(2, 2)$unsqueeze(2)
    r  <- torch::nnf_softplus(reg_raw$select(2, 3))$unsqueeze(2)
    reg <- torch::torch_cat(list(cx, cy, r), dim = 2)

    list(logit_stem = logit_stem, reg = reg)
  }
)

# Backward-compatible alias. Not used by the loader internally,
# but kept so external code can still instantiate CircleCNN().
CircleCNN <- CircleCNN_MobileNetV3Large

# ------------------------------------------------------------
# Helper: infer which architecture a checkpoint needs
# ------------------------------------------------------------
.infer_circlecnn_architecture <- function(model_state) {
  state_names <- names(model_state)

  if (is.null(state_names) || length(state_names) == 0L) {
    stop("Checkpoint model_state has no parameter names.")
  }

  if (!("adapter.weight" %in% state_names)) {
    stop("Checkpoint model_state has no 'adapter.weight' key. Can't infer in_channels.")
  }
  if (!("head_stem.weight" %in% state_names)) {
    stop("Checkpoint model_state has no 'head_stem.weight' key. Can't infer architecture.")
  }

  has_chan_loggain <- "chan_loggain" %in% state_names
  feat_dim <- as.integer(model_state[["head_stem.weight"]]$size()[2])

  if (has_chan_loggain || feat_dim == 1280L) {
    return("efficientnet_v2_s")
  }
  if (feat_dim == 960L) {
    return("mobilenet_v3_large")
  }

  stop(
    paste0(
      "Unsupported CircleCNN checkpoint. Could not infer architecture from state_dict ",
      "(has_chan_loggain=", has_chan_loggain,
      ", head_stem.in_features=", feat_dim, ")."
    )
  )
}

#' Load a CircleCNN checkpoint
#'
#' Reconstructs a checkpoint-compatible CircleCNN architecture from a saved
#' `.pt` file, loads the trained weights, moves the model to the selected
#' device, and returns a ready-to-use model object.
#'
#' The loader currently supports EfficientNetV2-S and MobileNetV3-Large based
#' CircleCNN checkpoints and infers the correct architecture from the checkpoint
#' state dictionary.
#'
#' @param ckpt_path Character path to a saved checkpoint file.
#' @param device Optional torch device. If `NULL`, CUDA is used when available,
#'   otherwise CPU.
#' @param use_pretrained Logical; whether to initialise the backbone with
#'   pretrained torchvision weights before loading the checkpoint state.
#'
#' @return A list with components:
#' \describe{
#'   \item{model}{A torch model in evaluation mode.}
#'   \item{device}{The torch device used.}
#'   \item{architecture}{The inferred model architecture.}
#'   \item{threshold}{Checkpoint threshold if present, otherwise `NULL`.}
#'   \item{window_size_m}{Checkpoint window size if present, otherwise `NULL`.}
#' }
#'
#' @examples
#' \dontrun{
#'   obj <- circlecnn_load_checkpoint("CNN_MobileNetV3Large_v1.pt")
#'   model <- obj$model
#' }

circlecnn_load_checkpoint <- function(ckpt_path,
                                      device = NULL,
                                      use_pretrained = FALSE) {

  if (!file.exists(ckpt_path)) {
    stop(paste0("Checkpoint not found: ", ckpt_path))
  }

  if (is.null(device)) {
    device <- torch::torch_device(if (torch::cuda_is_available()) "cuda" else "cpu")
  }

  ckpt <- torch::torch_load(ckpt_path, device = "cpu")

  if (is.null(ckpt$model_state)) {
    stop("Checkpoint does not contain $model_state.")
  }

  model_state <- ckpt$model_state
  arch <- .infer_circlecnn_architecture(model_state)

  # adapter.weight shape is [3, in_channels, 1, 1]
  in_ch <- as.integer(model_state[["adapter.weight"]]$size()[2])

  model <- switch(
    arch,
    efficientnet_v2_s = CircleCNN_EfficientNetV2S(
      in_channels = in_ch,
      use_pretrained = use_pretrained
    ),
    mobilenet_v3_large = CircleCNN_MobileNetV3Large(
      in_channels = in_ch,
      use_pretrained = use_pretrained
    ),
    stop(paste0("Unsupported CircleCNN architecture: ", arch))
  )

  model$load_state_dict(model_state)
  model$to(device = device)
  model$eval()

  list(
    model = model,
    device = device,
    architecture = arch,
    threshold = if (!is.null(ckpt$threshold)) ckpt$threshold else NULL,
    window_size_m = if (!is.null(ckpt$window_size_m)) ckpt$window_size_m else NULL
  )
}
