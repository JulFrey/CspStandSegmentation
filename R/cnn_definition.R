# ============================================================
# cnn_definition.R
# EfficientNetV2-S CircleCNN definition + checkpoint loader
#
# This file is intentionally "definition-only":
# - defines CircleCNN (nn_module)
# - defines circlecnn_load_checkpoint()
# - does NOT run training or load data
# ============================================================

if (!requireNamespace("torch", quietly = TRUE)) {
  stop("Package 'torch' is required. Install it first.")
}
if (!requireNamespace("torchvision", quietly = TRUE)) {
  stop("Package 'torchvision' is required. Install it first.")
}

# ------------------------------------------------------------
# CircleCNN (EfficientNetV2-S backbone)
# - Adapter 4 -> 3 channels for EfficientNet
# - Optional CoordConv grids added in forward() if input has C=2
# - Outputs: list(logit_stem=[B,1], reg=[B,3]) where reg=(cx_norm, cy_norm, r_norm)
# ------------------------------------------------------------
CircleCNN <- torch::nn_module(
  initialize = function(in_channels = 4,
                        dropout_p = 0.10,
                        use_pretrained = TRUE) {

    # 4-channel input: (density, zrel/verticality, xgrid, ygrid)
    self$adapter <- torch::nn_conv2d(in_channels, 3, kernel_size = 1, bias = FALSE)

    # Learnable gains for first 2 channels (density, zrel/verticality)
    # stored as log-gains => exp() makes them positive
    self$chan_loggain <- torch::nn_parameter(
      torch::torch_tensor(c(0, log(2)), dtype = torch::torch_float())
    )

    # EfficientNetV2-S backbone
    self$backbone <- torchvision::model_efficientnet_v2_s(pretrained = use_pretrained, progress = TRUE)

    self$drop <- torch::nn_dropout(dropout_p)

    feat_dim <- 1280
    self$head_stem <- torch::nn_linear(feat_dim, 1)
    self$head_reg  <- torch::nn_linear(feat_dim, 3)
  },

  forward = function(x) {
    x <- x$to(dtype = torch::torch_float())

    # Expect batched input: [B, C, H, W]
    B <- x$size(1)
    C <- x$size(2)
    H <- x$size(3)
    W <- x$size(4)

    # If input is [B,2,H,W], add CoordConv grids -> [B,4,H,W]
    if (C == 2) {
      dev <- x$device
      xlin <- torch::torch_linspace(-1, 1, steps = W, device = dev)$view(c(1, 1, 1, W))$expand(c(B, 1, H, W))
      ylin <- torch::torch_linspace(-1, 1, steps = H, device = dev)$view(c(1, 1, H, 1))$expand(c(B, 1, H, W))
      x <- torch::torch_cat(list(x, xlin, ylin), dim = 2)
      C <- 4
    }

    # Fail loudly if not 2 or 4
    if (C != 4) {
      stop(paste0("CircleCNN forward expected C=2 or C=4 channels, got C=", C))
    }

    # Apply learnable gains to channels 1:2 only
    gain <- torch::torch_exp(self$chan_loggain)$view(c(1, 2, 1, 1))  # [1,2,1,1]
    x_main <- x$narrow(2, 1, 2) * gain
    x_rest <- x$narrow(2, 3, 2)
    x <- torch::torch_cat(list(x_main, x_rest), dim = 2)

    # Adapter 4 -> 3
    x <- self$adapter(x)

    # EfficientNet features
    x <- self$backbone$features(x)
    x <- self$backbone$avgpool(x)

    # Flatten to [B, 1280]
    feat <- torch::torch_flatten(x, start_dim = 2)
    feat <- self$drop(feat)

    logit_stem <- self$head_stem(feat)
    reg_raw <- self$head_reg(feat)

    # cx, cy unconstrained; r positive
    cx <- reg_raw$select(2, 1)$unsqueeze(2)
    cy <- reg_raw$select(2, 2)$unsqueeze(2)
    r  <- torch::nnf_softplus(reg_raw$select(2, 3))$unsqueeze(2)
    reg <- torch::torch_cat(list(cx, cy, r), dim = 2)

    list(logit_stem = logit_stem, reg = reg)
  }
)

# ------------------------------------------------------------
# Helper: load a checkpoint (.pt) and return a ready-to-use model
# ------------------------------------------------------------
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
  if (!("adapter.weight" %in% names(ckpt$model_state))) {
    stop("Checkpoint model_state has no 'adapter.weight' key. Can't infer in_channels.")
  }

  # adapter.weight shape is [3, in_channels, 1, 1]
  in_ch <- as.integer(ckpt$model_state[["adapter.weight"]]$size()[2])

  model <- CircleCNN(in_channels = in_ch, use_pretrained = use_pretrained)
  model$load_state_dict(ckpt$model_state)
  model$to(device = device)
  model$eval()

  list(
    model = model,
    device = device,
    threshold = if (!is.null(ckpt$threshold)) ckpt$threshold else NULL,
    window_size_m = if (!is.null(ckpt$window_size_m)) ckpt$window_size_m else NULL
  )
}