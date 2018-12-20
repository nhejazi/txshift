#' NAME
#'
#' description
#'
#' @param fold ...
#' @param long_data ...
#' @param ... ...
#'
#' @importFrom origami training validation
#'
#' @export
#
cv_haldensify <- function(fold, long_data, lambda_seq) {
  # make training and validation folds
  train_set <- origami::training(long_data)
  valid_set <- origami::validation(long_data)

  # do stuff with HAL
  hal_lambda_seq <- hal9001::fit_hal(X = as.matrix(train_set[, -1]),
                                     Y = as.matrix(train_set[, 1]),
                                     #weights = ...,
                                     family = "binomial",
                                     fit_type = "glmnet",
                                     lambda = lambda_seq,
                                     use_min = TRUE,
                                     return_lasso = TRUE,
                                     yolo = FALSE)

  # get predictions on validation set to evaluate loss properly
  preds_hal_long <- predict(hal_lambda_seq,
                            new_data = as.matrix(valid_set[, -1]))


}

#' NAME
#'
#' description
#'
#' @param A ...
#' @param W ...
#' @param ... ...
#'
#' @importFrom origami make_folds cross_validate
#'
#' @export
#
haldensify <- function(A, W, ...) {
  # re-format input data into long hazards structure
  long_data <- format_long_hazards(A = A, W = W, n_bins = 10)

  # extract observation-level IDs and drop column
  obs_ids <- long_data$obs_id
  long_data[, obs_id := NULL]

  # make folds with origami
  folds <- origami::make_folds(long_data, cluster_ids = obs_ids)

  # call cross_validate on cv_density function...
  cv_haldensity <- origami::cross_validate(cv_fun = cv_haldensify,
                                           folds = folds,
                                           long_data = long_data,
                                           ...)
}

#' Generate long format hazards data for conditional density estimation
#'
#' description forthcoming
#'
#' @param A ...
#' @param W ...
#' @param type ...
#' @param n_bins ...
#' @param width ...
#'
#' @importFrom data.table as.data.table rbindlist setnames
#' @importFrom ggplot2 cut_interval cut_number cut_width
#' @importFrom assertthat assert_that
#'
#' @export
#
format_long_hazards <- function(A, W,
                                type = c("equal_range", "equal_mass",
                                         "equal_width"),
                                n_bins = 10, width = NULL) {
  # clean up arguments
  type <- match.arg(type)

  # get lower and upper limits of A
  bounds <- range(A)

  # set grid along A and find interval membership of observations along grid
  if (type == "equal_range") {
    bins <- ggplot2::cut_interval(A, n_bins, right = FALSE)
    bin_id <- as.numeric(bins)
  } else if (type == "equal_mass") {
    bins <- ggplot2::cut_number(A, n_bins, right = FALSE)
    bin_id <- as.numeric(bins)
  } else {
    assertthat::assert_that(!is.null(width))
    bins <- ggplot2::cut_width(A, width, closed = "left")
    bin_id <- as.numeric(bins)
  }

  df_all_obs <- vector("list", length(A))
  for (i in seq_along(A)) {
    # create repeating bin IDs for this subject (these map to intervals)
    all_bins <- matrix(seq_along(levels(bins)), ncol = 1)

    # create indicator and "turn on" indicator for interval membership
    bin_indicator <- rep(0, nrow(all_bins))
    bin_indicator[bin_id[i]] <- 1
    id <- rep(i, nrow(all_bins))

    # get correct value of baseline variables and repeat along intervals
    if (is.null(dim(W))) {
      # assume vector
      obs_w <- rep(W[i], nrow(all_bins))
      names_w <- "W"
    } else {
      # assume two-dimensional array
      obs_w <- rep(as.numeric(W[i, ]), nrow(all_bins))
      obs_w <- matrix(obs_w, ncol = ncol(W), byrow = TRUE)

      # use names from array if present
      if (is.null(names(W))) {
        names_w <- paste("W", seq_len(ncol(W)), sep = "_")
      } else {
        names_w <- names(W)
      }
    }

    # create data table with membership indicator and interval limits
    suppressWarnings(
      hazards_df <- data.table::as.data.table(cbind(id, bin_indicator,
                                                    all_bins, obs_w))
    )

    # give explicit names and add to appropriate position in list
    data.table::setnames(hazards_df, c("obs_id", "in_bin", "bin_id", names_w))
    df_all_obs[[i]] <- hazards_df
  }

  # combine observation-level hazards data into larger structure
  out <- data.table::rbindlist(df_all_obs)
  return(out)
}

