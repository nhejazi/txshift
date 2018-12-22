#' NAME
#'
#' description
#'
#' @param fold ...
#' @param long_data ...
#' @param wts ...
#' @param lambda_seq ...
#'
#' @importFrom stats aggregate plogis
#' @importFrom origami training validation fold_index
#' @importFrom hal9001 fit_hal make_design_matrix apply_copy_map
#' @importFrom assertthat assert_that
#'
#' @export
#
cv_haldensify <- function(fold, long_data, wts = rep(1, nrow(long_data)),
                          lambda_seq = exp(seq(-1, -13, length = 1000))) {
  # make training and validation folds
  train_set <- origami::training(long_data)
  valid_set <- origami::validation(long_data)

  # subset observation-level weights to the correct size
  wts_train <- wts[fold$training_set]
  wts_valid <- wts[fold$validation_set]

  # fit a HAL regression on the training set
  hal_fit_train <- hal9001::fit_hal(X = as.matrix(train_set[, -c(1, 2)]),
                                    Y = as.numeric(train_set$in_bin),
                                    fit_type = "glmnet",
                                    use_min = TRUE,
                                    family = "binomial",
                                    return_lasso = TRUE,
                                    lambda = lambda_seq,
                                    fit_glmnet = TRUE,
                                    standardize = FALSE,   # pass to glmnet
                                    weights = wts_train,   # pass to glmnet
                                    yolo = FALSE)

  # get intercept and coefficient fits for this value of lambda from glmnet
  alpha_hat <- hal_fit_train$glmnet_lasso$a0
  betas_hat <- hal_fit_train$glmnet_lasso$beta
  coefs_hat <- rbind(alpha_hat, betas_hat)

  # make design matrix for validation set manually
  pred_x_basis <- hal9001:::make_design_matrix(as.matrix(valid_set),
                                               hal_fit_train$basis_list)
  pred_x_basis <- hal9001:::apply_copy_map(pred_x_basis,
                                           hal_fit_train$copy_map)
  pred_x_basis <- cbind(rep(1, nrow(valid_set)), pred_x_basis)

  # manually predict along sequence of lambdas
  preds_logit <- pred_x_basis %*% coefs_hat
  preds <- stats::plogis(as.matrix(preds_logit))

  # initialize lists to store predictions and loss across individuals
  hazards_pred_list <- vector("list", length(unique(valid_set$obs_id)))

  # compute hazard for a given observation by looping over individuals
  for (id in unique(valid_set$obs_id)) {
    # get predictions for the current observation only
    idx_this_obs <- which(unique(valid_set$obs_id) %in% id)
    preds_this_obs <- matrix(preds[which(valid_set$obs_id == id), ],
                             ncol = ncol(preds))
    n_records_this_obs <- nrow(preds_this_obs)

    # NOTE: pred_hazard = (1 - pred) if 0 in this bin * pred if 1 in this bin
    if (no_records_this_obs > 1) {
      hazard_prefailure <- (1 - preds_this_obs[-n_records_this_obs, ])
      hazard_failure <- preds_this_obs[n_records_this_obs, ]
      hazards_pred <- rbind(hazard_prefailure, hazard_failure)
    } else {
      hazards_pred <- preds_this_obs
    }

    # check dimensions
    assertthat::assert_that(all(dim(preds_this_obs) == dim(hazards_pred)))

    # multiply hazards across rows to construct the individual-level hazard
    hazards_pred <- matrix(apply(hazards_pred, 2, prod), nrow = 1)
    hazards_pred_list[[idx_this_obs]] <- hazards_pred
  }

  # aggregate predictions across observations
  hazards_pred <- do.call(rbind, hazards_pred_list)

  # collapse weights to the observation level
  wts_valid_reduced <- stats::aggregate(wts_valid, list(valid_set$obs_id),
                                        unique)
  colnames(wts_valid_reduced) <- c("id", "weight")

  # construct output
  out <- list(preds = hazards_pred,
              ids = wts_valid_reduced$id,
              wts = wts_valid_reduced$weight,
              fold = origami::fold_index())
  return(out)
}

#' NAME
#'
#' description
#'
#' @param A ...
#' @param W ...
#' @param wts ...
#' @param grid_type ...
#' @param n_bins ...
#' @param width ...
#' @param lambda_seq ...
#'
#' @importFrom origami make_folds cross_validate
#' @importFrom hal9001 fit_hal
#'
#' @export
#
haldensify <- function(A, W, wts = rep(1, length(A)),
                       grid_type = c("equal_range", "equal_mass",
                                     "equal_width"),
                       n_bins = 10, width = NULL,
                       lambda_seq = exp(seq(-1, -13, length = 1000))
                      ) {
  # re-format input data into long hazards structure
  long_data <- format_long_hazards(A = A, W = W, wts = wts, type = grid_type,
                                   n_bins = n_bins, width = width)

  # extract wieghts from long format data structure
  wts_long <- long_data$wts
  long_data[, wts := NULL]

  # make folds with origami
  folds <- origami::make_folds(long_data, cluster_ids = long_data$obs_id)

  # call cross_validate on cv_density function...
  haldensity <- origami::cross_validate(cv_fun = cv_haldensify,
                                        folds = folds,
                                        long_data = long_data,
                                        wts = wts_long,
                                        lambda_seq = lambda_seq,
                                        use_future = FALSE,
                                        .combine = FALSE
                                       )

  # re-organize output from origami::cross_validate
  hazards_pred <- do.call(rbind, haldensity$preds)
  obs_ids <- do.call(c, haldensity$ids)
  obs_wts <- do.call(c, haldensity$wts)

  # compute loss for the given individual
  hazards_loss <- apply(hazards_pred, 2, function(x) {
                          pred_weighted <- x * obs_wts
                          loss_weighted <- -log(pred_weighted)
                          return(loss_weighted)
                       })

  # take column means to have average loss across sequence of lambdas
  loss_mean <- colMeans(hazards_loss)
  lambda_loss_min <- lambda_seq[which.min(loss_mean)]

  # fit a HAL regression on the full data set with the CV-selected lambda
  # NOTE: How should the data structure for this regression be formatted?
  hal_fit <- hal9001::fit_hal(X = as.matrix(long_data[, -c(1, 2)]),
                              Y = as.numeric(long_data$in_bin),
                              fit_type = "glmnet",
                              use_min = TRUE,
                              family = "binomial",
                              return_lasso = TRUE,
                              lambda = lambda_loss_min,
                              fit_glmnet = TRUE,
                              standardize = FALSE,   # pass to glmnet
                              weights = wts_long,    # pass to glmnet
                              yolo = FALSE)

  # predict conditional density estimate from HAL fit on full data
  density_pred <- predict(hal_fit, new_data = long_data)

  # TODO: construct output
}

#' Generate long format hazards data for conditional density estimation
#'
#' description forthcoming
#'
#' @param A ...
#' @param W ...
#' @param wts ...
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
format_long_hazards <- function(A, W, wts = rep(1, length(A)),
                                type = c("equal_range", "equal_mass",
                                         "equal_width"),
                                n_bins = 10, width = NULL) {
  # clean up arguments
  type <- match.arg(type)

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

  # create list to capture data tables for each observation
  df_all_obs <- vector("list", length(A))

  # loop over observations to create expanded set of records for each
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

    # get correct value of weights and repeat along intervals
    # NOTE: the weights are always a vector
    obs_wts <- rep(wts[i], nrow(all_bins))

    # create data table with membership indicator and interval limits
    suppressWarnings(
      hazards_df <- data.table::as.data.table(cbind(id, bin_indicator,
                                                    all_bins, obs_w,
                                                    obs_wts))
    )

    # trim records to simply end at the failure time for a given observation
    hazards_df_reduced <- hazards_df[seq_len(bin_id[i]), ]

    # give explicit names and add to appropriate position in list
    data.table::setnames(hazards_df_reduced,
                         c("obs_id", "in_bin", "bin_id", names_w, "wts"))
    df_all_obs[[i]] <- hazards_df_reduced
  }

  # combine observation-level hazards data into larger structure
  out <- data.table::rbindlist(df_all_obs)
  return(out)
}

