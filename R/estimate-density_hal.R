#' NAME
#'
#' description
#'
#' @param args...
#'
#' @importFrom origami training validation
#'
#' @export
#
cv_haldensify <- function(fold, long_hazard_data, ...) {
  # make training and validation folds
  train_set <- origami::training(long_hazard_data)
  valid_set <- origami::validation(long_hazard_data)

  # do stuff with HAL
}

#' NAME
#'
#' description
#'
#' @param args...
#'
#' @importFrom origami make_folds cross_validate
#'
#' @export
#
haldensify <- function(long_hazard_data) {
  # extract observation-level IDs and drop column
  subj_ids <- long_hazard_data$id
  long_hazard_data[, id := NULL]

  # make folds with origami
  folds <- origami::make_folds(long_hazard_data, cluster_ids = subj_ids)

  # call cross_validate on cv_density function...
}

#' Generate long format hazards data for conditional density estimation
#'
#' description forthcoming
#'
#' @param A ...
#' @param W ...
#' @param n_bins ...
#'
#' @importFrom data.table as.data.table rbindlist setnames
#' @importFrom assertthat assert_that
#'
#' @export
#
format_long_hazards <- function(A, W, n_bins = 10) {
  # get lower and upper limits of A
  bounds <- range(A)

  # set grid along A and find interval membership of observations along grid
  grids <- seq(from = bounds[1], to = bounds[2], length.out = n_bins)
  bin_id <- findInterval(A, grids)

  df_all_obs <- vector("list", length(A))
  for (i in seq_along(A)) {
    # create matrix of upper and lower interval limits
    bins <- matrix(grids[c(1:2, 2:3, 3:4, 4:5, 5:6, 6:7, 7:8, 8:9, 9:n_bins)],
                   ncol = 2, byrow = TRUE)

    # create indicator and "turn on" indicator for interval membership
    bin_indicator <- rep(0, nrow(bins))
    bin_indicator[bin_id[i]] <- 1
    id <- rep(i, nrow(bins))

    # get correct value of baseline variables and repeat along intervals
    if (is.null(dim(W))) {
      # assume vector
      obs_w <- rep(W[i], nrow(bins))
      names_w <- "W"
    } else {
      # assume two-dimensional array
      obs_w <- rep(as.numeric(W[i, ]), nrow(bins))
      obs_w <- matrix(obs_w, ncol = ncol(W), byrow = TRUE)

      # use names from array if present
      if (is.null(names(W))) {
        names_w <- paste("W", 1:ncol(W), sep = "_")
      } else {
        names_w <- names(W)
      }
    }

    # create data table with membership indicator and interval limits
    suppressWarnings(
      hazards_df <- data.table::as.data.table(cbind(id, bin_indicator,
                                                    bins, obs_w))
    )

    # give explicit names and add to appropriate position in list
    data.table::setnames(hazards_df, c("id", "haz_ind", "lower", "upper",
                                       names_w))
    df_all_obs[[i]] <- hazards_df
  }

  # combine observation-level hazards data into larger structure
  out <- data.table::rbindlist(df_all_obs)

  # check that output object has the correct dimensionality and return
  assertthat::assert_that(nrow(out) == length(A) * (length(grids) - 1))
  return(out)
}

