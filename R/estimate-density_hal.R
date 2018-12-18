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
cv_haldensify <- function() {
  ...
}

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
haldensify <- function() {
  ...
}

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
format_long_hazards <- function(A, W, bin_width) {
  bounds <- range(x)
  self$grids <- seq(from = bounds[1], to = bounds[2], by = bin_width)

  all_df <- list()
  b <- 1
  for (i in self$grids) {
    Y <- ( (i - .5 * self$bin_width <= x) & (x < i + .5 * self$bin_width)) + 0L
    all_df[[b]] <- data.frame(id = 1:length(x), Y = Y, box = i)
    b <- b + 1
  }
  all_df <- do.call("rbind", all_df)
  all_df$Y <- as.numeric(as.character(all_df$Y)) # turn factor to numeric
  # self$longiDF <- all_df[order( all_df[,1], all_df[,3] ),2:3]
  return(all_df[order(all_df[, 1], all_df[, 3]), 2:3])

  # NOTE: ONLY FOR USE WITH CONDENSIER, TO PASS IDs TO CV.GLMNET PROPERLY
  #repeated_obs_idx <- c(0, which(task$Y == 1))
  #ids_times_arg <- repeated_obs_idx - dplyr::lag(repeated_obs_idx)
  #ids_times_arg <- ids_times_arg[!is.na(ids_times_arg)]
  #ids_in <- seq_len(length(which(task$Y == 1)))
  #ids <- rep(ids_in, times = ids_times_arg)
  #foldid <- origami:::folds2foldvec(make_folds(cluster_ids = ids))
}

