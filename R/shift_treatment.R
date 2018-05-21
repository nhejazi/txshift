#' Treatment Shifting
#'
#'
#' @param A A \code{numeric} vector of observed treatment values.
#' @param W A \code{numeric} matrix of observed baseline covariate values.
#' @param delta A \code{numeric} indicating the magnitude of the shift to be
#'  computed for the treatment \code{A}.
#' @param type A \code{character} vector describing the type of shift to be
#'  performed on the treatment. Currently, this is limited only to "additive",
#'  though the relevant theory supports other forms of shifting as well.
#' @param direc A \code{character} vector describing the direction in which to
#'  shift the treatment --- currently, only "up" and "down" are supported.
#'
#' @keywords internal
#
tx_shift <- function(A,
                     W = NULL,
                     delta,
                     type = c("additive"),
                     direc = c("down", "up")) {
  # could support more than just additive shifts?
  if (type == "additive") {
    if (direc == "down") {
      a_shift <- A - delta
    }
    if (direc == "up") {
      a_shift <- A + delta
    }
  }
  # can support other types of treatment shifts
  # (e.g., multiplicative shifting?)
  return(a_shift)
}
