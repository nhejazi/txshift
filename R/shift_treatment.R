#' Shift Treatments
#'
#' description
#'
#' @param A ...
#' @param W ...
#' @param delta ...
#' @param type ...
#' @param direc ...
#'
#' @keywords internal
#
tx_shift <- function(A,
                     W = NULL,
                     delta,
                     type = "additive",
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

