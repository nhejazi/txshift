#' Shift Treatments


tx_shift <- function(a, w = NULL, delta, type = "additive",
                     direc = c("down", "up")) {
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

