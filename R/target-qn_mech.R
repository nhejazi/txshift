utils::globalVariables(c("key"))

#' Targeting for Part of EIF Corresponding to Weighted Censoring Distribution
#'
#' DESCRIPTION HERE
#'
#' @param Qn_shift ...
#' @param ipc_weights ...
#' @param data_in ...
#'
#' @importFrom data.table set
#' @importFrom stats optim
#'
#' @keywords internal
#'
#' @export
#'
#' @author Nima Hejazi
#' @author David Benkeser
#
target_qn <- function(Qn_shift, ipc_weights, data_in) {
    # figure out what rows are duplicated by removing Y from
    # data in and finding duplicated rows
    n_samp <- nrow(data_in)
    rm_col <- "Y"

    # NOTE: it seems like we will want to move this out of this inner function
    # as these duplicates are not going to change during the procedure.
    # Thus, we should be able to figure out who's duplicated just once and then 
    # keep updating values based on that.
    dups <- duplicated(data_in[ , !rm_col, by = key(data_in), with = FALSE])
    rev_dups <- duplicated(data_in[ , !rm_col, by = key(data_in),
                           with = FALSE], fromLast = TRUE)
    data.table::set(data_in, j = "any_dup", value = (dups | rev_dups))
    data.table::set(data_in, j = "dups", value = dups)
    data.table::set(data_in, j = "id", value = seq_len(n_samp))

    # if there are no duplicates then qn = 1/n * ipc_weight
    if (all(data_in$dups)) {
      qn_dens <- (1 / n_samp) * ipc_weights
    } else {
      # initialize an empty vector of all NaN's
      qn_dens <- rep(NaN, n_samp)
      # anyone who has a unique value gets 1/n * ipc_weight
      qn_dens[!data_in$any_dup] <- (1 / n_samp) *
          ipc_weights[ipc_weights > 0][!data_in$any_dup]
      # fill in values for duplicated folks
      # here the loop is over the indices of the folks who have
      # the first of several duplicated values
      for (i in which(data_in$any_dup & !data_in$dup)) {
        # indicator to be set to TRUE for anyone who has this particular (a,w)
        data_in$tmp_idx <- FALSE
        # find observations with this value of (A, W) by looping over folks that
        # have any duplicated value
        for (j in which(data_in$any_dup)) {
          if (identical(data_in[j, -c("Y", "dups", "id", "any_dup", "tmp_idx"),
                       with = FALSE], data_in[i,-c("Y", "dups", "id", "any_dup",
                                                   "tmp_idx"), with = FALSE])) {
            data_in$tmp_idx[j] <- TRUE
          }
        }
        # now we can compute the density for these folks as 1/n
        # times the sum of their ipc_weights
        qn_dens[data_in$tmp_idx] <- 1/n_samp *
            sum(ipc_weights[ipc_weights > 0][data_in$tmp_idx])
      }
    }

    # if no duplicates duplicated_idx == integer(0)
    if (all(data_in$dups)) {
      psi <- sum(Qn_shift * qn_dens)
    } else {
      psi <- sum(Qn_shift[!data_in$dup] * qn_dens[!data_in$dup])
    }

    # compute D^F_W part of the EIF
    Dfw_eif <- Qn_shift - psi

    # bookkeeping with zeros for both D^F_W and qn terms
    qn_with_zeros <- rep(0, length(ipc_weights))
    Dfw_with_zeros <- rep(0, length(ipc_weights))
    qn_with_zeros[ipc_weights != 0] <- qn_dens
    Dfw_with_zeros[ipc_weights != 0] <- Dfw_eif

    # solve for root of the score equation of epsilon (one option using optim)
    ipcw_nloglik <- function(eps_in, Dfw, qn, ipc_weights) {
      return(sum(- ipc_weights * log((1 + eps_in * Dfw) * qn)))
    }

    # optimize over the negative log-likelihood loss
    fluc_mod <- stats::optim(par = 0, fn = ipcw_nloglik, Dfw = Dfw_eif,
                             qn = qn_dens,
                             ipc_weights = ipc_weights[ipc_weights != 0],
                             method = "Brent",
                             lower = max(-1 / Dfw_eif[Dfw_eif > 0]),
                             upper = min(-1 / Dfw_eif[Dfw_eif < 0]))

    # # could try using Rsolnp -- doesn't work, no solution
    # constraint function
    # sum_to_1 <- function(eps_in, qn, Dfw, ipc_weights){
    #   sum_to_one <- sum(qn) - 1
    #   sum_to_one
    # }
    # # objective function
    # ipcw_nloglik <- function(eps_in, Dfw, qn, ipc_weights){
    #   return(sum(- ipc_weights * log((1 + eps_in * Dfw)*qn)))
    # }
    # fluc_mod <- solnp(pars = 0, fun = ipcw_nloglik, 
    #                   UB = min(-1/Dfw_eif[Dfw_eif < 0]),
    #                   LB = max(-1/Dfw_eif[Dfw_eif > 0]),
    #                   eqfun = sum_to_1, eqB = c(0), 
    #                   Dfw = Dfw_eif, qn = qn_dens, 
    #                   ipc_weights = ipc_weights[ipc_weights != 0])

    # update the estimated density qn based on epsilon
    qn_update <- (1 + fluc_mod$par * Dfw_eif) * qn_dens
    return(list(eps = fluc_mod$par, qn = qn_update))
}

