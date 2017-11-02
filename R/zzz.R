.onAttach <- function(...) {
  packageStartupMessage(paste("shifttx: Targeted Learning with",
                              "Stochastic Treatment Regimes"))
  packageStartupMessage("Version: ",
                        utils::packageDescription("shifttx")$Version)
}

