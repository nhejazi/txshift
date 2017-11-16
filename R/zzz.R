.onAttach <- function(...) {
  packageStartupMessage(paste("shifttx: Targeted Learning with",
                              "Stochastic Interventions"))
  packageStartupMessage("Version: ",
                        utils::packageDescription("shifttx")$Version)
}

