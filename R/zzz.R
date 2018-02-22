.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "txshift v", utils::packageDescription("txshift")$Version,
    ": Targeted Learning with Stochastic Interventions"
  ))
}
