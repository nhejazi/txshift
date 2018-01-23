.onAttach <- function(...) {
  packageStartupMessage(paste(
    "txshift: Targeted Learning with",
    "Stochastic Interventions"
  ))
  packageStartupMessage(
    "Version: ",
    utils::packageDescription("txshift")$Version
  )
}
