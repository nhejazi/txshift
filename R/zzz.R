.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "txshift v", utils::packageDescription("txshift")$Version,
    ": Targeted Learning of the Causal Effects of Stochastic Interventions"
  ))
}
