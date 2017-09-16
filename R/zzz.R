.onAttach <- function(...) {
  packageStartupMessage("shifttx: Estimate Causal Effects with Stochastic Treatments")
  packageStartupMessage("Version: ",
                        utils::packageDescription("shifttx")$Version)
}

