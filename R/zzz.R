.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "txshift v", utils::packageDescription("txshift")$Version,
    ": ", utils::packageDescription("txshift")$Title
  ))
}
