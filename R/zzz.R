.onAttach <- function(...) {
  packageStartupMessage(paste0(
    "this awesome package called txshift v", utils::packageDescription("txshift")$Version,
    ": ", utils::packageDescription("txshift")$Title
  ))
}
