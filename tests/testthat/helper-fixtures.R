test_fixture_path <- function(name) {
  testthat::test_path("fixtures", paste0(name, ".rds"))
}

load_test_fixture <- function(name) {
  path <- test_fixture_path(name)
  testthat::skip_if_not(
    file.exists(path),
    paste("Missing test fixture:", name)
  )
  readRDS(path)
}
