is_ci <- function() nzchar(Sys.getenv("CI"))

is_cran <- function() {
  tolower(Sys.getenv("NOT_CRAN", unset = "")) != "true"
}

if (!exists("gold_runs")) {
  if (!is_cran() & identical(Sys.getenv("TEST_VDIFFR"), "true")) {
    # load gold runs, possibly from cache
    .mc_options <- default_sampling()
    gold_runs <- load_gold(TRUE)
    options(.mc_options)
  } else {
    # load gold runs, use cache if available, but do short
    # sampling if not
    if (!is_cran() | is_ci()) {
      .mc_options <- very_fast_sampling()
    } else {
      .mc_options <- fake_sampling()
    }
    gold_runs <- load_gold(TRUE)
    options(.mc_options)
  }
}
