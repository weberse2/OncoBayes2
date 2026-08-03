# Full-sampling variant of start-example.R, swapped in at docs-build time so
# the pkgdown website renders examples with the package default sampling.
# These plain `#` comment lines are NOT emitted by roxygen2; only the `#'`
# lines below become part of the rendered @examples. We still define
# .user_mc_options (with no slim override) so the paired stop-example.R
# restore via options(.user_mc_options) stays valid.
#' @examples
#' .user_mc_options <- options()
#'
