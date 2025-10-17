#' @section Group and strata definitions:
#'
#' The groups and strata as defined when running the \code{blrm_exnex}
#' analysis cannot be changed at a later stage. As a result no
#' evaluations can be performed for groups which have not been present
#' in the data set used for running the analysis. However, it is
#' admissible to code the group (and/or stratum) column as a
#' \code{factor} which contains empty levels. These groups are thus
#' not contained in the fitting data set and they are assigned by
#' default to the first stratum. In addition priors must be setup for
#' these groups (and/or strata). These empty group (and/or strata)
#' levels are then allowed in subsequent evaluations. This enables the
#' evaluation of the hierarchical model in terms of representing a
#' prior for future groups.
