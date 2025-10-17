#' @param x Character giving the parameter name to be mapped to the
#'     x-axis.  This also supports 'tidy' parameter selection by
#'     specifying \code{x = vars(...)}, where \code{...} is specified
#'     the same way as in \code{\link[dplyr:select]{dplyr::select()}}
#'     and similar functions. Examples of using \code{x} in this way
#'     can be found in the examples. For \code{blrm_trial} methods, it
#'     defaults to the first entry in \code{summary(blrm_trial,
#'     "drug_info")$drug_name}.
#' @param group Grouping variable(s) whose levels will be mapped to
#'     different facets of the plot. \code{group} can be a character
#'     vector, tidy parameter(s) of the form \code{group = vars(...)},
#'     or a formula to be passed directly to
#'     \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}. For
#'     \code{blrm_trial} methods, it defaults to \code{group_id}, plus
#'     all entries of \code{summary(blrm_trial,
#'     "drug_info")$drug_name} except the first, which is mapped to
#'     \code{x}.
