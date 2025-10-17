#' @examples
#' ## Example from Neuenschwander, B., et al. (2009). Stats in Medicine
#'
#' dref <- 50
#'
#' ## Since there is no prior information the hierarchical model
#' ## is not used in this example by setting tau to (almost) 0.
#' blrmfit <- blrm_exnex(
#'   cbind(num_toxicities, num_patients - num_toxicities) ~
#'       1 + log(drug_A / dref) |
#'       0 |
#'       group_id,
#'   data = hist_SA,
#'   prior_EX_mu_comp = mixmvnorm(c(1, logit(1 / 2), log(1), diag(c(2^2, 1)))),
#'   ## Setting prior_tau_dist=NULL disables the hierarchical prior which is
#'   ## not required in this example as we analyze a single trial.
#'   prior_tau_dist = NULL,
#'   prior_PD = FALSE
#' )
