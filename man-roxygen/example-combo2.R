#' @examples
#' dref <- c(6, 960)
#'
#' num_comp   <- 2 # two investigational drugs
#' num_inter  <- 1 # one drug-drug interaction needs to be modeled
#' num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
#' num_strata <- 1 # no stratification needed
#'
#' blrmfit <- blrm_exnex(
#'   cbind(num_toxicities, num_patients - num_toxicities) ~
#'     1 + I(log(drug_A / dref[1])) |
#'       1 + I(log(drug_B / dref[2])) |
#'       0 + I(drug_A / dref[1] * drug_B / dref[2]) |
#'       group_id,
#'   data = codata_combo2,
#'   prior_EX_mu_comp  = list(mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
#'                            mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))),
#'   prior_EX_tau_comp = list(mixmvnorm(c(1,
#'                                        log(0.250), log(0.125),
#'                                        diag(c(log(4)/1.96, log(4)/1.96)^2))),
#'                            mixmvnorm(c(1,
#'                                        log(0.250), log(0.125),
#'                                        diag(c(log(4)/1.96, log(4)/1.96)^2)))),
#'   prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
#'   prior_EX_tau_inter = mixmvnorm(c(1, log(0.125), (log(4) / 1.96)^2)),
#'   prior_is_EXNEX_comp = rep(FALSE, num_comp),
#'   prior_is_EXNEX_inter = rep(FALSE, num_inter),
#'   prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
#'   prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
#'   prior_tau_dist = 1
#' )
