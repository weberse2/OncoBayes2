#' Extract Diagnostic Quantities of \pkg{OncoBayes2} Models
#'
#' Extract quantities that can be used to diagnose sampling behavior
#' of the algorithms applied by \pkg{Stan} at the back-end of
#' \pkg{OncoBayes2}.
#'
#' @name diagnostic-quantities
#' @aliases log_posterior nuts_params rhat neff_ratio
#'
#' @param object A `blrmfit` or `blrmtrial` object.
#' @param pars An optional character vector of parameter names.
#'   For `nuts_params` these will be NUTS sampler parameter
#'   names rather than model parameters. If `pars` is omitted
#'   all parameters are included.
#' @param ... Arguments passed to individual methods.
#'
#' @return The exact form of the output depends on the method.
#'
#' @details For more details see
#'   [bayesplot::bayesplot-extractors()].
#'
#' @template start-example
#' @examples
#' example_model("single_agent", silent = TRUE)
#'
#' head(log_posterior(blrmfit))
#'
#' np <- nuts_params(blrmfit)
#' str(np)
#' # extract the number of divergence transitions
#' sum(subset(np, Parameter == "divergent__")$Value)
#'
#' head(rhat(blrmfit))
#' head(neff_ratio(blrmfit))
#'
#' @template stop-example
#'
NULL

#' @rdname diagnostic-quantities
#' @method log_posterior blrmfit
#' @importFrom bayesplot log_posterior
#' @export log_posterior
#' @export
log_posterior.blrmfit <- function(object, ...) {
  .contains_draws(object)
  lp <- posterior::as_draws_df(as_draws_array(object, variable = "lp__"))
  data.frame(
    Chain = lp$.chain,
    Iteration = lp$.iteration,
    Value = lp$lp__
  )
}

#' @rdname diagnostic-quantities
#' @method nuts_params blrmfit
#' @importFrom bayesplot nuts_params
#' @export nuts_params
#' @export
nuts_params.blrmfit <- function(object, pars = NULL, ...) {
  .contains_draws(object)
  .draws_diag_to_nuts_params(object$draws_diag, pars = pars)
}

#' @rdname diagnostic-quantities
#' @method rhat blrmfit
#' @importFrom bayesplot rhat
#' @export rhat
#' @export
rhat.blrmfit <- function(object, pars = NULL, ...) {
  .contains_draws(object)
  draws <- as_draws_array(object, variable = pars)
  out <- posterior::summarise_draws(draws, rhat = posterior::rhat)
  setNames(out$rhat, out$variable)
}

#' @rdname diagnostic-quantities
#' @method neff_ratio blrmfit
#' @importFrom bayesplot neff_ratio
#' @export neff_ratio
#' @export
neff_ratio.blrmfit <- function(object, pars = NULL, ...) {
  .contains_draws(object)
  draws <- as_draws_array(object, variable = pars)
  out <- posterior::summarise_draws(draws, ess_bulk = posterior::ess_bulk)
  setNames(out$ess_bulk / posterior::ndraws(draws), out$variable)
}

## --- internal

.draws_diag_to_nuts_params <- function(draws_diag, pars = NULL) {
  if (is.null(draws_diag)) {
    stop("The model does not contain sampler diagnostics.", call. = FALSE)
  }
  draws_diag <- posterior::subset_draws(draws_diag, variable = pars)
  diag_df <- posterior::as_draws_df(draws_diag)
  variables <- posterior::variables(draws_diag)
  out <- lapply(variables, function(variable) {
    data.frame(
      Chain = diag_df$.chain,
      Iteration = diag_df$.iteration,
      Parameter = variable,
      Value = diag_df[[variable]]
    )
  })
  out <- do.call(rbind, out)
  out$Parameter <- factor(out$Parameter, levels = variables)
  rownames(out) <- NULL
  out
}

.contains_draws <- function(object) {
  assert_that(
    nsamples(object) > 0,
    msg = "The model does not contain posterior draws."
  )
}
