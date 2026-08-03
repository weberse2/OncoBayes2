# Code adapted from
# https://raw.githubusercontent.com/paul-buerkner/brms/refs/heads/master/R/posterior.R

#'
#' Transform `blrmfit` or `blrm_trial` to `draws` objects
#'
#' @description
#' Transform a `blrmfit` or `blrm_trial` object to a format supported by the
#' \pkg{posterior} package.
#'
#' @aliases as_draws as_draws_matrix as_draws_array as_draws_df
#' @aliases as_draws_rvars as_draws_list
#'
#' @param x A `blrmfit` or `blrm_trial` object.
#' @param variable A character vector providing the variables to
#'   extract.  By default, all variables are extracted.
#' @param regex Logical; Should variable be treated as a (vector of)
#'   regular expressions? Any variable in `x` matching at least
#'   one of the regular expressions will be selected. Defaults to
#'   `FALSE`.
#' @param inc_warmup Should warmup draws be included? Defaults to
#'   `FALSE`.
#' @param ... Arguments passed to individual methods (if applicable).
#'
#' @details To subset iterations, chains, or draws, use the
#'   [posterior::subset_draws()] method after
#'   transforming the input object to a `draws` object.
#'
#' @seealso [posterior::draws()]
#'   [posterior::subset_draws()]
#'
#' @template start-example
#' @examples
#' # fit an example model. See documentation for "combo2" example
#' example_model("combo2")
#'
#' post <- as_draws(blrmfit)
#'
#' @template stop-example
#'
#' @name draws-OncoBayes2
NULL

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws
#' @method as_draws blrmfit
#' @export
#' @export as_draws
as_draws.blrmfit <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  # draws_list is the fastest format to convert to at the moment
  .as_draws_conversion(
    x,
    as_draws_list,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix blrmfit
#' @export
#' @export as_draws_matrix
as_draws_matrix.blrmfit <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_matrix,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_array
#' @method as_draws_array blrmfit
#' @export
#' @export as_draws_array
as_draws_array.blrmfit <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_array,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_df
#' @method as_draws_df blrmfit
#' @export
#' @export as_draws_df
as_draws_df.blrmfit <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_df,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_list
#' @method as_draws_list blrmfit
#' @export
#' @export as_draws_list
as_draws_list.blrmfit <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_list,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars blrmfit
#' @export
#' @export as_draws_rvars
as_draws_rvars.blrmfit <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .as_draws_conversion(
    x,
    as_draws_rvars,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @keywords internal
.as_draws_conversion <- function(
  x,
  draws_converter,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  if (!is.logical(inc_warmup) || length(inc_warmup) != 1L || is.na(inc_warmup)) {
    stop("`inc_warmup` must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(x$draws)) {
    stop("The model does not contain posterior draws.", call. = FALSE)
  }
  out <- draws_converter(.stored_blrmfit_draws(x, inc_warmup = inc_warmup), ...)
  # subset variables
  subset_draws(out, variable = variable, regex = regex)
}

#' @keywords internal
.stored_blrmfit_draws <- function(x, inc_warmup = FALSE) {
  if (!isTRUE(inc_warmup)) {
    return(x$draws)
  }
  if (is.null(x$draws_warmup)) {
    stop(
      "Warmup draws were not saved. Refit with save_warmup = TRUE ",
      "to use inc_warmup = TRUE.",
      call. = FALSE
    )
  }
  posterior::bind_draws(x$draws_warmup, x$draws, along = "iteration")
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws
#' @method as_draws blrm_trial
#' @export
#' @export as_draws
as_draws.blrm_trial <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws(
    x$blrmfit,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_matrix
#' @method as_draws_matrix blrm_trial
#' @export
#' @export as_draws_matrix
as_draws_matrix.blrm_trial <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_matrix(
    x$blrmfit,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_array
#' @method as_draws_array blrm_trial
#' @export
#' @export as_draws_array
as_draws_array.blrm_trial <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_array(
    x$blrmfit,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_df
#' @method as_draws_df blrm_trial
#' @export
#' @export as_draws_df
as_draws_df.blrm_trial <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_df(
    x$blrmfit,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_list
#' @method as_draws_list blrm_trial
#' @export
#' @export as_draws_list
as_draws_list.blrm_trial <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_list(
    x$blrmfit,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}

#' @rdname draws-OncoBayes2
#' @importFrom posterior as_draws_rvars
#' @method as_draws_rvars blrm_trial
#' @export
#' @export as_draws_rvars
as_draws_rvars.blrm_trial <- function(
  x,
  variable = NULL,
  regex = FALSE,
  inc_warmup = FALSE,
  ...
) {
  .assert_is_blrm_trial_and_prior_is_set(x)
  as_draws_rvars(
    x$blrmfit,
    variable = variable,
    regex = regex,
    inc_warmup = inc_warmup,
    ...
  )
}
