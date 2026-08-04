#' @keywords internal
.blrmfit_draws_from_backend <- function(
  backend,
  stanfit,
  cmdstanfit = NULL,
  variable_label_specs = NULL,
  iter,
  warmup,
  chains,
  thin,
  save_warmup,
  exclude_variables = character()
) {
  fit <- switch(
    backend,
    rstan = stanfit,
    cmdstanr = cmdstanfit,
    stop("Unknown MCMC backend: ", backend, call. = FALSE)
  )
  draws <- posterior::as_draws_array(fit)
  draws <- .exclude_draws_variables(draws, exclude_variables)
  draws <- .relabel_draws_variables(draws, variable_label_specs)
  draws_warmup <- NULL
  if (isTRUE(save_warmup)) {
    draws_all <- .as_draws_array_with_warmup(fit, backend)
    draws_all <- .exclude_draws_variables(draws_all, exclude_variables)
    draws_all <- .relabel_draws_variables(draws_all, variable_label_specs)
    draws_warmup <- .split_warmup_draws(draws_all, draws)
  }

  draws_diag <- .as_draws_diag(fit, backend, inc_warmup = FALSE)
  draws_warmup_diag <- NULL
  if (isTRUE(save_warmup)) {
    draws_diag_all <- .as_draws_diag(fit, backend, inc_warmup = TRUE)
    draws_warmup_diag <- .split_warmup_draws(draws_diag_all, draws_diag)
  }

  list(
    draws = draws,
    draws_warmup = draws_warmup,
    draws_diag = draws_diag,
    draws_warmup_diag = draws_warmup_diag,
    metadata_mcmc = list(
      iter = iter,
      warmup = warmup,
      chains = chains,
      thin = thin,
      save_warmup = save_warmup,
      algorithm = "NUTS",
      backend = backend,
      backend_versions = .mcmc_backend_versions(backend)
    ),
    backend = backend
  )
}

#' @keywords internal
.exclude_draws_variables <- function(draws, exclude_variables = character()) {
  if (!length(exclude_variables)) {
    return(draws)
  }
  variables <- posterior::variables(draws)
  exclude_pattern <- paste0(
    "^(",
    paste(
      gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", exclude_variables),
      collapse = "|"
    ),
    ")(\\[|$)"
  )
  keep_variables <- variables[!grepl(exclude_pattern, variables)]
  posterior::subset_draws(draws, variable = keep_variables)
}

#' @keywords internal
.as_draws_array_with_warmup <- function(fit, backend) {
  switch(
    backend,
    rstan = posterior::as_draws_array(rstan::extract(
      fit,
      permuted = FALSE,
      inc_warmup = TRUE
    )),
    cmdstanr = fit$draws(format = "draws_array", inc_warmup = TRUE),
    stop("Unknown MCMC backend: ", backend, call. = FALSE)
  )
}

#' @keywords internal
.relabel_draws_variables <- function(draws, variable_label_specs) {
  if (is.null(variable_label_specs)) {
    return(draws)
  }
  variables <- posterior::variables(draws)
  for (spec in variable_label_specs) {
    variables <- do.call(
      .label_index_names,
      c(list(x = variables, par = spec$par), spec$factors)
    )
  }
  posterior::variables(draws) <- variables
  draws
}

#' @keywords internal
.label_index_names <- function(x, par, ...) {
  idx <- grep(paste0("^", par, "\\["), x)
  if (!length(idx)) {
    return(x)
  }

  fct <- list(...)
  idx_str <- t(sapply(
    strsplit(gsub("(.*)\\[([0-9,]*)\\]$", "\\2", x[idx]), ","),
    as.numeric
  ))
  if (length(fct) == 1) {
    idx_str <- matrix(idx_str, ncol = 1)
  }
  ni <- ncol(idx_str)
  stopifnot(ni == length(fct))

  label_str <- idx_str
  for (i in seq_len(ni)) {
    label_str[, i] <- levels(fct[[i]])[idx_str[, i]]
  }
  x[idx] <- paste0(
    par,
    "[",
    do.call(paste, c(as.data.frame(label_str), list(sep = ","))),
    "]"
  )
  x
}

#' @keywords internal
.as_draws_diag <- function(fit, backend, inc_warmup = FALSE) {
  if (!isTRUE(inc_warmup)) {
    return(.nuts_params_to_draws_array(bayesplot::nuts_params(fit)))
  }
  switch(
    backend,
    rstan = .nuts_params_to_draws_array(bayesplot::nuts_params(
      fit,
      inc_warmup = TRUE
    )),
    cmdstanr = fit$sampler_diagnostics(
      format = "draws_array",
      inc_warmup = TRUE
    ),
    stop("Unknown MCMC backend: ", backend, call. = FALSE)
  )
}

#' @keywords internal
.nuts_params_to_draws_array <- function(nuts_params) {
  stopifnot(all(c("Chain", "Iteration", "Parameter", "Value") %in%
    names(nuts_params)))
  chains <- sort(unique(nuts_params$Chain))
  iterations <- sort(unique(nuts_params$Iteration))
  variables <- levels(nuts_params$Parameter)
  if (is.null(variables)) {
    variables <- unique(as.character(nuts_params$Parameter))
  }
  x <- array(
    NA_real_,
    dim = c(length(iterations), length(chains), length(variables)),
    dimnames = list(
      NULL,
      paste0("chain:", chains),
      variables
    )
  )
  chain_id <- match(nuts_params$Chain, chains)
  iteration_id <- match(nuts_params$Iteration, iterations)
  variable_id <- match(as.character(nuts_params$Parameter), variables)
  x[cbind(iteration_id, chain_id, variable_id)] <- nuts_params$Value
  if (anyNA(x)) {
    stop("Incomplete NUTS diagnostic table.", call. = FALSE)
  }
  posterior::as_draws_array(x)
}

#' @keywords internal
.split_warmup_draws <- function(draws_all, draws_sampling) {
  n_warmup <- posterior::niterations(draws_all) -
    posterior::niterations(draws_sampling)
  stopifnot(n_warmup > 0)
  posterior::subset_draws(draws_all, iteration = seq_len(n_warmup))
}

#' @keywords internal
.mcmc_backend_versions <- function(backend) {
  versions <- list()
  if (requireNamespace("posterior", quietly = TRUE)) {
    versions$posterior <- as.character(utils::packageVersion("posterior"))
  }
  if (requireNamespace("rstan", quietly = TRUE)) {
    versions$rstan <- as.character(utils::packageVersion("rstan"))
  }
  if (identical(backend, "cmdstanr") && requireNamespace("cmdstanr", quietly = TRUE)) {
    versions$cmdstanr <- as.character(utils::packageVersion("cmdstanr"))
    versions$cmdstan <- tryCatch(
      cmdstanr::cmdstan_version(),
      error = function(e) NA_character_
    )
  }
  versions
}
