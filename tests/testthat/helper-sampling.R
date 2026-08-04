library(checkmate)
library(assertthat)
library(Formula)
library(abind)
library(rstan)

set_sampling_default <- function(
  iter,
  warmup,
  chains,
  cores = getOption("mc.cores", 1),
  save_warmup = FALSE,
  backend = "rstan",
  control = list()
) {
  options(
    OncoBayes2.MC.iter = iter,
    OncoBayes2.MC.warmup = warmup,
    OncoBayes2.MC.chains = chains,
    mc.cores = cores,
    OncoBayes2.MC.save_warmup = save_warmup,
    OncoBayes2.MC.backend = backend,
    OncoBayes2.MC.control = control
  )
}

very_fast_sampling <- function() {
  message("Tests running with very fast sampling")
  ## note: 250 warmups are needed to get Stans NUTS adaptation to work
  set_sampling_default(500, 250, 1, 1, control = list(adapt_delta = 0.85))
}

fake_sampling <- function() {
  message("Tests running with fake sampling")
  set_sampling_default(10, 2, 1, 1, control = list(adapt_delta = 0.85))
}

default_sampling <- function() {
  set_sampling_default(2000, 1000, 4, 1)
}

fast_sampling <- function() {
  message("Tests running with fast sampling")
  set_sampling_default(1000, 500, 1, 1, control = list(adapt_delta = 0.85))
}

run_example <- function(example) {
  env <- new.env()
  suppressWarnings(example_model(example, env, silent = TRUE))
  invisible(env)
}

## Skip a test unless both the cmdstanr package and a usable CmdStan
## installation are available. skip_if_not_installed("cmdstanr") only checks
## for the R package, but the cmdstanr backend additionally needs a compiled
## CmdStan toolchain (absent on e.g. CI), without which model compilation
## fails with "CmdStan path has not been set yet".
skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  cmdstan_available <- tryCatch(
    {
      ver <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
      !is.null(ver) && !is.na(ver)
    },
    error = function(e) FALSE
  )
  if (!isTRUE(cmdstan_available)) {
    testthat::skip("CmdStan toolchain not available")
  }
}


## set up slim sampling in case we are on CRAN
if (identical(Sys.getenv("NOT_CRAN"), "true")) {
  very_fast_sampling()
} else {
  fake_sampling()
}


## takes a flat array of a mixture mvn and returns an RBesT mixture
## mvn
array2mix <- function(a, p) {
  Nc <- dim(a)[1]
  l <- dim(a)[2]
  mix <- list()
  for (i in seq_len(Nc)) {
    w <- a[i, 1]
    m <- a[i, 2:(2 + p - 1)]
    s <- a[i, (2 + p):(2 + p + p - 1)]
    if (p > 1) {
      rho <- a[i, (2 + p + p):l]
      R <- diag(1, p, p)
      R[lower.tri(R)] <- rho
      R[upper.tri(R)] <- t(R)[upper.tri(R)]
    } else {
      R <- diag(1, 1, 1)
    }
    mix[[i]] <- unname(c(w, m, diag(s, nrow = p) %*% R %*% diag(s, nrow = p)))
  }
  do.call(RBesT::mixmvnorm, mix)
}

## takes a fitted blrmfit object and returns the blrmfit object with a
## posterior containing one draw with all slopes and intercept set to
## the prior mean (exactly).
sample_prior_mean <- function(blrmfit) {
  ps <- prior_summary(blrmfit)

  num_comp <- dim(ps$EX_mu_log_beta)[2]
  has_inter <- ps$has_inter
  num_inter <- 0
  if (has_inter) {
    num_inter <- length(grep("^m", dimnames(ps$EX_mu_eta)$prior))
  }

  EX_m_intercept <- apply(
    ps$EX_mu_log_beta["weight", , , drop = FALSE] *
      ps$EX_mu_log_beta["m_intercept", , , drop = FALSE],
    c(1, 2),
    sum
  )
  EX_m_log_slope <- apply(
    ps$EX_mu_log_beta["weight", , , drop = FALSE] *
      ps$EX_mu_log_beta["m_log_slope", , , drop = FALSE],
    c(1, 2),
    sum
  )

  variables <- character()
  values <- numeric()
  for (group in levels(blrmfit$group_fct)) {
    for (component_id in seq_len(num_comp)) {
      component <- levels(blrmfit$labels$component)[component_id]
      variables <- c(
        variables,
        paste0("beta_group[", group, ",", component, ",intercept]"),
        paste0("beta_group[", group, ",", component, ",slope]")
      )
      values <- c(
        values,
        EX_m_intercept[component_id],
        exp(EX_m_log_slope[component_id])
      )
    }
  }
  if (has_inter) {
    EX_mu_eta <- array(
      summary(array2mix(ps$EX_mu_eta, num_inter))$mean,
      num_inter
    )
    for (group in levels(blrmfit$group_fct)) {
      for (inter_id in seq_len(num_inter)) {
        inter <- levels(blrmfit$labels$param_eta)[inter_id]
        variables <- c(
          variables,
          paste0("eta_group[", group, ",", inter, "]")
        )
        values <- c(
          values,
          EX_mu_eta[inter_id]
        )
      }
    }
  }

  blrmfit$draws <- posterior::as_draws_array(
    array(
      values,
      dim = c(1, 1, length(values)),
      dimnames = list(
        NULL,
        "chain:1",
        variables
      )
    )
  )
  blrmfit$draws_warmup <- NULL
  blrmfit$draws_diag <- NULL
  blrmfit$draws_warmup_diag <- NULL
  blrmfit$metadata_mcmc <- modifyList(
    blrmfit$metadata_mcmc,
    list(
      iter = posterior::niterations(blrmfit$draws),
      warmup = 0L,
      chains = posterior::nchains(blrmfit$draws),
      save_warmup = FALSE,
      algorithm = "Fixed_param"
    ),
    keep.null = TRUE
  )
  blrmfit$stanfit <- NULL
  blrmfit
}
