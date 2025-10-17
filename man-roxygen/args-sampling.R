#' @param iter number of iterations (including warmup).
#' @param warmup number of warmup iterations.
#' @param save_warmup save warmup samples (\code{TRUE} /
#'     \code{FALSE}). Only if set to \code{TRUE}, then all random
#'     variables are saved in the posterior. This substantially
#'     increases the storage needs of the posterior.
#' @param thin period of saving samples.
#' @param init positive number to specify uniform range on
#'     unconstrained space for random initialization. See
#'     \code{\link[rstan:stan]{stan}}.
#' @param chains number of Markov chains.
#' @param cores number of cores for parallel sampling of chains.
#' @param control additional sampler parameters for NuTS algorithm.
#' @param backend sets Stan backend to be used. Possible choices are
#'     \code{"rstan"} (default) or \code{"cmdstanr"}.
