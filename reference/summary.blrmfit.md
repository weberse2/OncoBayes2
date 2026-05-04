# Summarise model results

Provides model summaries for
[`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md)
and
[`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
analyses.

## Usage

``` r
# S3 method for class 'blrmfit'
summary(
  object,
  newdata,
  transform = !predictive,
  prob = 0.95,
  interval_prob,
  predictive = FALSE,
  ...
)
```

## Arguments

- object:

  fitted model object

- newdata:

  optional data frame specifying for what to predict; if missing, then
  the data of the input model `object` is used

- transform:

  logical (defaults to `FALSE`) indicating if the linear predictor on
  the logit link scale is transformed with `inv_logit` to the 0-1
  response scale.

- prob:

  central probability mass to report, i.e. the quantiles 0.5-prob/2 and
  0.5+prob/2 are displayed. Multiple central widths can be specified.

- interval_prob:

  optional vector of sorted quantiles for which the interval
  probabilities are calculated

- predictive:

  logical indicates if the posterior predictive is being summarized.
  Defaults to `FALSE`.

- ...:

  not used in this function

## Value

Returns a `data.frame` of the key summaries of the posterior mean,
standard deviation, central probability interval, median and optional
interval probabilities. Each row of the `data.frame` corresponds to the
respective input data which is by default the same data set as used for
the
[`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md)
analysis or the data specified in the `newdata` argument.

## Details

The calculated posterior summaries are returned as a `data.frame` and
contain optional interval probabilites for the specified vector of
sorted quantiles. These summaries are calculated on the response scale
by default and can be obtained on the link scale when setting
`transform=FALSE`.

When the results are requested for the predictive distribution with
`predictive=TRUE`, then the link scale refers to the total counts while
the transformed scale divides the (predictive) counts by the number of
trials such that results are on the 0-1 scale.

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

example_model("single_agent", silent = TRUE)
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

## obtain underdosing (0-0.16], target dosing (0.16-0.33] and
## overdosing (0.33-1] probabilities
summary(blrmfit, interval_prob = c(0, 0.16, 0.33, 1))
#>          mean          sd         2.5%          50%      97.5% [0,0.16]
#> 1 0.000666800 0.001343123 1.083755e-07 0.0003054777 0.00357141      1.0
#> 2 0.004758683 0.006701482 8.491255e-06 0.0034943188 0.01880681      1.0
#> 3 0.023789085 0.024556450 2.342204e-04 0.0224469994 0.06749307      1.0
#> 4 0.120178943 0.101934619 6.572024e-03 0.1020491566 0.27082706      0.6
#> 5 0.637061758 0.145485587 4.716297e-01 0.6164261785 0.84390315      0.0
#>   (0.16,0.33] (0.33,1]
#> 1         0.0        0
#> 2         0.0        0
#> 3         0.0        0
#> 4         0.4        0
#> 5         0.0        1

## obtain predictive distribution for respective cohorts and
## calculate probability for no event, 1 event or >1 event
## note that this does the calculation for the cohort sizes
## as put into the data-set
summary(blrmfit, interval_prob = c(-1, 0, 1, 10), predictive = TRUE)
#>         mean         sd 2.5% 50% 97.5%    (-1,0]       (0,1]       (1,10]
#> 1 0.00200040 0.04481972    0   0     0 0.9980058 0.001988017 6.187214e-06
#> 2 0.01903473 0.13938859    0   0     0 0.9813390 0.018291972 3.690539e-04
#> 3 0.11894542 0.35632873    0   0     1 0.8915458 0.098539818 9.914384e-03
#> 4 0.48071577 0.73154840    0   0     2 0.6418672 0.254452634 1.036802e-01
#> 5 1.27412352 0.70747936    0   1     2 0.1507736 0.424329248 4.248971e-01

## to obtain the predictive for a cohort-size of 6 for all patients
## in the data-set one would need to use the newdata argument, e.g.
summary(blrmfit,
  newdata = transform(hist_SA, num_patients = 6),
  interval_prob = c(-1, 0, 1, 10), predictive = TRUE
)
#>        mean         sd 2.5% 50% 97.5%      (-1,0]       (0,1]       (1,10]
#> 1 0.0040008 0.06361477    0   0     0 0.996030050 0.003939273 3.067682e-05
#> 2 0.0285521 0.17213017    0   0     1 0.972370719 0.026729239 9.000420e-04
#> 3 0.1427345 0.39448765    0   0     1 0.872722591 0.112939247 1.433816e-02
#> 4 0.7210737 0.95653739    0   0     3 0.545916684 0.261360807 1.927225e-01
#> 5 3.8223706 1.39955992    1   4     6 0.007521125 0.046926422 9.455525e-01

## Recover user set sampling defaults
options(.user_mc_options)
```
