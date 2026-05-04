# Bayesian Logistic Regression Model for N-compounds with EXNEX

Bayesian Logistic Regression Model (BLRM) for N compounds using
EXchangability and NonEXchangability (EXNEX) modeling.

## Usage

``` r
blrm_exnex(
  formula,
  data,
  prior_EX_mu_comp,
  prior_EX_mu_mean_comp,
  prior_EX_mu_sd_comp,
  prior_EX_tau_comp,
  prior_EX_tau_mean_comp,
  prior_EX_tau_sd_comp,
  prior_EX_corr_eta_comp,
  prior_EX_mu_inter,
  prior_EX_mu_mean_inter,
  prior_EX_mu_sd_inter,
  prior_EX_tau_inter,
  prior_EX_tau_mean_inter,
  prior_EX_tau_sd_inter,
  prior_EX_corr_eta_inter,
  prior_is_EXNEX_inter,
  prior_is_EXNEX_comp,
  prior_EX_prob_comp,
  prior_EX_prob_inter,
  prior_NEX_mu_comp,
  prior_NEX_mu_mean_comp,
  prior_NEX_mu_sd_comp,
  prior_NEX_mu_inter,
  prior_NEX_mu_mean_inter,
  prior_NEX_mu_sd_inter,
  prior_tau_dist,
  sample_map = FALSE,
  iter = getOption("OncoBayes2.MC.iter", 2000),
  warmup = getOption("OncoBayes2.MC.warmup", 1000),
  save_warmup = getOption("OncoBayes2.MC.save_warmup", TRUE),
  thin = getOption("OncoBayes2.MC.thin", 1),
  init = getOption("OncoBayes2.MC.init", 0.5),
  chains = getOption("OncoBayes2.MC.chains", 4),
  cores = getOption("mc.cores", 1L),
  control = getOption("OncoBayes2.MC.control", list()),
  backend = getOption("OncoBayes2.MC.backend", "rstan"),
  prior_PD = FALSE,
  verbose = FALSE
)

# S3 method for class 'blrmfit'
print(x, ..., prob = 0.95, digits = 2)
```

## Arguments

- formula:

  the model formula describing the linear predictors of the model. The
  lhs of the formula is a two-column matrix which are the number of
  occured events and the number of times no event occured. The rhs of
  the formula defines the linear predictors for the marginal models for
  each drug component, then the interaction model and at last the
  grouping and optional stratum factors of the models. These elements of
  the formula are separated by a vertical bar. The marginal models must
  follow a intercept and slope form while the interaction model must not
  include an interaction term. See the examples below for an example
  instantiation.

- data:

  optional data frame containing the variables of the model. If not
  found in `data`, the variables are taken from `environment(formula)`.

- prior_EX_mu_comp:

  List of bivariate normal mixture priors for intercept and slope
  parameters \\\boldsymbol\mu_i = (\mu\_{\alpha i}, \mu\_{\beta i})\\ of
  each component. In case of a single drug model, then a mixture prior
  is accepted as well.

- prior_EX_mu_mean_comp, prior_EX_mu_sd_comp:

  **\[deprecated\]** Please use `prior_EX_mu_comp` instead. Mean and sd
  for the prior on the mean parameters \\\boldsymbol\mu_i =
  (\mu\_{\alpha i}, \mu\_{\beta i})\\ of each component. Two column
  matrix (intercept, log-slope) with one row per component.

- prior_EX_tau_comp:

  List of bivariate normal mixture priors for heterogeniety parameter
  \\\boldsymbol\tau\_{si} = (\tau\_{\alpha s i}, \tau\_{\beta s i})\\ of
  each stratum. If no differential discounting is required (i.e. if
  there is only one stratum \\s = 1\\), then it suffices to provide a
  bivariate normal mixture prior instead of a list with just one
  element.

- prior_EX_tau_mean_comp, prior_EX_tau_sd_comp:

  **\[deprecated\]** Please use `prior_EX_tau_comp` instead. Prior mean
  and sd for heterogeniety parameter \\\boldsymbol\tau\_{si} =
  (\tau\_{\alpha s i}, \tau\_{\beta s i})\\ of each stratum. If no
  differential discounting is required (i.e. if there is only one
  stratum \\s = 1\\), then it is a two-column matrix (intercept,
  log-slope) with one row per component. Otherwise it is a
  three-dimensional array whose first dimension indexes the strata,
  second dimension indexes the components, and third dimension of length
  two for (intercept, log-slope).

- prior_EX_corr_eta_comp:

  Prior LKJ correlation parameter for each component given as numeric
  vector. If missing, then a 1 is assumed corresponding to a marginal
  uniform prior of the correlation.

- prior_EX_mu_inter:

  Multivariate normal mixture prior for interaction parameter vector
  \\\boldsymbol{\mu\_{\eta}}\\. Dimension must correspond to the number
  of interactions.

- prior_EX_mu_mean_inter, prior_EX_mu_sd_inter:

  **\[deprecated\]** Please use `prior_EX_mu_inter` instead. Prior mean
  and sd for population mean parameters \\\mu\_{\eta k}\\ of each
  interaction parameter. Vector of length equal to the number of
  interactions.

- prior_EX_tau_inter:

  List of multivariate normal mixture priors for heterogeniety
  interaction parameter vector \\\boldsymbol{\tau\_{\eta s}}\\ of each
  stratum. If no differential discounting is required (i.e. if there is
  only one stratum \\s = 1\\), then it suffices to provide a mixture
  prior instead of a list with just one element.

- prior_EX_tau_mean_inter, prior_EX_tau_sd_inter:

  **\[deprecated\]** Please use `prior_EX_tau_inter` instead. Prior mean
  and sd for heterogeniety parameter \\\tau\_{\eta s k}\\ of each
  stratum. Matrix with one column per interaction and one row per
  stratum.

- prior_EX_corr_eta_inter:

  Prior LKJ correlation parameter for interaction given as numeric. If
  missing, then a 1 is assumed corresponding to a marginal uniform prior
  of the correlations.

- prior_is_EXNEX_inter:

  Defines if non-exchangability is admitted for a given interaction
  parameter. Logical vector of length equal to the number of
  interactions. If missing `FALSE` is assumed for all interactions.

- prior_is_EXNEX_comp:

  Defines if non-exchangability is admitted for a given component.
  Logical vector of length equal to the number of components. If missing
  `TRUE` is assumed for all components.

- prior_EX_prob_comp:

  Prior probability \\p\_{ij}\\ for exchangability of each component per
  group. Matrix with one column per component and one row per group.
  Values must lie in \\\[0-1\]\\ range.

- prior_EX_prob_inter:

  Prior probability \\p\_{\eta k j}\\ for exchangability of each
  interaction per group. Matrix with one column per interaction and one
  row per group. Values must lie in \\\[0-1\]\\ range.

- prior_NEX_mu_comp:

  List of bivariate normal mixture priors \\\boldsymbol m\_{ij}\\ and
  \\\boldsymbol s\_{ij} = \text{diag}(\boldsymbol S\_{ij})\\ of each
  component for non-exchangable case. If missing set to the same prior
  as given for the EX part. It is required that the specification be the
  same across groups j.

- prior_NEX_mu_mean_comp, prior_NEX_mu_sd_comp:

  **\[deprecated\]** Please use `prior_NEX_mu_comp` instead. Prior mean
  \\\boldsymbol m\_{ij}\\ and sd \\\boldsymbol s\_{ij} =
  \text{diag}(\boldsymbol S\_{ij})\\ of each component for
  non-exchangable case. Two column matrix (intercept, log-slope) with
  one row per component. If missing set to the same prior as given for
  the EX part. It is required that the specification be the same across
  groups j.

- prior_NEX_mu_inter:

  Multivariate normal mixture prior (mean \\m\_{\eta k j}\\, sd
  \\s\_{\eta k j}\\ and covariance) for the interaction parameter vector
  for non-exchangable case. Dimension must correspond to the number of
  interactions. If missing set to the same prior as given for the EX
  part.

- prior_NEX_mu_mean_inter, prior_NEX_mu_sd_inter:

  **\[deprecated\]** Please use `prior_NEX_mu_inter` instead. Prior mean
  \\m\_{\eta k j}\\ and sd \\s\_{\eta k j}\\ for each interaction
  parameter for non-exchangable case. Vector of length equal to the
  number of interactions. If missing set to the same prior as given for
  the EX part.

- prior_tau_dist:

  Defines the distribution used for heterogeniety parameters. Choices
  are 0=fixed to it's mean, 1=log-normal, 2=truncated normal or `NULL`
  shutting off the hierarchical structure of the model.

- sample_map:

  Logical flag (defaults to `FALSE`) controlling inclusion of MAP priors
  for each stratum defined as part of the generated posterior. If set to
  `TRUE` then the posterior samples will contain `map_log_beta` and
  `map_eta` variables.

- iter:

  number of iterations (including warmup).

- warmup:

  number of warmup iterations.

- save_warmup:

  save warmup samples (`TRUE` / `FALSE`). Only if set to `TRUE`, then
  all random variables are saved in the posterior. This substantially
  increases the storage needs of the posterior.

- thin:

  period of saving samples.

- init:

  positive number to specify uniform range on unconstrained space for
  random initialization. See
  [`stan`](https://mc-stan.org/rstan/reference/stan.html).

- chains:

  number of Markov chains.

- cores:

  number of cores for parallel sampling of chains.

- control:

  additional sampler parameters for NuTS algorithm.

- backend:

  sets Stan backend to be used. Possible choices are `"rstan"` (default)
  or `"cmdstanr"`.

- prior_PD:

  Logical flag (defaults to `FALSE`) indicating if to sample the prior
  predictive distribution instead of conditioning on the data.

- verbose:

  Logical flag (defaults to `FALSE`) controlling if additional output
  like stan progress is reported.

- x:

  `blrmfit` object to print

- ...:

  not used in this function

- prob:

  central probability mass to report, i.e. the quantiles 0.5-prob/2 and
  0.5+prob/2 are displayed. Multiple central widths can be specified.

- digits:

  number of digits to show

## Value

The function returns a S3 object of type `blrmfit`.

## Details

`blrm_exnex` is a flexible function for Bayesian meta-analytic modeling
of binomial count data. In particular, it is designed to model counts of
the number of observed dose limiting toxicities (DLTs) by dose, for
guiding dose-escalation studies in Oncology. To accommodate dose
escalation over more than one agent, the dose may consist of
combinations of study drugs, with any number of treatment components.

In the simplest case, the aim is to model the probability \\\pi\\ that a
patient experiences a DLT, by complementing the binomial likelihood with
a monotone logistic regression

\$\$\text{logit}\\\pi(d) = \log\\\alpha + \beta \\ t(d),\$\$

where \\\beta \> 0\\. Most typically, \\d\\ represents the dose, and
\\t(d)\\ is an appropriate transformation, such as \\t(d) = \log (d \big
/ d^\*)\\. A joint prior on \\\boldsymbol \theta = (\log\\\alpha,
\log\\\beta)\\ completes the model and ensures monotonicity \\\beta \>
0\\.

Many extensions are possible. The function supports general combination
regimens, and also provides framework for Bayesian meta-analysis of
dose-toxicity data from multiple historical and concurrent sources.

For an example of a single-agent trial refer to
[`example-single-agent()`](https://opensource.nibr.com/OncoBayes2/reference/example-single-agent.md).

## Functions

- `print(blrmfit)`: print function.

## Combination of two treatments

For a combination of two treatment components, the basic modeling
framework is that the DLT rate \\\pi(d_1,d_2)\\ is comprised of (1) a
"no-interaction" baseline model \\\tilde \pi(d_1,d_2)\\ driven by the
single-agent toxicity of each component, and (2) optional interaction
terms \\\gamma(d_1,d_2)\\ representing synergy or antagonism between the
drugs. On the log-odds scale,

\$\$\text{logit} \\\pi(d_1,d_2) = \text{logit} \\ \tilde \pi(d_1,d_2) +
\eta \\ \gamma(d_1,d_2). \$\$

The "no interaction" part \\\tilde \pi(d_1,d_2)\\ represents the
probability of a DLT triggered by either treatment component acting
*independently*. That is, \$\$ \tilde \pi(d_1,d_2) = 1- (1 -
\pi_1(d_1))(1 - \pi_2(d_2)). \$\$ In simple terms, P(no DLT for
combination) = P(no DLT for drug 1) \* P(no DLT from drug 2). To
complete this part, the treatment components can then be modeled with
monotone logistic regressions as before.

\$\$\text{logit} \\ \pi_i(d_i) = \log\\ \alpha_i + \beta_i \\
t(d_i),\$\$

where \\t(d_i)\\ is a monotone transformation of the doses of the
respective drug component \$i\$, such as \\t(d_i) = \log (d_i \big /
d_i^\*)\\.

The inclusion of an interaction term \\\gamma(d_1,d_2)\\ allows DLT
rates above or below the "no-interaction" rate. The magnitude of the
interaction term may also be made dependent on the doses (or other
covariates) through regression. As an example, one could let

\$\$\gamma(d_1, d_2) = \frac{d_1}{d_1^\*} \frac{d_1}{d_2^\*}.\$\$

The specific functional form is specified in the usual notation for a
design matrix. The interaction model must respect the constraint that
whenever any dose approaches zero, then the interaction term must vanish
as well. Therefore, the interaction model must not include an intercept
term which would violate this consistency requirement. A dual
combination example can be found in
[`example-combo2()`](https://opensource.nibr.com/OncoBayes2/reference/example-combo2.md).

## General combinations

The model is extended to general combination treatments consisting of
\\N\\ components by expressing the probability \\\pi\\ on the logit
scale as

\$\$ \text{logit} \\ \pi(d_1,\ldots,d_N) = \text{logit} \Bigl( 1 -
\prod\_{i = 1}^N ( 1 - \pi_i(d_i) ) \Bigr) + \sum\_{k=1}^K \eta_k \\
\gamma_k(d_1,\ldots,d_N), \$\$

Multiple drug-drug interactions among the \\N\\ components are now
possible, and are represented through the \\K\\ interaction terms
\\\gamma_k\\.

Regression models can be again be specified for each \\\pi_i\\ and
\\\gamma_k\\, such as

\$\$ \text{logit}\\ \pi_i(d_i) = \log\\ \alpha_i + \beta_i \\ t(d_i)
\$\$

Interactions for some subset \\I(k) \subset \\1,\ldots,N \\\\ of the
treatment components can be modeled with regression as well, for example
on products of doses,

\$\$ \gamma_k(d_1,\ldots,d_N) = \prod\_{i \in I(k)}
\frac{d_i}{d_i^\*}.\$\$

For example, \\I(k) = \\1,2,3\\\\ results in the three-way interaction
term

\$\$ \frac{d_1}{d_1^\*} \frac{d_2}{d_2^\*} \frac{d_3}{d_3^\*} \$\$

for drugs 1, 2, and 3.

For a triple combination example please refer to
[`example-combo3()`](https://opensource.nibr.com/OncoBayes2/reference/example-combo3.md).

## Meta-analytic framework

Information on the toxicity of a drug may be available from multiple
studies or sources. Furthermore, one may wish to stratify observations
within a single study (for example into groups of patients corresponding
to different geographic regions, or multiple dosing `dose_info`
corresponding to different schedules).

`blrm_exnex` provides tools for robust Bayesian hierarchical modeling to
jointly model data from multiple sources. An additional index \\j=1,
\ldots, J\\ on the parameters and observations denotes the \\J\\ groups.
The resulting model allows the DLT rate to differ across the groups. The
general \\N\\-component model becomes

\$\$ \text{logit} \\ \pi_j(d_1,\ldots,d_N) = \text{logit} \Bigl( 1 -
\prod\_{i = 1}^N ( 1 - \pi\_{ij}(d_i) ) \Bigr) + \sum\_{k=1}^K
\eta\_{kj} \\ \gamma\_{k}(d_1,\ldots,d_N), \$\$

for groups \\j = 1,\ldots,J\\. The component toxicities \\\pi\_{ij}\\
and interaction terms \\\gamma\_{k}\\ are modelled, as before, through
regression. For example, \\\pi\_{ij}\\ could be a logistic regression on
\\t(d_i) = \log(d_i/d_i^\*)\\ with intercept and log-slope \\\boldsymbol
\theta\_{ij}\\, and \\\gamma\_{k}\\ regressed with coefficient
\\\eta\_{kj}\\ on a product \\\prod\_{i\in I(k)} (d_i/d_i^\*)\\ for some
subset \\I(k)\\ of components.

Thus, for \\j=1,\ldots,J\\, we now have group-specific parameters
\\\boldsymbol\theta\_{ij} = (\log\\ \alpha\_{ij}, \log\\ \beta\_{ij})\\
and \\\boldsymbol\nu\_{j} = (\eta\_{1j}, \ldots, \eta\_{Kj})\\ for each
component \\i=1,\ldots,N\\ and interaction \\k=1,\ldots,K\\.

The structure of the prior on
\\(\boldsymbol\theta\_{i1},\ldots,\boldsymbol\theta\_{iJ})\\ and
\\(\boldsymbol\nu\_{1}, \ldots, \boldsymbol\nu\_{J})\\ determines how
much information will be shared across groups \\j\\. Several modeling
choices are available in the function.

- *EX (Full exchangeability):* One can assume the parameters are
  conditionally exchangeable given hyperparameters

  \$\$\boldsymbol \theta\_{ij} \sim \text{N}\bigl( \boldsymbol
  \mu\_{\boldsymbol \theta i}, \boldsymbol \Sigma\_{\boldsymbol \theta
  i} \bigr), \$\$

  independently across groups \\j = 1,\ldots, J\\ and treatment
  components \\i=1,\ldots,N\\. The covariance matrix \\\boldsymbol
  \Sigma\_{\boldsymbol \theta i}\\ captures the patterns of cross-group
  heterogeneity, and is parametrized with standard deviations
  \\\boldsymbol \tau\_{\boldsymbol\theta i} = (\tau\_{\alpha i},
  \tau\_{\beta i})\\ and the correlation \\\rho_i\\. Similarly for the
  interactions, the fully-exchangeable model is

  \$\$\boldsymbol \nu\_{j} \sim \text{N}\bigl( \boldsymbol
  \mu\_{\boldsymbol \nu}, \boldsymbol \Sigma\_{\boldsymbol \nu}
  \bigr)\$\$

  for groups \\j = 1,\ldots, J\\ and interactions \\k=1,\ldots,K\\, and
  the prior on the covariance matrix \\\boldsymbol \Sigma\_{\boldsymbol
  \nu}\\ captures the amount of heterogeneity expected in the
  interaction terms a-priori. The covariance is again parametrized with
  standard deviations \\(\tau\_{\eta 1}, \ldots, \tau\_{\eta K})\\ and
  its correlation matrix.

- *Differential discounting:* For one or more of the groups
  \\j=1,\ldots,J\\, larger deviations of \\\boldsymbol\theta\_{ij}\\ may
  be expected from the mean \\\boldsymbol\mu_i\\, or of the interactions
  \\\eta\_{kj}\\ from the mean \\\mu\_{\eta,k}\\. Such differential
  heterogeneity can be modeled by mapping the groups \\j = 1,\ldots,J\\
  to *strata* through \\s_j \in \\1,\ldots,S\\\\, and modifying the
  model specification to \$\$\boldsymbol \theta\_{ij} \sim
  \text{N}\bigl( \boldsymbol \mu\_{\boldsymbol \theta i}, \boldsymbol
  \Sigma\_{\boldsymbol \theta ij} \bigr), \$\$ where \$\$\boldsymbol
  \Sigma\_{\boldsymbol \theta ij} = \left( \begin{array}{cc}
  \tau^2\_{\alpha s_j i} & \rho_i \tau\_{\alpha s_j i} \tau\_{\beta s_j
  i}\\ \rho_i \tau\_{\alpha s_j i} \tau\_{\beta s_j i} & \tau^2\_{\beta
  s_j i} \end{array} \right).\$\$ For the interactions, the model
  becomes \$\$\boldsymbol \nu\_{j} \sim \text{N}\bigl( \boldsymbol
  \mu\_{\boldsymbol \nu}, \boldsymbol \Sigma\_{\boldsymbol \nu j}
  \bigr),\$\$ where the covariance matrix \\\boldsymbol
  \Sigma\_{\boldsymbol \nu j}\\ is modelled as stratum specific standard
  deviations \\(\tau\_{\eta 1 s_j}, \ldots, \tau\_{\eta K s_j})\\ and a
  stratum independent correlation matrix. Each stratum \\s=1,\ldots,S\\
  then corresponds to its own set of standard deviations \\\tau\\
  leading to different discounting per stratum. Independent priors are
  specified for the component parameters \\\tau\_{\alpha s i}\\ and
  \\\tau\_{\beta s i}\\ and for the interaction parameters \\\tau\_{\eta
  s k}\\ for each stratum \\s=1,\ldots,S\\. Inference for strata \\s\\
  where the prior is centered on larger values of the \\\tau\\
  parameters will exhibit less shrinkage towards the the means,
  \\\boldsymbol\mu\_{\boldsymbol \theta i}\\ and \\\boldsymbol
  \mu\_{\boldsymbol \nu}\\ respectively.

- *EXNEX (Partial exchangeability):* Another mechansim for increasing
  robustness is to introduce mixture priors for the group-specific
  parameters, where one mixture component is shared across groups, and
  the other is group-specific. The result, known as an
  EXchangeable-NonEXchangeable (EXNEX) type prior, has a form

  \$\$\boldsymbol \theta\_{ij} \sim p\_{\boldsymbol \theta ij}\\
  \text{N}\bigl( \boldsymbol \mu\_{\boldsymbol \theta i}, \boldsymbol
  \Sigma\_{\boldsymbol \theta i} \bigr) +(1-p\_{\boldsymbol \theta
  ij})\\ \text{N}\bigl(\boldsymbol m\_{\boldsymbol \theta ij},
  \boldsymbol S\_{\boldsymbol \theta ij}\bigr)\$\$

  when applied to the treatment-component parameters, and

  \$\$\boldsymbol \nu\_{kj} \sim p\_{\boldsymbol \nu\_{kj}}
  \\\text{N}\bigl(\mu\_{\boldsymbol \nu}, \boldsymbol
  \Sigma\_{\boldsymbol \nu}\bigr)\_k + (1-p\_{\boldsymbol \nu\_{kj}})\\
  \text{N}(m\_{\boldsymbol \nu\_{kj}}, s^2\_{\boldsymbol \nu\_{kj}})\$\$

  when applied to the interaction parameters. The *exchangeability
  weights* \\p\_{\boldsymbol \theta ij}\\ and \\p\_{\boldsymbol
  \nu\_{kj}}\\ are fixed constants in the interval \\\[0,1\]\\ that
  control the degree to which inference for group \\j\\ is informed by
  the exchangeable mixture components. Larger values for the weights
  correspond to greater exchange of information, while smaller values
  increase robustness in case of outlying observations in individual
  groups \\j\\.

## References

Neuenschwander, B., Roychoudhury, S., & Schmidli, H. (2016). On the use
of co-data in clinical trials. *Statistics in Biopharmaceutical
Research*, 8(3), 345-354.

Neuenschwander, B., Wandel, S., Roychoudhury, S., & Bailey, S. (2016).
Robust exchangeability designs for early phase clinical trials with
multiple strata. *Pharmaceutical statistics*, 15(2), 123-134.

Neuenschwander, B., Branson, M., & Gsponer, T. (2008). Critical aspects
of the Bayesian approach to phase I cancer trials. *Statistics in
medicine*, 27(13), 2420-2439.

Neuenschwander, B., Matano, A., Tang, Z., Roychoudhury, S., Wandel, S.
Bailey, Stuart. (2014). A Bayesian Industry Approach to Phase I
Combination Trials in Oncology. In *Statistical methods in drug
combination studies* (Vol. 69). CRC Press.

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

# fit an example model. See documentation for "combo3" example
example_model("combo3")
#> Running combo3 example:
#>  ## example combo3
#> 
#>  library(abind)
#> 
#>  dref <- c(500, 500, 1000)
#>  num_comp <- 3
#>  num_inter <- choose(3, 2) + 1
#>  num_strata <- nlevels(hist_combo3$stratum_id)
#>  num_groups <- nlevels(hist_combo3$group_id)
#> 
#>  blrmfit <- blrm_exnex(
#>    cbind(num_toxicities, num_patients - num_toxicities) ~
#>      1 + I(log(drug_A / dref[1])) |
#>        1 + I(log(drug_B / dref[2])) |
#>        1 + I(log(drug_C / dref[3])) |
#>        0
#>        + I(drug_A / dref[1] * drug_B / dref[2])
#>          + I(drug_A / dref[1] * drug_C / dref[3])
#>          + I(drug_B / dref[2] * drug_C / dref[3])
#>          + I(drug_A / dref[1] * drug_B / dref[2] * drug_C / dref[3]) |
#>        stratum_id / group_id,
#>    data = hist_combo3,
#>    prior_EX_mu_comp = replicate(num_comp, mixmvnorm(c(1, logit(1/3), 0, diag(c(2^2, 1)))), FALSE),
#>    prior_EX_tau_comp = list(replicate(num_comp,
#>                                       mixmvnorm(c(1, log(c(0.25, 0.125)),
#>                                                 diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE),
#>                             replicate(num_comp,
#>                                       mixmvnorm(c(1, log(2 * c(0.25, 0.125)),
#>                                                 diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE)),
#>    prior_EX_mu_inter = mixmvnorm(c(1, rep.int(0, num_inter),
#>                                       diag((rep.int(sqrt(2) / 2, num_inter))^2))),
#>    prior_EX_tau_inter = replicate(num_strata,
#>                                   mixmvnorm(c(1, rep.int(log(0.25), num_inter),
#>                                               diag((rep.int(log(2) / 1.96, num_inter))^2))), FALSE),
#>    prior_EX_prob_comp = matrix(0.9, nrow = num_groups, ncol = num_comp),
#>    prior_EX_prob_inter = matrix(1.0, nrow = num_groups, ncol = num_inter),
#>    prior_is_EXNEX_comp = rep(TRUE, num_comp),
#>    prior_is_EXNEX_inter = rep(FALSE, num_inter),
#>    prior_tau_dist = 1,
#>    prior_PD = FALSE
#>  )
#> Warning: There were 2 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# print a summary of the prior
prior_summary(blrmfit, digits = 3)
#> Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability
#> 
#> Mixture configuration
#> ---------------------
#> EXNEX components : 3 
#> component
#> I(log(drug_A/dref[1])) I(log(drug_B/dref[2])) I(log(drug_C/dref[3])) 
#>                      1                      1                      1 
#> 
#> EXNEX interactions: 0 
#> interaction
#>                  I(drug_A/dref[1] * drug_B/dref[2]) 
#>                                                   0 
#>                  I(drug_A/dref[1] * drug_C/dref[3]) 
#>                                                   0 
#>                  I(drug_B/dref[2] * drug_C/dref[3]) 
#>                                                   0 
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 
#>                                                   0 
#> 
#> Prior probability for exchangeability per group
#>             component
#> group        I(log(drug_A/dref[1])) I(log(drug_B/dref[2]))
#>   Combo                         0.9                    0.9
#>   HistAgent1                    0.9                    0.9
#>   HistAgent2                    0.9                    0.9
#>             component
#> group        I(log(drug_C/dref[3]))
#>   Combo                         0.9
#>   HistAgent1                    0.9
#>   HistAgent2                    0.9
#> 
#>             interaction
#> group        I(drug_A/dref[1] * drug_B/dref[2])
#>   Combo                                       1
#>   HistAgent1                                  1
#>   HistAgent2                                  1
#>             interaction
#> group        I(drug_A/dref[1] * drug_C/dref[3])
#>   Combo                                       1
#>   HistAgent1                                  1
#>   HistAgent2                                  1
#>             interaction
#> group        I(drug_B/dref[2] * drug_C/dref[3])
#>   Combo                                       1
#>   HistAgent1                                  1
#>   HistAgent2                                  1
#>             interaction
#> group        I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])
#>   Combo                                                        1
#>   HistAgent1                                                   1
#>   HistAgent2                                                   1
#> 
#> EXchangable hyperparameter priors
#> ---------------------------------
#> Component parameters
#> Mean mu_log_beta
#>                               prior weight m_intercept m_log_slope s_intercept s_log_slope    rho
#> component              mix                                                                       
#> I(log(drug_A/dref[1])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_B/dref[2])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_C/dref[3])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> 
#> Heterogeneity tau_log_beta (log-normal)
#>                                         prior weight m_tau_intercept m_tau_log_slope s_tau_intercept s_tau_log_slope    rho
#> stratum   component              mix                                                                                       
#> stratum_1 I(log(drug_A/dref[1])) comp_1        1.000          -1.386          -2.079           0.707           0.707  0.000
#>           I(log(drug_B/dref[2])) comp_1        1.000          -1.386          -2.079           0.707           0.707  0.000
#>           I(log(drug_C/dref[3])) comp_1        1.000          -1.386          -2.079           0.707           0.707  0.000
#> stratum_2 I(log(drug_A/dref[1])) comp_1        1.000          -0.693          -1.386           0.707           0.707  0.000
#>           I(log(drug_B/dref[2])) comp_1        1.000          -0.693          -1.386           0.707           0.707  0.000
#>           I(log(drug_C/dref[3])) comp_1        1.000          -0.693          -1.386           0.707           0.707  0.000
#> 
#> Correlation LKJ
#> component
#> I(log(drug_A/dref[1])) I(log(drug_B/dref[2])) I(log(drug_C/dref[3])) 
#>                      1                      1                      1 
#> 
#> Interaction parameters
#> Mean mu_eta
#>       prior     w  m[1]  m[2]  m[3]  m[4]  s[1]  s[2]  s[3]  s[4] rho[2,1] rho[3,1] rho[4,1] rho[3,2] rho[4,2] rho[4,3]
#> mix                                                                                                                    
#> comp1       1.000 0.000 0.000 0.000 0.000 0.707 0.707 0.707 0.707    0.000    0.000    0.000    0.000    0.000    0.000
#> 
#> Heterogeneity tau_eta (log-normal)
#>                 prior      w   m[1]   m[2]   m[3]   m[4]   s[1]   s[2]   s[3]   s[4] rho[2,1] rho[3,1] rho[4,1] rho[3,2] rho[4,2] rho[4,3]
#> stratum   mix                                                                                                                             
#> stratum_1 comp1        1.000 -1.386 -1.386 -1.386 -1.386  0.354  0.354  0.354  0.354    0.000    0.000    0.000    0.000    0.000    0.000
#> stratum_2 comp1        1.000 -1.386 -1.386 -1.386 -1.386  0.354  0.354  0.354  0.354    0.000    0.000    0.000    0.000    0.000    0.000
#> 
#> Correlation LKJ
#> interaction 
#>           1 
#> 
#> NonEXchangable priors
#> ---------------------
#> Component parameters
#> Mean mu_log_beta
#>                               prior weight m_intercept m_log_slope s_intercept s_log_slope    rho
#> component              mix                                                                       
#> I(log(drug_A/dref[1])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_B/dref[2])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> I(log(drug_C/dref[3])) comp_1        1.000      -0.693       0.000       2.000       1.000  0.000
#> 
#> Interaction parameters
#> Mean mu_eta
#>       prior     w  m[1]  m[2]  m[3]  m[4]  s[1]  s[2]  s[3]  s[4] rho[2,1] rho[3,1] rho[4,1] rho[3,2] rho[4,2] rho[4,3]
#> mix                                                                                                                    
#> comp1       1.000 0.000 0.000 0.000 0.000 0.707 0.707 0.707 0.707    0.000    0.000    0.000    0.000    0.000    0.000

# print a summary of the posterior (model parameters)
print(blrmfit)
#> Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability
#> 
#> Number of observations: 18 
#> Number of groups      : 3 
#> Number of strata      : 2 
#> Number of components  : 3 
#> Number of interactions: 4 
#> EXNEX components      : 3 
#> EXNEX interactions    : 0 
#> 
#> Observations per group:
#>        Group n Stratum n_total
#> 1      Combo 3     BID      17
#> 2 HistAgent1 7      QD      52
#> 3 HistAgent2 8      QD      32
#> 
#> Groups per stratum:
#>   Stratum Groups n_total
#> 1     BID      1      17
#> 2      QD      2      84
#> 
#> Component posterior:
#> Population mean posterior mu_log_beta
#> intercept:
#>                         mean se_mean  sd 2.5%   50% 97.5% n_eff Rhat
#> I(log(drug_A/dref[1])) -1.44    0.37 1.2 -3.2 -1.06  0.22    10 0.90
#> I(log(drug_B/dref[2])) -3.21    0.25 0.8 -4.5 -3.01 -2.24    10 0.93
#> I(log(drug_C/dref[3])) -0.53    0.22 0.7 -1.5 -0.31  0.45    10 0.90
#> log-slope:
#>                        mean se_mean   sd  2.5%  50% 97.5% n_eff Rhat
#> I(log(drug_A/dref[1])) 0.13    0.28 0.89 -1.42 0.12  1.30    10 0.91
#> I(log(drug_B/dref[2])) 0.13    0.19 0.60 -0.75 0.42  0.89    10 1.03
#> I(log(drug_C/dref[3])) 0.40    0.33 1.05 -1.39 0.61  1.68    10 0.90
#> 
#> Population heterogeniety posterior tau_log_beta
#> intercept:
#>                            mean se_mean   sd  2.5%  50% 97.5% n_eff Rhat
#> BID,I(log(drug_A/dref[1])) 0.30   0.039 0.12 0.151 0.29  0.48    10 0.90
#> BID,I(log(drug_B/dref[2])) 0.58   0.184 0.58 0.043 0.26  1.49    10 1.08
#> BID,I(log(drug_C/dref[3])) 0.44   0.110 0.35 0.072 0.33  1.06    10 1.04
#> QD,I(log(drug_A/dref[1]))  0.48   0.086 0.27 0.190 0.32  0.91    10 0.89
#> QD,I(log(drug_B/dref[2]))  0.50   0.126 0.40 0.164 0.26  1.14    10 0.96
#> QD,I(log(drug_C/dref[3]))  0.74   0.181 0.57 0.164 0.63  1.78    10 0.92
#> log-slope:
#>                            mean se_mean   sd  2.5%   50% 97.5% n_eff Rhat
#> BID,I(log(drug_A/dref[1])) 0.15   0.036 0.11 0.045 0.126  0.40  10.0 0.94
#> BID,I(log(drug_B/dref[2])) 0.19   0.073 0.23 0.025 0.048  0.62  10.0 0.96
#> BID,I(log(drug_C/dref[3])) 0.18   0.032 0.08 0.092 0.151  0.29   6.4 1.20
#> QD,I(log(drug_A/dref[1]))  0.29   0.052 0.16 0.106 0.277  0.50  10.0 0.95
#> QD,I(log(drug_B/dref[2]))  0.48   0.147 0.46 0.047 0.355  1.24  10.0 1.01
#> QD,I(log(drug_C/dref[3]))  0.28   0.049 0.16 0.117 0.238  0.58  10.0 0.90
#> 
#> Population correlation posterior rho_log_beta
#>                          mean se_mean   sd  2.5%   50% 97.5% n_eff Rhat
#> I(log(drug_A/dref[1])) -0.048    0.21 0.68 -0.88 -0.25  0.92    10  1.1
#> I(log(drug_B/dref[2]))  0.093    0.24 0.75 -0.90  0.23  0.97    10  1.0
#> I(log(drug_C/dref[3]))  0.136    0.15 0.49 -0.69  0.33  0.75    10  0.9
#> 
#> Interaction model posterior:
#> Population mean posterior mu_eta
#>                                                      mean se_mean   sd  2.5%
#> I(drug_A/dref[1] * drug_B/dref[2])                  -0.47    0.22 0.69 -1.36
#> I(drug_A/dref[1] * drug_C/dref[3])                   0.18    0.22 0.70 -0.81
#> I(drug_B/dref[2] * drug_C/dref[3])                   0.11    0.10 0.32 -0.49
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.12    0.14 0.44 -0.48
#>                                                       50% 97.5% n_eff Rhat
#> I(drug_A/dref[1] * drug_B/dref[2])                  -0.61  0.73    10 0.91
#> I(drug_A/dref[1] * drug_C/dref[3])                   0.23  1.08    10 1.01
#> I(drug_B/dref[2] * drug_C/dref[3])                   0.17  0.46    10 0.92
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.11  0.72    10 0.93
#> 
#> Population heterogeniety posterior tau_eta
#>                                                         mean se_mean    sd 2.5%
#> BID,I(drug_A/dref[1] * drug_B/dref[2])                  0.25   0.021 0.061 0.16
#> BID,I(drug_A/dref[1] * drug_C/dref[3])                  0.23   0.026 0.082 0.14
#> BID,I(drug_B/dref[2] * drug_C/dref[3])                  0.35   0.043 0.136 0.18
#> BID,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 0.22   0.017 0.055 0.16
#> QD,I(drug_A/dref[1] * drug_B/dref[2])                   0.25   0.025 0.080 0.16
#> QD,I(drug_A/dref[1] * drug_C/dref[3])                   0.23   0.053 0.168 0.11
#> QD,I(drug_B/dref[2] * drug_C/dref[3])                   0.28   0.020 0.063 0.20
#> QD,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.21   0.040 0.103 0.15
#>                                                          50% 97.5% n_eff Rhat
#> BID,I(drug_A/dref[1] * drug_B/dref[2])                  0.24  0.34   8.1 1.12
#> BID,I(drug_A/dref[1] * drug_C/dref[3])                  0.21  0.37  10.0 0.94
#> BID,I(drug_B/dref[2] * drug_C/dref[3])                  0.34  0.62  10.0 1.25
#> BID,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 0.21  0.33  10.0 1.41
#> QD,I(drug_A/dref[1] * drug_B/dref[2])                   0.25  0.39  10.0 0.90
#> QD,I(drug_A/dref[1] * drug_C/dref[3])                   0.16  0.60  10.0 0.98
#> QD,I(drug_B/dref[2] * drug_C/dref[3])                   0.26  0.37  10.0 0.93
#> QD,I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  0.17  0.42   6.6 1.34
#> 
#> Population correlation posterior Sigma_corr_eta
#>                                                                                                           mean
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.000
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.072
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.044
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.031
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    0.072
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.000
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.110
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.045
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.044
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.110
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.000
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.037
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  -0.031
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                   0.045
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  -0.037
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.000
#>                                                                                                         se_mean
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                       NaN
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   1.1e-01
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   5.1e-02
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  5.3e-02
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   1.1e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   1.2e-17
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   1.7e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  8.5e-02
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   5.1e-02
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   1.7e-01
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   2.3e-17
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  2.0e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  5.3e-02
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  8.5e-02
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  2.0e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 3.4e-17
#>                                                                                                              sd
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                   0.0e+00
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   3.5e-01
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   1.6e-01
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  1.7e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   3.5e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   3.7e-17
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   5.5e-01
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  2.7e-01
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   1.6e-01
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   5.5e-01
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   7.4e-17
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  6.4e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  1.7e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  2.7e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  6.4e-01
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 8.3e-17
#>                                                                                                          2.5%
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.00
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.47
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.34
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.28
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.47
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.00
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.86
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.33
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   -0.34
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.86
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.00
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.80
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  -0.28
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  -0.33
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  -0.80
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.00
#>                                                                                                            50%
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.000
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.106
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                    0.019
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.012
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    0.106
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.000
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   -0.198
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.026
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    0.019
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   -0.198
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.000
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  -0.280
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                   0.012
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                   0.026
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  -0.280
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.000
#>                                                                                                         97.5%
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    1.00
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.61
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                    0.11
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.17
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    0.61
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    1.00
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    0.76
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.36
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                    0.11
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                    0.76
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                    1.00
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                   0.91
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                   0.17
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                   0.36
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                   0.91
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])  1.00
#>                                                                                                         n_eff
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                     NaN
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                      10
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                      10
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                     10
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                      10
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                      10
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                      10
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                     10
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                      10
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                      10
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                      10
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                     10
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                     10
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                     10
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                     10
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])     6
#>                                                                                                         Rhat
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2])                                    NaN
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_C/dref[3])                                   0.93
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_B/dref[2] * drug_C/dref[3])                                   0.90
#> I(drug_A/dref[1] * drug_B/dref[2]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  0.91
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   0.93
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   0.89
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   0.90
#> I(drug_A/dref[1] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  1.03
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                                   0.90
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                                   0.90
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                                   0.89
#> I(drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3])                  0.90
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2])                  0.91
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_C/dref[3])                  1.03
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_B/dref[2] * drug_C/dref[3])                  0.90
#> I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]),I(drug_A/dref[1] * drug_B/dref[2] * drug_C/dref[3]) 0.89
#> Warning: Parts of the model have not converged (some Rhats are > 1.1).
#> Be careful when analysing the results! It is recommend to run
#> more iterations and/or setting stronger priors.

# summary of posterior for DLT rate by dose for observed covariate levels
summ <- summary(blrmfit, interval_prob = c(0, 0.16, 0.33, 1))
print(cbind(hist_combo3, summ))
#>    stratum_id   group_id drug_A drug_B drug_C num_toxicities num_patients
#> 1         BID      Combo    400    800     80              0            4
#> 2         BID      Combo    400    800    160              1            8
#> 3         BID      Combo    400    800    240              2            5
#> 4          QD HistAgent1      0    240      0              0            5
#> 5          QD HistAgent1      0    400      0              0           10
#> 6          QD HistAgent1      0    800      0              1           11
#> 7          QD HistAgent1    240     80      0              0            4
#> 8          QD HistAgent1    400    400      0              1            6
#> 9          QD HistAgent1    400    800      0              0           11
#> 10         QD HistAgent1    400   1000      0              1            5
#> 11         QD HistAgent2      0      0     80              0            3
#> 12         QD HistAgent2      0      0    160              0            3
#> 13         QD HistAgent2      0      0    320              0            6
#> 14         QD HistAgent2      0      0    480              0            3
#> 15         QD HistAgent2      0      0    640              1            5
#> 16         QD HistAgent2    400      0    160              1            6
#> 17         QD HistAgent2    400      0    320              2            5
#> 18         QD HistAgent2    400      0    240              1            1
#>          mean         sd         2.5%         50%      97.5% [0,0.16]
#> 1  0.20371098 0.07733550 9.620129e-02 0.201822163 0.31968567      0.3
#> 2  0.22361330 0.08007056 1.098705e-01 0.227943497 0.33955953      0.2
#> 3  0.24476545 0.08390434 1.264024e-01 0.253043507 0.36250684      0.2
#> 4  0.01703650 0.01085987 2.704102e-03 0.015709715 0.03293509      1.0
#> 5  0.03253753 0.02215977 4.940107e-03 0.029305529 0.07052722      1.0
#> 6  0.08482887 0.06988129 9.477103e-03 0.080687466 0.19800189      0.8
#> 7  0.08146756 0.09038718 1.074730e-02 0.039476231 0.26672691      0.8
#> 8  0.13560027 0.07384235 5.530077e-02 0.109265821 0.25538344      0.6
#> 9  0.12164949 0.05012800 6.007321e-02 0.110126570 0.21385390      0.8
#> 10 0.11992114 0.06138381 5.630877e-02 0.094002190 0.23095714      0.7
#> 11 0.02545156 0.04149528 2.894899e-05 0.004400239 0.11155608      1.0
#> 12 0.04024950 0.06449306 3.633123e-04 0.012684074 0.17771327      0.9
#> 13 0.06864600 0.09319206 3.864447e-03 0.034094033 0.27053479      0.9
#> 14 0.10557835 0.10764986 1.051573e-02 0.075278930 0.33633626      0.8
#> 15 0.15816416 0.11867879 3.671271e-02 0.128517617 0.39144890      0.6
#> 16 0.30206204 0.15060725 1.601939e-01 0.256111022 0.59069877      0.1
#> 17 0.32930774 0.13994434 1.998291e-01 0.289977332 0.59712922      0.0
#> 18 0.31525645 0.14517816 1.756684e-01 0.264551860 0.59125717      0.0
#>    (0.16,0.33] (0.33,1]
#> 1          0.7      0.0
#> 2          0.7      0.1
#> 3          0.7      0.1
#> 4          0.0      0.0
#> 5          0.0      0.0
#> 6          0.2      0.0
#> 7          0.2      0.0
#> 8          0.4      0.0
#> 9          0.2      0.0
#> 10         0.3      0.0
#> 11         0.0      0.0
#> 12         0.1      0.0
#> 13         0.1      0.0
#> 14         0.1      0.1
#> 15         0.3      0.1
#> 16         0.6      0.3
#> 17         0.7      0.3
#> 18         0.7      0.3

# summary of posterior for DLT rate by dose for new set of covariate levels
newdata <- expand.grid(
  stratum_id = "BID", group_id = "Combo",
  drug_A = 400, drug_B = 800, drug_C = c(320, 400, 600, 800),
  stringsAsFactors = FALSE
)
summ_pred <- summary(blrmfit, newdata = newdata, interval_prob = c(0, 0.16, 0.33, 1))
print(cbind(newdata, summ_pred))
#>   stratum_id group_id drug_A drug_B drug_C      mean         sd      2.5%
#> 1        BID    Combo    400    800    320 0.2680163 0.08890877 0.1474356
#> 2        BID    Combo    400    800    400 0.2944229 0.09540909 0.1669803
#> 3        BID    Combo    400    800    600 0.3784127 0.12889709 0.2023964
#> 4        BID    Combo    400    800    800 0.4778682 0.18929560 0.2046738
#>         50%     97.5% [0,0.16] (0.16,0.33] (0.33,1]
#> 1 0.2663310 0.3939148      0.1         0.6      0.3
#> 2 0.2912059 0.4390051      0.0         0.6      0.4
#> 3 0.3844232 0.5885953      0.0         0.3      0.7
#> 4 0.4780968 0.7661702      0.0         0.2      0.8

# update the model after observing additional data
newdata$num_patients <- rep(3, nrow(newdata))
newdata$num_toxicities <- c(0, 1, 2, 2)
library(dplyr)
#> 
#> Attaching package: ‘dplyr’
#> The following objects are masked from ‘package:stats’:
#> 
#>     filter, lag
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, setdiff, setequal, union
blrmfit_new <- update(blrmfit,
  data = rbind(hist_combo3, newdata) %>%
    arrange(stratum_id, group_id)
)
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

# updated posterior summary
summ_upd <- summary(blrmfit_new, newdata = newdata, interval_prob = c(0, 0.16, 0.33, 1))
print(cbind(newdata, summ_upd))
#>   stratum_id group_id drug_A drug_B drug_C num_patients num_toxicities
#> 1        BID    Combo    400    800    320            3              0
#> 2        BID    Combo    400    800    400            3              1
#> 3        BID    Combo    400    800    600            3              2
#> 4        BID    Combo    400    800    800            3              2
#>        mean         sd      2.5%       50%     97.5% [0,0.16] (0.16,0.33]
#> 1 0.2866255 0.08415192 0.1530756 0.3084664 0.3915255      0.2         0.5
#> 2 0.3401901 0.09704336 0.1798541 0.3576501 0.4586301      0.0         0.4
#> 3 0.4992189 0.13063222 0.2828894 0.5138069 0.6708180      0.0         0.1
#> 4 0.6575665 0.14815459 0.4446452 0.6602006 0.8403893      0.0         0.0
#>   (0.33,1]
#> 1      0.3
#> 2      0.6
#> 3      0.9
#> 4      1.0
## Recover user set sampling defaults
options(.user_mc_options)
```
