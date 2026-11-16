##################################################################################
##  Maximum-likelihood fit of the permeability model of Barbour et al. (2026)   ##
##  with the BFGS optimiser and the analytic gradient of the log-likelihood.    ##
##################################################################################
#
#  The function fit_perm() fits the same model as permeability::fitPermeability():
#
#  fitPermeability() uses optim() with
#   its defaults (Nelder-Mead, 500 iterations); fit_perm() uses BFGS with the
#   analytic gradient (relative tolerance 1e-14, up to 5,000 iterations) and
#   reports whether the fit converged.
#

# Log-likelihood (kappaLogLikelihood of the package) and its analytic gradient
#---------------------------------------------------------------------------------
# X        standardised design matrix (n steps x k parameters)
# crossed  logical, TRUE if the observed step crossed
# p_stay   proportion of null steps that did not cross
# Returns two functions of the coefficients b: fn (log-likelihood), gr (gradient)
perm_ll_grad <- function(X, crossed, p_stay) {
  lp <- log(p_stay)
  list(
    fn = function(b) {
      k <- as.vector(exp(X %*% b)); u <- k * lp            # u = log(p_stay ^ kappa)
      sum(log(-expm1(u[crossed]))) + sum(u[!crossed])
    },
    gr = function(b) {
      k <- as.vector(exp(X %*% b)); u <- k * lp
      w <- ifelse(crossed, -exp(u) * u / (-expm1(u)), u)   # d logLik / d (X %*% b)
      as.vector(crossprod(X, w))
    })
}


# Fit the permeability model with BFGS
#---------------------------------------------------------------------------------
# formula  one-sided formula for log(kappa), e.g. ~ country * traffic_class
# data     crossing table (see the end of this file)
# start    optional starting values (standardised scale, one per parameter)
# maxit, reltol  passed to optim()
# grad_tol       largest absolute gradient accepted for a converged fit
fit_perm <- function(formula, data, start = NULL,
                     maxit = 5000, reltol = 1e-14, grad_tol = 1e-3) {
  if (!requireNamespace("permeability", quietly = TRUE))
    stop("The 'permeability' package (Barbour et al. 2026) is needed.")
  data <- as.data.frame(data)                 # the package needs a data.frame, not a tibble
  need <- c("crossed", "null.crossed", "null.stayed")
  if (!all(need %in% names(data)))
    stop("Columns missing in data: ", paste(setdiff(need, names(data)), collapse = ", "))

  # Design matrix, standardised as in fitPermeability()
  X.raw <- model.matrix(formula, data = data)
  if (nrow(X.raw) != nrow(data))
    stop("Missing values in the covariates of the formula: remove those rows first.")
  X <- permeability:::scaleX(X.raw)

  crossed <- data$crossed == 1
  p_stay  <- data$null.stayed / (data$null.crossed + data$null.stayed)
  if (any(p_stay[!crossed] == 0) || any(p_stay[crossed] == 1))
    stop("The log-likelihood is -Inf for ", sum(p_stay[!crossed] == 0),
         " non-crossing step(s) with null.stayed = 0 and ", sum(p_stay[crossed] == 1),
         " crossing step(s) with null.crossed = 0. Handle these steps before fitting.")

  f <- perm_ll_grad(unclass(X)[, , drop = FALSE], crossed, p_stay)

  # Starting values: 0 for the covariates and a rough estimate for the intercept
  b0 <- setNames(rep(0, ncol(X)), colnames(X))
  b0[1] <- log(-log(1 - mean(crossed)) / mean(-log(p_stay)))
  if (!is.finite(b0[1])) b0[1] <- 0
  if (!is.null(start)) {
    if (length(start) != length(b0)) stop("`start` needs ", length(b0), " values.")
    b0[] <- start
  }

  # Maximise the log-likelihood (fnscale = -1)
  o <- optim(b0, f$fn, f$gr, method = "BFGS",
             control = list(fnscale = -1, maxit = maxit, reltol = reltol))

  # Standard errors from the Hessian at the maximum
  H  <- optimHess(o$par, f$fn, f$gr, control = list(fnscale = -1))
  se <- tryCatch(sqrt(diag(solve(-H))), error = function(e) rep(NA_real_, length(o$par)))

  coef <- data.frame(estimate = o$par, se = se)
  coef$z_score <- coef$estimate / coef$se
  coef$p_value <- 2 * (1 - pnorm(abs(coef$z_score)))
  o$hessian <- H

  structure(list(coef.table   = coef,
                 logLik       = o$value,
                 AIC          = -2 * o$value + 2 * nrow(coef),
                 formula      = formula,
                 optim.output = o,
                 data         = data,
                 converged    = o$convergence == 0 && max(abs(f$gr(o$par))) < grad_tol),
            class = "permfit")
}

##################################################################################
# TO FIT THE MODELS WITH THIS FUNCTION YOU NEED:
#
#   - R and the 'permeability' package (Barbour et al. 2026).
#   - A crossing table built with the package (prepPermeability() and then
#     buildCrossingTable()), with one row per step and the columns
#       crossed       1 if the observed step crossed the barrier, 0 if not
#       null.crossed  number of null steps that crossed
#       null.stayed   number of null steps that did not cross
#     plus the covariates used in the formula.
#
# The output has the same elements as the output of fitPermeability()
# (coef.table, logLik, AIC, formula, optim.output, data) plus `converged`,
# so the log-likelihood and AIC of the two functions can be compared directly
# for the same formula and data.
#
##################################################################################

##################################################################################
# EXAMPLE
##################################################################################
#
# fit_perm() needs two arguments:
#   formula  a one-sided formula with the covariates of kappa
#   data     the crossing table, one row per step, with the columns
#            crossed, null.crossed, null.stayed and the covariates of the formula

  library(permeability)

  # 1. Crossing table ------------------------------------------------------------
  # You could run this With real data, data provide by Barbour et al. 2026 with package:
  #   permdata       <- prepPermeability(move.df = ..., barrier = ..., ...)
  #   crossing_table <- buildCrossingTable(permdata = permdata, n.null = 60)
  # with the covariates of each step (or of the barrier segment of each step)
  # joined to it as extra columns.
  #
  # Here, this INVENTED table shows the columns that fit_perm() uses two countries data pop1, pop2
  # (2,000 steps, 60 null steps per step):
  set.seed(1)
  n <- 2000
  crossing_table <- data.frame(
    ID            = sample(paste0("individual", 1:10), n, replace = TRUE),   # not used by fit_perm()
    country       = factor(sample(c("pop1", "pop2"), n, replace = TRUE), levels = c("pop1", "pop2")),
    traffic_class = factor(sample(c("A", "B", "C"), n, replace = TRUE), levels = c("A", "B", "C")),
    log_aadt      = runif(n, 3, 4.5),               # continuous covariate
    null.crossed  = rbinom(n, 60, 0.1))             # null steps that crossed
  crossing_table$null.stayed <- 60 - crossing_table$null.crossed       # null steps that did not cross
  # observed crossings, simulated with kappa = 0.1 (pop1) and 0.4 (pop2)
  p_stay <- crossing_table$null.stayed / 60
  kappa  <- ifelse(crossing_table$country == "pop2", 0.4, 0.1)
  crossing_table$crossed <- rbinom(n, 1, 1 - p_stay^kappa)             # 1 = crossed, 0 = not
  head(crossing_table)

  # 2. Fit -----------------------------------------------------------------------
  fit0 <- fit_perm(~ 1,                       data = crossing_table)   # one kappa for all steps
  fit1 <- fit_perm(~ country,                 data = crossing_table)   # categorical covariate
  fit2 <- fit_perm(~ country + log_aadt,      data = crossing_table)   # + continuous covariate
  fit3 <- fit_perm(~ country * traffic_class, data = crossing_table)   # interaction

  # 3. Output --------------------------------------------------------------------
  fit1$coef.table    # estimate, se, z_score, p_value (standardised scale of scaleX)
  fit1$logLik
  fit1$AIC
  fit1$converged     # must be TRUE

  # Model comparison
  data.frame(model     = c("~ 1", "~ country", "~ country + log_aadt", "~ country * traffic_class"),
             AIC       = c(fit0$AIC, fit1$AIC, fit2$AIC, fit3$AIC),
             converged = c(fit0$converged, fit1$converged, fit2$converged, fit3$converged))

