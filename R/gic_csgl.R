#' Calculates generalized information criterion
#'
#' Calculates GIC as defined in Lin (2014).
#'
#' @param Z Matrix of log OTU proportions
#' @param y Outcome vector
#' @param betahat Estimated coefficient vector
#' @param int Estimated intercept; default = 0
#'
#' @return GIC
calc.gic <- function(Z, y, betahat, int = 0) {
  n <- nrow(Z)
  p <- ncol(Z)
  s.lam <- sum(betahat != 0)
  sig2.lam <- sum( (y - int - Z %*% betahat)^2 )/n
  gic <- log(sig2.lam) + (s.lam-1) * (log(log(n))/n)*log(max(p,n))
  return(gic)
}

#' Fits CSGL, choosing lambda by GIC
#'
#' Fits CSGL for a sequence of lambdas and chooses lambda based on the GIC.
#'
#' @param Z Matrix of log OTU proportions
#' @param y Outcome vector
#' @param groups Vector indicating group membership of each OTU
#' @param mu Augmented Lagrangian parameter
#' @param theta Defines convex combination between L1 and L2 penalties
#'     (theta = 1 is L1 only, theta = 0 is L2 only)
#' @param lam.seq Sequence of lambda values to consider
#' @param nlam Number of lambda values to try
#' @param min.frac Minimum value of the penalty parameter, as a fraction of the maximum value
#' @param thresh Threshold for convergence
#' @param maxit Maximum number of iterations
#' @param std Logical flag for variable standardization before fitting the model
#' @param verbose Logical flag for printing updates during model fitting
#'
#' @return A list with components
#' \item{gic.int}{Intercept for lambda chosen by GIC}
#' \item{gic.beta}{Coefficient vector for lambda chosen by GIC}
#' \item{fullfit}{Full model fit (all lambda values)}
#' \item{gic}{Vector of GIC values for each lambda}
#'
#' @export
gic.csgl <- function(Z, y, groups, mu = 1, theta = 0.95, lam.seq = NULL, nlam = 25,
                     min.frac = 0.1, thresh = 1e-7, maxit = 1e4, std = T, verbose = F) {
  if (!is.null(lam.seq)) {
    nlam <- length(lam.seq)
  }
  mod <- csgl(Z = Z, y = y, groups = groups, mu = mu, theta = theta,
              lam.seq = lam.seq, nlam = nlam, min.frac = min.frac,
              thresh = thresh, maxit = maxit, std = std, verbose = verbose)
  gics <- sapply(1:nlam, FUN = function(i) {
    calc.gic(Z, y, mod$beta[,i], mod$int[i])
  })
  gic.beta <- mod$beta[, which.min(gics)]
  gic.int <- mod$int[which.min(gics)]
  return(list(gic.int = gic.int, gic.beta = gic.beta, fullfit = mod, gic = gics))
}

