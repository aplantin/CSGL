#' Calculates maximum lambda value
#'
#' Calculates maximum lambda value as smallest lambda for which
#' all coefficients are zero.
#'
#' @param Z Matrix of log OTU proportions
#' @param y Outcome vector
#' @param n Number of subjects
#' @param p Number of OTUs
#' @param H Matrix I_p - (1/p) 1_p 1_p'
#' @param mu Augmented Lagrangian parameter
#' @param theta Defines convex combination between L1 and L2 penalties
#'     (theta = 1 is L1 only, theta = 0 is L2 only)
#' @param groups Vector indicating group membership of each OTU
#' @param thresh Threshold for convergence
#' @param erel Internal parameter for ADMM convergence
#' @param maxit Maximum number of iterations
#'
#' @importFrom stats sd
#'
#' @return Maximum lambda value
calc.lammax <- function(Z, y, n, p, H, mu, theta, groups, thresh, erel, maxit) {
  current.max <- max(t(Z) %*% y)
  beta <- rep(0, p)
  allzero <- all(beta == 0)
  while (allzero) {
    beta <- csgl.fit(Z = Z, y = y, n = n, p = p, H = H, mu = mu,
                     lambda = current.max, theta = theta,
                     groups = groups, thresh = thresh,
                     erel = erel, maxit = maxit)
    allzero = all(beta == 0)
    current.max = current.max*0.9
  }
  lam.max = current.max/0.8
  return(lam.max)
}

#' Compositional sparse group lasso
#'
#' Fits CSGL for a sequence of lambda values.
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
#' @return List with components
#' \item{int}{Vector of intercepts, one per lambda}
#' \item{beta}{Matrix of coefficients, one column per lambda}
#' \item{lambdas}{Sequence of lambda values}
#'
#' @export
csgl <- function(Z, y, groups, mu = 1, theta = 0.95, lam.seq = NULL, nlam = 25,
                 min.frac = 0.1, thresh = 1e-7, maxit = 1e4, std = T, verbose = F) {
  n = nrow(Z)
  p = ncol(Z)
  if (std) {
    y <- scale(y, scale = FALSE)
    Z <- scale(Z, scale = apply(Z, 2, sd)*sqrt(nrow(Z) - 1))
    fac <- 1/attr(Z, "scaled:scale")
  } else {
    fac <- rep(1, p)
  }
  erel = 1e-4
  H <- diag(p) - matrix(1/p, nrow = p, ncol = p)

  ## get max lam
  if (is.null(lam.seq)) {
    lam.max <- calc.lammax(Z = Z, y = y, n = n, p = p, H = H, mu = mu,
                           theta = theta, groups = groups, thresh  = thresh,
                           erel = erel, maxit  = maxit)
    lam.seq <- exp(seq(from = log(lam.max), to = log(lam.max*min.frac),
                       length.out = nlam))
  } else {
    nlam = length(lam.seq)
  }

  ## main fit
  betas <- matrix(nrow = p, ncol = nlam)
  for (i in 1:nlam) {
    if (verbose) {
      print(paste("This is lambda", i, "with value", lam.seq[i]))
    }
    betas[,i] <- csgl.fit(Z = Z, y  = y, n =  n, p = p, H = H, mu = mu,
                          lambda = lam.seq[i], theta = theta, groups = groups,
                          thresh = thresh, erel = erel, maxit = maxit)
  }

  ## rescale if standardized
  if (std) {
    ints <- c()
    for (i in 1:nlam) {
      betas[,i] <- betas[,i] * fac
      ints[i] <- attr(y, "scaled:center") -
        drop(crossprod(betas[,i], attr(Z, "scaled:center")))
    }
  } else {
    ints = NA
  }
  return(list(int = ints, beta = betas, lambdas = lam.seq))
}