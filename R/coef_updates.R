#' Update beta vector
#'
#' This function uses the SGL package to minimize the augmented Lagrangian
#' with respect to beta.
#'
#' @param y Outcome vector
#' @param Z Matrix of log OTU proportions
#' @param n Number of subjects (nrow(Z))
#' @param p Number of OTUs (ncol(Z))
#' @param mu Augmented Lagrangian parameter
#' @param alpha p-vector of coefficients
#' @param xi p-vector of Lagrange multipliers
#' @param lambda Regularization parameter
#' @param theta Defines convex combination between L1 and L2 penalties (theta = 1 is L1 only)
#' @param groups Vector indicating group membership of each OTU
#'
#' @importFrom SGL SGL
#'
#' @return Updated beta vector
update.beta <- function(y, Z, n, p, mu, alpha, xi, lambda, theta, groups) {
  y.tilde <- c(y, sqrt(n*mu) * (xi - alpha))
  Z.tilde <- rbind(Z, -sqrt(n*mu)*diag(p))
  beta.new <- SGL(data = list(x = Z.tilde, y = y.tilde),
                  index = factor(groups), type = "linear", alpha = theta,
                  lambdas = lambda, standardize = FALSE)$beta
  return(beta.new)
}

#' Update alpha vector
#'
#' This function updates alpha by projecting (beta + xi) onto the subspace
#' defined by sum(alpha) = 0.
#'
#' @param H Matrix I_p - (1/p) 1_p 1_p'
#' @param beta p-vector of coefficients
#' @param xi p-vector of Lagrange multipliers
#' @param p Number of OTUs
#'
#' @return Updated alpha vector
update.alpha <- function(H, beta, xi, p) {
  alpha <- H %*% (beta + xi)
  return(alpha)
}

#' Update xi vector
#'
#' This function updates the Lagrange multipliers xi.
#'
#' @param alpha p-vector of coefficients
#' @param beta p-vector of coefficients
#' @param xi p-vector of Lagrange multipliers
#'
#' @return Updated xi vector
update.xi <- function(alpha, beta, xi) {
  xi.new = xi + beta - alpha
  return(xi.new)
}