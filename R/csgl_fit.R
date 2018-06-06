#' Fits CSGL
#'
#' Function to fit CSGL. Interior function for csgl (sequence of lambdas),
#' gic.csgl (choosing best lambda by GIC), and cv.csgl (choosing best
#' lambda by k-fold cross-validation).
#'
#' @param Z Matrix of log OTU proportions
#' @param y Outcome vector
#' @param n Number of subjects
#' @param p Number of OTUs
#' @param H Matrix I_p - (1/p) 1_p 1_p'
#' @param mu Augmented Lagrangian parameter
#' @param lambda Regularization parameter
#' @param theta Defines convex combination between L1 and L2 penalties
#'     (theta = 1 is L1 only, theta = 0 is L2 only)
#' @param groups Vector indicating group membership of each OTU
#' @param thresh Threshold for convergence
#' @param erel Internal parameter for ADMM convergence
#' @param maxit Maximum number of iterations
#'
#' @return Fitted coefficients beta
#'
csgl.fit <- function(Z, y, n, p, H, mu, lambda, theta, groups, thresh, erel, maxit) {
  alpha = rep(0, p)
  beta = rep(0, p)
  xi = rep(0, p)

  converged = FALSE
  k=0
  while (k < maxit & !converged) {
    ## updates
    beta <- update.beta(y = y, Z = Z, n = n, p = p, mu = mu, alpha = alpha,
                        xi = xi, lambda = lambda, theta = theta, groups = groups)
    alpha0 <- alpha; alpha <- update.alpha(H, beta, xi, p)
    xi <- update.xi(alpha, beta, xi)

    ## check convergence
    k <- k + 1
    epri <- sqrt(p) * thresh + erel * max(sqrt(sum(beta^2)), sqrt(sum(alpha^2)))
    edual <- sqrt(n) * thresh + erel * mu*sqrt(sum(xi^2))
    sk <- sqrt(sum((mu*(alpha0 - alpha))^2))
    rk <- sqrt(sum((beta - alpha)^2))
    converged = (sk <= edual & rk <= epri)
  }
  return(beta)
}
