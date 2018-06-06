#' Gets folds for cross-validation
#'
#' @param n Number of subjects
#' @param nfolds Desired number of folds
#'
#' @return List of folds
#'
getFolds <- function(n, nfolds) {
  reord <- sample(1:n)
  folds <- list()
  perFold = floor(n/nfolds)
  current = 0
  for (i in 1:nfolds) {
    folds[[i]] <- reord[(current + 1:perFold)]
    current = current + perFold
  }
  if (n %% nfolds != 0) {
    for (i in 1:(n %% nfolds)) {
      folds[[i]] <- c(folds[[i]], reord[n - i + 1])
    }
  }
  return(folds)
}

#' Fits CSGL, choosing lambda by cross-validation
#'
#' Fits CSGL for a sequence of lambda values and chooses lambda based on k-fold cross-validation.
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
#' @param nfolds Number of folds for cross-validation
#' @param thresh Threshold for convergence
#' @param maxit Maximum number of iterations
#' @param std Logical flag for variable standardization before fitting the model
#' @param verbose Logical flag for printing updates during model fitting
#'
#' @return A list with components
#' \item{cv.int}{Intercept for lambda chosen by CV}
#' \item{cv.beta}{Coefficient vector for lambda chosen by CV}
#' \item{cv.lam}{Lambda chosen by CV}
#' \item{lamseq}{Sequence of lambda values tried}
#' \item{pe}{Cross-validated prediction error for each lambda}
#'
#' @export
cv.csgl <- function(Z, y, groups, mu = 1, theta = 0.95, lam.seq = NULL,
                    nlam = 25, min.frac = 0.1, nfolds = 10, thresh = 1e-7,
                    maxit = 1e4, std = T, verbose = F) {
  n = nrow(Z)
  p = ncol(Z)

  if (std) {
    y.std <- scale(y, scale = FALSE)
    Z.std <- scale(Z, scale = apply(Z, 2, sd)*sqrt(nrow(Z) - 1))
    fac.std <- 1/attr(Z.std, "scaled:scale")
  } else {
    y.std = y
    Z.std = Z
    fac.std <- rep(1, p)
  }

  erel = 1e-4
  H <- diag(p) - matrix(1/p, nrow = p, ncol = p)

  ## get max lam and sequence of lambdas
  if (is.null(lam.seq)) {
    lam.max <- calc.lammax(Z = Z.std, y = y.std, n = n, p = p, H  = H, mu = mu,
                           theta = theta, groups = groups, thresh = thresh,
                           erel = erel, maxit = maxit)
    lam.seq <- exp(seq(from = log(lam.max), to = log(lam.max*min.frac),
                       length.out = nlam))
  } else {
    nlam <- length(lam.seq)
  }

  ## cross-validated PE
  folds <- getFolds(n, nfolds)
  pe <- matrix(nrow = n, ncol = nlam)

  ## fit for each fold
  for (i in 1:nfolds) {
    print(paste("This is fold", i))

    Z.test <- Z[folds[[i]], ]
    y.test <- y[folds[[i]]]
    Z.train <- Z[-folds[[i]], ]
    y.train <- y[-folds[[i]]]

    if (std) {
      y.train <- scale(y.train, scale = FALSE)
      Z.train <- scale(Z.train, scale = apply(Z.train, 2, sd)*sqrt(nrow(Z.train) - 1))
      fac <- 1/attr(Z.train, "scaled:scale")
    } else {
      fac <- rep(1, p)
    }

    betas <- matrix(nrow = p, ncol = nlam)
    for (j in 1:nlam) {
      betas[,j] <- csgl.fit(Z = Z.train, y = y.train, n = length(y.train),
                            p = p, H = H, mu = mu, lambda = lam.seq[j],
                            theta = theta, groups = groups, thresh = thresh,
                            erel = erel, maxit = maxit)
    }

    ## rescale if standardized, and calculate y - yhat
    if (std) {
      ints <- c()
      for (j in 1:nlam) {
        betas[,j] <- betas[,j] * fac
        ints[j] <- attr(y.train, "scaled:center") -
          drop(crossprod(betas[,j], attr(Z.train, "scaled:center")))
        pe[folds[[i]], j] <- y.test - ints[j] - Z.test %*% betas[,j]
      }
    } else {
      ints = NA
      pe[folds[[i]], j] <- y.test - Z.test %*% betas[,j]
    }
  }
  pe.res <- apply(pe, 2, FUN = function(x) sum(x^2)/n)
  cv.lam <- lam.seq[which.min(pe.res)]
  cv.beta <- csgl.fit(Z = Z.std, y = y.std, n = n, p = p, H = H, mu = mu,
                      lambda = cv.lam, theta = theta, groups = groups,
                      thresh = thresh, erel = erel, maxit = maxit)
  if (std) {
    cv.beta <- cv.beta * fac.std
    cv.int <- attr(y.std, "scaled:center") -
      drop(crossprod(cv.beta, attr(Z.std, "scaled:center")))
  } else {
    cv.int <- NA
  }
  return(list(cv.int = cv.int, cv.beta = cv.beta, cv.lam = cv.lam,
              lamseq = lam.seq, pe = pe.res))
}
