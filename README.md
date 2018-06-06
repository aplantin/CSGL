## CSGL
Fits the Compositional Sparse Group Lasso, an L1- and L2-penalized linear log-contrast model. 

## Installation instructions: 

    library(devtools)
    install_github("aplantin/CSGL")
    library(CSGL) 

## Sample usage with simulated data  

    # Function to generate compositional data 
    library(mvtnorm) 
    gen.comp.data <- function(p, n, tb, rho, sigma, seed){
        set.seed(seed)
        theta <- c(rep(log(0.5*p), 5), rep(0, (p-5)))
        H <- abs(outer(1:p, 1:p, "-"))
        Sigma <- rho^H
        W <- rmvnorm(n, theta, Sigma)
        X <- W
        for(i in 1:n){ X[i,] <- exp(W[i,])/sum(exp(W[i,])) }
        Z <- log(X)
        err <- rnorm(n, mean=0, sd=sigma)
        y <- Z%*%tb + err
        return(list(X=X, Z=Z, y=y))
    }
    
    # Generate a small compositional dataset 
    dat <- gen.comp.data(p = 30, n = 50, tb = c(1, -1, rep(0, 28)), 
                         rho = 0.2, sigma = 0.5, seed = 1)
    
    # Run CSGL with a (very short) sequence of lambda values 
    res <- csgl(dat$Z, dat$y, groups = rep(1:6, each = 5), mu = 1, theta = 0.95, nlam = 3, min.frac = 0.1)
    res$beta 
    
    # Run CSGL with a sequence of lambda values, choosing lambda by 5-fold cross-validation 
    res.cv <- cv.csgl(dat$Z, dat$y, groups = rep(1:6, each = 5), mu = 1, theta = 0.95, nlam = 10, min.frac = 0.1, nfolds = 5)
    res.cv$cv.beta
    
    # Run CSGL with a sequence of lambda values, choosing lambda by GIC 
    res.gic <- gic.csgl(dat$Z, dat$y, groups = rep(1:6, each = 5), mu = 1, theta = 0.95, nlam = 10, min.frac = 0.1) 
    res.gic$gic.beta
    
