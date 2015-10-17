EP.fn <- function(object) {
  ## Algorithm 3.5 in Rasmussen
  # Variable naming conventions:
  # Cavity params are suffixed _cav
  # Site/local params are suffixed _loc
  # Global params are suffixed _g

  n       = nrow(object@X)            # Number of samples
  nu_loc  = matrix(0, nrow=n, ncol=1) # Local nu;     nx1; \~\nu
  tau_loc = matrix(0, nrow=n, ncol=1) # Local tau;    nx1; \~\tau
  mu_g    = matrix(0, nrow=n, ncol=1) # Global mu;    nx1; \mu
  Sigma_g = object@K                  # Global Sigma; nxn; \Sigma

  # Cavity parameters used to calculate Z
  tau_cav = matrix(0, nrow=n, ncol=1) # \tau_{-i}
  mu_cav  = matrix(0, nrow=n, ncol=1) # \mu_{-i}

  Y = ifelse(object@Y, 1, -1)
  converged = FALSE
  while(!converged) {
    tau_loc_old = tau_loc
    for (i in sample(n)) { # Updating in a random order can improve convergence
      var_g   = Sigma_g[i,i]
      # Calculate cavity parameters; \tau_{-i}, \nu_{-i} (Eq. 3.56)
      tau_cav[i] = 1/var_g - tau_loc[i]
      nu_cav  = mu_g[i]/var_g - nu_loc[i]

      var_cav = 1/tau_cav[i]
      mu_cav[i]  = nu_cav/tau_cav[i]
      # Params for gaussian marginal (\hat{q}) approximating poduct of
      # cavity distribution and exact likelihood (3.57)
      z       = (Y[i]*mu_cav[i]) / sqrt(1 + var_cav)
      Z_loc   = pnorm(z)
      mu_hat  = mu_cav[i] +
                (Y[i]*var_cav*dnorm(z)) /
                (Z_loc*sqrt(1 + var_cav))
      var_hat = var_cav -
                (var_cav^2 * dnorm(z))/((1+var_cav)*Z_loc) *
                (z + dnorm(z)/Z_loc)

      # Update local params
      tau_loc_delta = 1/var_hat - tau_cav[i] - tau_loc[i]
      tau_loc[i]    = tau_loc[i] + tau_loc_delta
      nu_loc[i]     = mu_hat/var_hat - nu_cav
      Sigma_g       = Sigma_g - (
                        tau_loc_delta / (1 + tau_loc_delta*var_g)
                      ) * crossprod(Sigma_g[i, , drop = FALSE])
      mu_g          = Sigma_g %*% nu_loc

      if (Z_loc==0) {
        stop("Cannot determine site parameters here; local marginal likelihood is 0")
      }
    }

    rootS_loc = diag(as.vector(sqrt(tau_loc)))
    R         = chol(diag(n) + rootS_loc %*% object@K %*% rootS_loc)
    V         = backsolve(R, rootS_loc %*% object@K, transpose=TRUE)
    Sigma_g   = object@K - crossprod(V)
    mu_g      = Sigma_g %*% nu_loc

    converged = sqrt(sum((tau_loc - tau_loc_old)^2)) < object@tol
  }

  ## TODO: Can I move logZ_EP calculation to the getLml() method?

  ## Calculate approximate log likelihood; log Z_{EP}. Uses Section
  ## 3.6.3, Rasmussen.
  # Eq 3.73
  logZ_EPa = 0.5*sum(log(1+ tau_loc/tau_cav)) - sum(log(diag(R)))

  # Term 3 in Eq 3.65; does not need modifying for speed/stability.
  logZ_EPb = sum(log(pnorm(Y*mu_cav / sqrt(1+1/tau_cav))))

  # Eq 3.74
  logZ_EPc = 0.5 * t(nu_loc) %*% (
    object@K -
    crossprod(forwardsolve(t(R), rootS_loc%*%object@K)) -
    diag((1/c(tau_cav + tau_loc)))
  ) %*% nu_loc

  # Eq 3.75
  logZ_EPd = 0.5 * crossprod(
    (mu_cav*tau_cav)/(tau_cav + tau_loc),
    mu_cav*tau_loc - 2*nu_loc)

  # Calculate gradients of log marginal likelihood with respect to
  # hyperparameters.
  #
  # Implementation of Algorithm 5.2 from Rasmussen (Secion 5.5)

  # TODO: Check these are the same, and which is faster.
  #R         = chol(diag(nrow(X)) + rootS_loc %*% K %*% rootS_loc)
  #R2        = chol(diag(nrow(X)) +
  #                 tcrossprod(sqrt(object@siteParams$tau_loc))*K)
  b = nu_loc -
      rootS_loc %*% backsolve(R,
        backsolve(R, rootS_loc %*% object@K %*% nu_loc,
                  transpose=TRUE))
  F_ = tcrossprod(b) -
       crossprod(backsolve(R, rootS_loc, transpose=TRUE))

  C   <- getCovarFun(object)
  dk  <- getKernelGrad(C)
  nHP <- length(unlist(getHP(C)))
  dK  <- array(0, dim=c(nHP, nrow(X), nrow(X)))
  for (i in 1:nrow(X)) {
    for (j in 1:i) {
      gradList <- dk(X[i,], X[j,])
      gradVec  <- unlist(gradList)
      dK[,i,j] <- dK[,j,i] <- gradVec
    }
  }

  dlml = vapply(1:nHP, function(hpInd) {
    sum(diag(F_ %*% dK[hpInd,,])) # TODO: Only compute diagonal
  }, numeric(1)) * 0.5

  return(list(tau_loc = as.numeric(tau_loc),
              nu_loc  = as.numeric(nu_loc),
              lml     = as.numeric(logZ_EPa + logZ_EPb + logZ_EPc + logZ_EPd),
              dlml    = dlml))
}
