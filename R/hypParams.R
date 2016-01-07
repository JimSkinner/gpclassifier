#' Return maximum likelihood covariance function hyperparameters
#' @import optimx
#' @import memoise
#' @import mvtnorm
setGeneric(name="tune",
           def=function(object) standardGeneric("tune"))

tune_theta <- function(object) {
  N <- nrow(object@X)
  D <- ncol(object@X)

  converged = FALSE
  iter = 0
  while (!converged) {
    # Coordiante descend theta_2
    theta2 <- tune_theta_gp(object)

    # Coordiante descend theta_1
    FL        <- tune_theta_fl(object)
    theta1    <- FL$theta1
    transform <- FL$transform
    V         <- FL$V

    # converged?
    iter = iter + 1
    if (iter > 10) break()
  }
  return(list(theta1=theta1, theta2=theta2, transform=transform, V=V))
}

tune_theta_fl <- function(object) {
  pars <- gpc_ppca(object)
  return(list(theta1   =pars[c('W','sigSq','mu')],
              transform=pars$transform,
              V        =pars$V))
}

gpc_ppca <- function(object) {
  N = nrow(object@X)
  D = ncol(object@X)
  K = object@phi1$K

  nParams = (
    D*K + # W
    1 +   # log(sigma^2)
    D     # mu
  )

  # Log likelihood
  L  = function(params) {
    W     = matrix(params[1:(D*K)], nrow=D, ncol=K)
    sigSq = exp(params[D*K + 1])
    mu    = params[(1:D)+D*K+1]

    # It turns out that the numeric stability gained by inverting C with a
    # cholesky decomposition is very important; enough to make the difference
    # between the optimisation succeeding/failing.
    M  = crossprod(W) + sigSq*diag(K)
    tryCatch({
      R2 = chol(M)
    }, error = function(err) {
      browser()
    })
    LinvWt = forwardsolve(t(R2), t(W))

    # A stable way of calculating the log determinant of C is also important.
    logDetC = sum(log(c(svd(W)$d^2, rep(0, D-K)) + rep(sigSq, D)))

    .a = -N*D*0.5*log(2*pi)
    .b = -N*0.5*logDetC
    .c = -0.5*sum(apply(X, 1, function(x) {
      crossprod(x-mu) - crossprod(LinvWt%*%(x-mu))
    }))/sigSq
    logLik_ppca = .a + .b + .c

    ## Calculate GPC likelihood
    M = crossprod(W) + sigSq*diag(K)
    V = tcrossprod(sweep(X, 2, mu), solve(M, t(W)))
    object@V   <- V
    ep         <- EP(object)
    logLik_gpc <- ep$lml #+ sum(apply(V, 1, function(v) log(dmvnorm(v))))

    return(logLik_ppca + logLik_gpc)
  }

  if (length(object@theta1)) { # theta1 has elements; use these as initial pars
    theta0 = unlist(object@theta1[c('W','sigSq','mu')])
  } else { # Come up with initial parameter set
    #theta0 = c(eigen(crossprod(X))$vectors[,1:K], # W_0
    theta0 = c(rnorm(D*K),
               1, # sigSq_0
               rep(0,D)) # mu_0
  }
  attempts = 10
  for (attempt in 1:attempts) {
    print(attempt)
    optObj = optimx(par= theta0,
                    fn = L,
                    method="BFGS",
                    control=list(kkt=FALSE,
                                 starttests=FALSE,
                                 trace=0,
                                 maximize=TRUE,
                                 follow.on=FALSE,
                                 save.failures=FALSE,
                                 usenumDeriv=FALSE))
    par = coef(optObj)
    if (is.na(par[1])) {
      theta0 = rnorm(D*K + D + 1)
    } else {
      break
    }
  }

  W     = matrix(par[1:(D*K)], nrow=D, ncol=K)
  sigSq = exp(par[D*K + 1])
  mu    = par[(1:D)+D*K+1]

  # Linear transformation mapping raw data to feature space.
  M         = crossprod(W) + sigSq*diag(K)
  MinvW     = solve(M, t(W))
  transform = function(X_) tcrossprod(sweep(X_, 2, mu), MinvW)

  # Feature representation of X
  V         = transform(X)

  return(list(W        =W,
              sigSq    =sigSq,
              transform=transform,
              mu       =mu,
              V        =V))

}

tune_theta_gp <- function(object) {
  ## Information shared between objective and gradient functions
  C = getCovarFun(object)
  relist.template <- getHP(C)

  currTheta = NA # TODO: I don't think this is needed?

  # TODO: Onlt remember the last 4 or so?
  shared_function <- memoise(function(theta) { # TODO: Move all this out into a new function that calculates dlml
    setHP(C)        <- relist(theta, relist.template)
    object@covarFun <- C
    object@K        <- kernelMatrix(getKernel(C), object@V)

    ep             <- EP(object)
    object@tau_loc <- ep$tau_loc
    object@nu_loc  <- ep$nu_loc
    object@lml     <- ep$lml
    object@dlml    <- ep$dlml

    currTheta      <- theta
    return(object)
  })

  objective_function = function(theta) {  # Function to optimise
    object <- shared_function(theta)
    object@lml
  }

  gradient_function = function(theta) {  # Gradient of objective function
    object <- shared_function(theta)
    object@dlml
  }

  theta0 = unlist(getHP(C))

  optObj = optimx(par     = theta0,
                  fn      = objective_function,
                  gr      = gradient_function,
                  method  = 'CG',
                  hessian = FALSE,
                  control = list(maximize=TRUE, starttests=FALSE,
                                 kkt=FALSE))
  return(setHP(C) <- relist(coef(optObj), relist.template))
}

setMethod(f         = "tune",
          signature = "GPC",
          def       = tune_theta)
