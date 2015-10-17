library("optimx")

tuneHPs.fn <- function(object) {

  ## Information shared between objective and gradient functions
  C = getCovarFun(object)
  relist.template <- getHP(C)

  currTheta = NA

  library("memoise") # TODO: Onlt remember the last 4 or so?
  shared_function <- memoize(function(theta) { # TODO: Move all this out into a new function that calculates dlml
    setHP(C)        <- relist(theta, relist.template)
    object@covarFun <- C
    object@K        <- kernelMatrix(getKernel(C), object@X)

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

  #browser()

#  x = seq(-2, 2, length=10)
#  y1  = apply(rbind(x,0,0)+theta0, 2, objective_function)
#  y2  = apply(rbind(0,x,0)+theta0, 2, objective_function)
#  y3  = apply(rbind(0,0,x)+theta0, 2, objective_function)
#
#  dy1 = apply(rbind(x,0,0)+theta0, 2, gradient_function)[1,]
#  dy2 = apply(rbind(0,x,0)+theta0, 2, gradient_function)[2,]
#  dy3 = apply(rbind(0,0,x)+theta0, 2, gradient_function)[3,]
#
#  par(mfcol=c(3,1))
#  plot(x, y1-mean(y1), type='l', ylim=range(y1-mean(y1), dy1))
#  lines(x, dy1, col='red')
#  lines(range(x), c(0,0), lty=2, col='grey')
#  plot(x, y2-mean(y2), type='l', ylim=range(y2-mean(y2), dy2))
#  lines(x, dy2, col='red')
#  lines(range(x), c(0,0), lty=2, col='grey')
#  plot(x, y3-mean(y3), type='l', ylim=range(y3-mean(y3), dy3))
#  lines(x, dy3, col='red')
#  lines(range(x), c(0,0), lty=2, col='grey')
#  browser()
  optObj = optimx(par     = theta0,
                  fn      = objective_function,
                  gr      = gradient_function,
                  method  = 'CG',
                  hessian = FALSE,
                  control = list(maximize=TRUE, starttests=FALSE,
                                 kkt=FALSE))
  return(setHP(C) <- relist(coef(optObj), relist.template))
}

