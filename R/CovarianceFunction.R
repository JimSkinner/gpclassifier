library(kernlab) # TODO: For package, where do these go?
library(functional)

CovarFun <- setClass(
  "CovarFun",
  slots    = c(k ="function", # Kernel
               dk="function", # Kernel gradient
               hp="list"),    # Hyperparameters
)

## Constructors
covarFun <- CovarFun

setMethod(f          = "initialize",
          signature  = "CovarFun",
          definition = function(.Object, k, dk, hp) {
            if (!missing(k))  .Object@k  <- k
            if (!missing(dk)) .Object@dk <- dk
            if (!missing(hp)) .Object@hp <- hp
            validObject(.Object)
            return(.Object)
          }
)

## Getters
setGeneric(name = "getKernel",
           def  = function(object) standardGeneric("getKernel"))

setGeneric(name = "getKernelGrad",
           def  = function(object) standardGeneric("getKernelGrad"))

setGeneric(name = "getHP",
           def  = function(object) standardGeneric("getHP"))

setMethod(f         = "getKernel",
          signature = "CovarFun",
          def       = function(object) {
            return(new("kernel", .Data=Curry(object@k, object), kpar=object@hp))
          }
)

setMethod(f         = "getKernelGrad",
          signature = "CovarFun",
          def       = function(object) return(Curry(object@dk, .Object=object))
)

setMethod(f         = "getHP",
          signature = "CovarFun",
          def       = function(object) {object@hp}
)

## Setters
setGeneric(name = "setHP<-",
           def  = function(object, value) {
             standardGeneric("setHP<-")
           }
)

setReplaceMethod(f         = "setHP",
                 signature = "CovarFun",
                 def       = function(object, value) {
                   object@hp <- value
                   validObject(object)
                   return(object)
                 }
)

###################################
## More specific class constructors
###################################

## General covariance function k_f which augments some covariance function k
## with two extra hyperparameters lsf and lsn such that:
## k_f(x, y) = exp(lsf)*k(x, y) + exp(lsn)*I[x=y]
##
## Automatically provides derivatives of lsf, lsn. Used for neater covariance
## funtion specifications. lsf, lsn stand for log(sigma^2_f) and log(sigma^2_n),
## where 'f' labels the magnitude parameter for the latent function, whilst 'n'
## labels the magnitude parameter for the noise.
##
## lsf - log sigma^2_f
## lsn - log sigma^2_n
covarFun.LatentPlusNoise <- function(x, dk=NA, hp=NA) {
  if (is(x, "CovarFun")) {
    k_f  = x@k
    dk_f = x@dk
    hp_f = x@hp
  } else if (is(x, "function") & is(dk, "function") & is(hp, "list")) {
    k_f  = x
    dk_f = dk
    hp_f = hp
  } else {
    stop(paste("Wrong method signature; must be (CovarFun) or (function,",
               "function, list)"))
  }

  # Augment hyperparameters with log signal variance, log noise variance
  hp = c(hp_f, lsf=0, lsn=0)

  # Augment covariance function with scaling of latent function + homoskedastic noise
  k = function(.Object, x, y) {
    exp(.Object@hp$lsf)*k_f(.Object, x, y) + exp(.Object@hp$lsn)*all(x==y)
  }

  # Augment covariance function derivatives with derivatives for lsf, lsn
  dk = function(.Object, x, y) {
    c(dk_f(.Object, x, y),
      "lsf" = exp(.Object@hp$lsf) * k_f(.Object, x, y),
      "lsn" = exp(.Object@hp$lsn) * all(x==y)
      #"lsn" = exp(.Object@hp$lsn) * isTRUE(all.equal(x, y))
    )
  }

  return(covarFun(k, dk, hp))
}

## TODO: Move CovarFun constructors to a new file
## TODO: I can save time by caching the exponentiated hyperparameters.

## TODO: Generalise this to any PosSemDef L

## Squared Exponential covariance function.
## ll  - log length scale
##
## ll must be of length 1 or the same length as the input vectors. If length
## 1 then this is an isotropic SE covar fun, else there is a length scale
## for each dimension.
covarFun.SE <- function(ll=0) {
  k = function(.Object, x, y) { # TODO: Robust against ll=-Inf?
    # This is safter but slower! :(
    #if (!(length(.Object@hp$ll) %in% c(1, length(x)))) {
    #  stop(paste("The length of the 'll' hyperparameter in the SE covariance",
    #             "function must be equal to 1 or the length of the input.",
    #             "length(ll) =", length(.Object@hp$ll), ", length(x) =",
    #             length(x)))
    #}
    exp(-0.5 * sum(exp(-2*.Object@hp$ll) * (x-y)^2))
  }
  dk = function(.Object, x, y) {
    if (length(.Object@hp$ll) == 1) { # Isotropic SE
      list(ll=.Object@k(.Object, x, y) * exp(-2*.Object@hp$ll)*crossprod(x-y))
    } else { # 1 param per dimension
      list(ll=.Object@k(.Object, x, y) * exp(-2*.Object@hp$ll)*(x-y)^2)
    }
  }
  hp = list(ll=ll)
  return(covarFun(k, dk, hp))
}

## TODO: Having classes for covariance function k and dk would be a LOT safer.
