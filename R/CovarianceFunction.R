library(kernlab) # TODO: For package, where do these go?
library(functional)

#' Class CovarFun.
#'
#' Class \code{CovarFun} defines a covariance function to be used as part of
#' a \code{GPC} gaussian process classifier. A \code{CovarFun} object is made
#' up of a kernel k(x,y), a set of hyperparameters used in the kernel, and
#' a function returning the gradient of the kernel with respect to each of the
#' hyperparameters.
#' @export
CovarFun <- setClass(
  "CovarFun",
  slots    = c(k ="function", # Kernel
               dk="function", # Kernel gradient
               hp="list"),    # Hyperparameters
)

## Constructors
covarFun <- CovarFun

#' Construtor method of CovarFun class
#'
#' Creates a new CovarFun object intended to be used inside a \code{GPC}
#' Gaussian Process Classifier.
#'
#' A \code{CovarFun} object extends the kernel object which supplies a kernel
#' function along with a list of hyperparameters. \code{CovarFun} also supplies
#' a function returning the gradient of the kernel with respect to the
#' hyperparameters, such that the hyperparameters may be tuned by the \code{GPC}
#' class.
#'
#' @param k A kernel function (object, x, y) -> numeric() which, given data x
#' and y returns an inner product. Kernel hyperparameters may be accessed with
#'  object@@hp.
#' @param hp A list of kernel hyperparameters.
#' @param dk A function (object, x, y) -> list() which returns the gradient of
#'  k with respect to the hyprparameters in the form of a lsit of the same
#'  shape as the kyperparameters
#'
#' @examples
#' # Isotropic squared exponential covariance function with log length scale
#' # (ll) hyperparameter
#' k  = exp(-0.5 * sum(exp(-2*.Object@@hp$ll) * (x-y)^2))
#' hp = list(ll=0)
#' dk = list(ll=.Object@@k(.Object, x, y) * exp(-2*.Object@@hp$ll)*crossprod(x-y))
#' C  = CovarFun(k, hp, dk)
#' @export
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

#' Return kernel function k:x,y -> numeric(). Differs from the kernel function
#' specified when constructing the \code{CovarFun}, since the kernel function
#' returned only requires parameters x,y, not object.
setGeneric(name = "getKernel",
           def  = function(object) standardGeneric("getKernel"))

#' Return kernel gradient function k:x,y -> list(). As with
#' getKernel(CovarFun), this differs from the kernel gradient function
#' specified when constructing the \code{CovarFun}, since the function
#' returned only requires parameters x,y, not object.
setGeneric(name = "getKernelGrad",
           def  = function(object) standardGeneric("getKernelGrad"))

#' Return the list of hyperparameters.
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

#' Change the \code{CovarFun} hyperparameter list to the list supplied.
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

#' Augment \code{CovarFun} with latent function and noise hyperparameters
#'
#' Creates a covariance function k which augments some covariance function k_f
#' with two extra hyperparameters lsf and lsn such that:
#' k(x, y) = exp(lsf)*k_f(x, y) + exp(lsn)*I[x=y]
#' Automatically provides derivatives of lsf, lsn. Used for neater covariance
#' funtion specifications. lsf, lsn stand for log(sigma^2_f) and log(sigma^2_n),
#' where 'f' labels the magnitude parameter and kernel for the latent
#' function, whilst 'n' labels the magnitude parameter for the noise.
#'
#' @param covarFun \code{CovarFun} object to augment.
#' @return covarFun augmented with lsf, lsn hyperparameterss for signal and
#'  noise magnitude.
#' @examples
#' C = covarFun.LatendPlusNoise(covarFun.SE(0))
#' @export
covarFun.LatentPlusNoise <- function(covarFun) {
  stopifnot(is(covarFun, "CovarFun"))
  k_f  = covarFun@k
  dk_f = covarFun@dk
  hp_f = covarFun@hp

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

#' Squared Exponential covariance function.
#' @param ll log length scale hyperparameter.  \code{ll} must be of length 1 or
#' the same length as the input vectors. If length 1 then this is an
#' isotropic SE covar fun, else there is a length scale for each dimension
#' of the input.
#' @return \code{CovarFun} for a squared exponential covariance function
#' @examples
#' C <- covarFun.se(0)
#' @export
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
