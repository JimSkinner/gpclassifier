source("CovarianceFunction.R")
source("siteParams.R")
source("hypParams.R")

## TODO: Keep documented generics here, move (undocumented) methods into a new
## file.

#' Class GPC.
#'
#' Class \code{GPC} defines a Gaussian Process Classifier.
#' @export
GPC <- setClass(
  "GPC",
  slots = c(
    covarFun = "CovarFun",     # Covariance function
    X        = "matrix",       # Input data
    Y        = "logical",      # Binary labels
    K        = "kernelMatrix", # Gram matrix for data and covariance function
    lml      = "numeric",      # (approximate) log marginal likelihood Z_EP
    dlml     = "numeric",      # Derivatives of lml wrt parameters
    nu_loc   = "numeric",      # Site parameters
    tau_loc  = "numeric",      # ""
    tol      = "numeric"       # Tolerance for site param convergence criteria
  ),
  prototype  = prototype(nu_loc   = numeric(0),
                         tau_loc  = numeric(0),
                         tol      = 2^-10,
                         covarFun = covarFun.SE()),
  validity   = function(object) {
    if (nrow(object@X) != length(object@Y)) {
      return("Y must have the same number of elements as X (samples in rows)")
    } else if (length(object@nu_loc)  != nrow(object@X) |
               length(object@tau_loc) != nrow(object@X)) {
      return(paste("Site parameters nu_loc & tau_loc must be of the same",
                   "length as the numer of samples"))
    }
    return(TRUE)
  },
)

#' Construtor method of GPC class
#'
#' Creates a new GPC object used to predict labels for new input data.
#'
#' Covariance function hyperparameters are selected automatically through
#' maximum likelihood. This class is implemented using the Expectation
#' Propagation approximation detailed in (Gaussian Processes for Machine
#' Learning, Rasmussen and Williams, 2006).
#'
#' @param X Matrix of input data; sample in rows.
#' @param Y Logical vector of binary labels.
#' @param covarFun Covariance function to use. Must be of class CovarFun. If
#'  omitted, the squared exponential covariance function is used by default.
#'
#' @return S4 object of class GPC, where covarance function hyperparameters
#'  have been set to their maximum likelihood estimates.
#'
#' @examples
#' # Create synthetic dataset
#' X <- matrix(rnorm(60), ncol=2)
#' Y <- rowSums(X^2) < 1
#'
#' # New GPX Object with default squared exponential covariance function.
#' gpc <- GPC(X, Y)
#'
#' # Predict labels for new data
#' Xst <- matrix(rnorm(60), ncol=2)
#' Yst <- predict(gpc, Xst)
#' @export
setMethod(f          = "initialize",
          signature  = "GPC",
          definition = function(.Object, X, Y, covarFun) {
            updateList = list()
            if (!missing(X))        {updateList$X = X}
            if (!missing(Y))        {updateList$Y = Y}
            if (!missing(covarFun)) {updateList$covarFun = covarFun}

            update(.Object) <- updateList
            validObject(.Object) # Maybe unneccesary
            return(.Object)
          }
)

#######################
## Methods
#######################

#' Calculate site parameters, likelihood and likelihood gradient using
#' Expectation Propagation
#' @export
setGeneric(name="EP",
           def=function(object) standardGeneric("EP"))
setMethod(f         = "EP",
          signature = "GPC",
          def       = GPC.EP)

#' Return maximum likelihood covariance function hyperparameters
#' @export
setGeneric(name="hpTune",
           def=function(object) standardGeneric("hpTune"))
setMethod(f         = "hpTune",
          signature = "GPC",
          def       = GPC.hpTune)

#' Update \code{X}, \code{Y} or \code{covarFun}
#'
#' Update the data \code{X}, labels \code{Y} or covariance function
#' \code{covarFun}, causing a recalculation of hyperparameters,
#' site parameters and covariance matrix.
#' @export
setGeneric(name="update<-",
           def=function(object, value) standardGeneric("update<-"))

setReplaceMethod(
  f         = "update",
  signature = c("GPC", "list"),
  def       = function(object, value) {
    if ('X' %in% names(value)) {object@X <- value$X}
    if ('Y' %in% names(value)) {object@Y <- value$Y}
    if ('covarFun' %in% names(value)) {object@covarFun <- value$covarFun}

    if (nrow(object@X) == 0) { # No data
      object@tau_loc <- numeric(0)
      object@nu_loc  <- numeric(0)
      object@lml     <- numeric(0)
      object@dlml    <- numeric(0)
      object@K       <- as.kernelMatrix(matrix(numeric(0)))
    } else {
      setHP(object@covarFun) <- tuneHPs(object)
      object@K        <- kernelMatrix(getKernel(object@covarFun), object@X)

      ep              <- EP(object)
      object@tau_loc  <- ep$tau_loc
      object@nu_loc   <- ep$nu_loc
      object@lml      <- ep$lml
      object@dlml     <- ep$dlml
    }
    return(object)
  }
)

setMethod(f         = "predict",
          signature = "GPC",
          def       = function(object, Xst) {
            # If new points unspecified, return predictions on training data
            if (missing(Xst)) {
              Xst <- object@X
            } else if (!is.matrix(Xst) || (ncol(Xst)!=ncol(object@X))) {
              stop(paste("Xst must be a matrix with the same number of columns",
                         "as the original data"))
            }

            # Algorithm 3.6, Rasmussen GP book.
            # Syntax used here: Xst, yst correspond to X*, y*; the data and
            # predictive probability of a new sample.
            rootS_loc = diag(as.vector(sqrt(object@tau_loc)))
            R = chol(diag(nrow(object@X)) + rootS_loc %*% object@K %*% rootS_loc)
            z = rootS_loc %*% backsolve(R,
              forwardsolve(t(R), rootS_loc %*% object@K %*% object@nu_loc))

            Yst = numeric(nrow(Xst))
            for (xind in 1:nrow(Xst)) { # TODO: Vectorize
              xst = Xst[xind,,drop=FALSE]

              ## Kernel matrix between training samples and test samples
              Kst = kernelMatrix(getKernel(object@covarFun), xst, object@X)

              # Approximate latent variable for test input vector x*
              fst = Kst %*% (object@nu_loc-z)

              # Variance of f*. v is just a working variable.
              # TODO: Can make crossprod(rootS_loc, Kst) more efficient bc diag
              v = forwardsolve(t(R), tcrossprod(rootS_loc, Kst))
              fst_var = kernelMatrix(getKernel(object@covarFun), xst) - crossprod(v)

              # Predictive class probability
              Yst[xind] = pnorm(fst/sqrt(1+fst_var))
            }
            return(Yst)
          }
)

setMethod(f="fitted", signature="GPC", def=function(object) predict(.Object))

#######################
## Accessors
#######################

## Getters

#' Return log marginal likelihood.
#' @export
setGeneric(name = "getLml",
           def  = function(object) standardGeneric("getLml"))

#' Return partial derivatives of the log marginal likelihood
#' with respect to each of the hyperparameters.
#' @export
setGeneric(name = "getDLml",
           def  = function(object) standardGeneric("getDLml"))

#' GPC Return the covariance matrix; the matrix of inner products
#' between each pair of data points in the space induced by the covariance
#' function.
#' @export
setGeneric(name = "getK",
           def  = function(object) standardGeneric("getK"))

#' GPC Return the covariance function (class \code{CovarFun}) used.
#' @export
setGeneric(name = "getCovarFun",
           def  = function(object) standardGeneric("getCovarFun"))

#setGeneric(name = "getHP", # Do not need: generic is set in CovarianceFunction
#           def  = function(object) standardGeneric("getHP"))

## TODO: Move these out to a new file GPC-methods.R
setMethod(f         = "getLml",
          signature = "GPC",
          def       = function(object) {object@lml}
)

setMethod(f         = "getDLml",
          signature = "GPC",
          def       = function(object) {object@dlml}
)

setMethod(f         = "getK",
          signature = "GPC",
          def       = function(object) {object@K}
)

setMethod(f         = "getCovarFun",
          signature = "GPC",
          def       = function(object) {object@covarFun}
)

setMethod(f         = "getHP",
          signature = "GPC",
          def       = function(object) {getHP(getCovarFun(object))}
)

## Setters

## TODO: Implement? Maybe not; just have the user construct a new GPC obj
## since that is what will pretty much happen anyway.

#setGeneric(name = "setData<-",
#           def  = function(object, value) standardGeneric("setData<-"))
##setGeneric(name = "setHP<-", ## Again; do not need bc of covarfun
##           def  = function(object, value) standardGeneric("setHP<-"))
#setGeneric(name = "setCovarFun<-",
#           def  = function(object, value) standardGeneric("setCovarFun<-"))
#setGeneric(name = "setSiteParams<-",
#           def  = function(object, value) standardGeneric("setSiteParams<-"))
#setGeneric(name = "setK<-",
#           def  = function(object, value) standardGeneric("setK<-"))
#
#setReplaceMethod(f         = "setData",
#                 signature = c("GPC", "list"),
#                 def       = function(object, value) {
#                   stopifnot(setequal(names(value), c("X", "Y")))
#                   update(object) <- value
#                   validObject(object)
#                   return(object)
#                 }
#)
#
#setReplaceMethod(f         = "setHP",
#                 signature = c("GPC", "list"),
#                 def       = function(object, value) {
#                   if (!setequal(value, getHP(object))) {
#                     newCovarFun = setHP(object@covarFun) <- value
#                     update(object) <- list(covarFun=newCovarFun)
#                     validObject(object)
#                   }
#                   return(object)
#                 }
#)
#
#setReplaceMethod(f         = "setCovarFun",
#                 signature = c("GPC", "CovarFun"),
#                 def       = function(object, value) {
#                   if (!identical(value, object@covarFun)) {
#                     object@covarFun <- value
#                     validObject(object@covarFun)
#                     update(object)
#                   }
#                   return(object)
#                 }
#)
#
#setReplaceMethod(f         = "setSiteParams",
#                 signature = c("GPC", "list"),
#                 def       = function(object, value) {
#                   object@nu_loc  <- value$nu_loc
#                   object@tau_loc <- value$tau_loc
#                   validObject(object)
#                   return(object)
#                 }
#)
#
#setReplaceMethod(f         = "setK",
#                 signature = c("GPC", "matrix"),
#                 def       = function(object, value) {
#                   object@K <- as.kernelMatrix(value)
#                   validObject(object)
#                   return(object)
#                 }
#)
