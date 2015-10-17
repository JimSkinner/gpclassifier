\name{EP.fn}
\alias{EP.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EP.fn(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object) 
{
    n = nrow(object@X)
    nu_loc = matrix(0, nrow = n, ncol = 1)
    tau_loc = matrix(0, nrow = n, ncol = 1)
    mu_g = matrix(0, nrow = n, ncol = 1)
    Sigma_g = object@K
    tau_cav = matrix(0, nrow = n, ncol = 1)
    mu_cav = matrix(0, nrow = n, ncol = 1)
    Y = ifelse(object@Y, 1, -1)
    converged = FALSE
    while (!converged) {
        tau_loc_old = tau_loc
        for (i in sample(n)) {
            var_g = Sigma_g[i, i]
            tau_cav[i] = 1/var_g - tau_loc[i]
            nu_cav = mu_g[i]/var_g - nu_loc[i]
            var_cav = 1/tau_cav[i]
            mu_cav[i] = nu_cav/tau_cav[i]
            z = (Y[i] * mu_cav[i])/sqrt(1 + var_cav)
            Z_loc = pnorm(z)
            mu_hat = mu_cav[i] + (Y[i] * var_cav * dnorm(z))/(Z_loc * 
                sqrt(1 + var_cav))
            var_hat = var_cav - (var_cav^2 * dnorm(z))/((1 + 
                var_cav) * Z_loc) * (z + dnorm(z)/Z_loc)
            tau_loc_delta = 1/var_hat - tau_cav[i] - tau_loc[i]
            tau_loc[i] = tau_loc[i] + tau_loc_delta
            nu_loc[i] = mu_hat/var_hat - nu_cav
            Sigma_g = Sigma_g - (tau_loc_delta/(1 + tau_loc_delta * 
                var_g)) * crossprod(Sigma_g[i, , drop = FALSE])
            mu_g = Sigma_g \%*\% nu_loc
            if (Z_loc == 0) {
                stop("Cannot determine site parameters here; local marginal likelihood is 0")
            }
        }
        rootS_loc = diag(as.vector(sqrt(tau_loc)))
        R = chol(diag(n) + rootS_loc \%*\% object@K \%*\% rootS_loc)
        V = backsolve(R, rootS_loc \%*\% object@K, transpose = TRUE)
        Sigma_g = object@K - crossprod(V)
        mu_g = Sigma_g \%*\% nu_loc
        converged = sqrt(sum((tau_loc - tau_loc_old)^2)) < object@tol
    }
    logZ_EPa = 0.5 * sum(log(1 + tau_loc/tau_cav)) - sum(log(diag(R)))
    logZ_EPb = sum(log(pnorm(Y * mu_cav/sqrt(1 + 1/tau_cav))))
    logZ_EPc = 0.5 * t(nu_loc) \%*\% (object@K - crossprod(forwardsolve(t(R), 
        rootS_loc \%*\% object@K)) - diag((1/c(tau_cav + tau_loc)))) \%*\% 
        nu_loc
    logZ_EPd = 0.5 * crossprod((mu_cav * tau_cav)/(tau_cav + 
        tau_loc), mu_cav * tau_loc - 2 * nu_loc)
    b = nu_loc - rootS_loc \%*\% backsolve(R, backsolve(R, rootS_loc \%*\% 
        object@K \%*\% nu_loc, transpose = TRUE))
    F_ = tcrossprod(b) - crossprod(backsolve(R, rootS_loc, transpose = TRUE))
    C <- getCovarFun(object)
    dk <- getKernelGrad(C)
    nHP <- length(unlist(getHP(C)))
    dK <- array(0, dim = c(nHP, nrow(X), nrow(X)))
    for (i in 1:nrow(X)) {
        for (j in 1:i) {
            gradList <- dk(X[i, ], X[j, ])
            gradVec <- unlist(gradList)
            dK[, i, j] <- dK[, j, i] <- gradVec
        }
    }
    dlml = vapply(1:nHP, function(hpInd) {
        sum(diag(F_ \%*\% dK[hpInd, , ]))
    }, numeric(1)) * 0.5
    return(list(tau_loc = as.numeric(tau_loc), nu_loc = as.numeric(nu_loc), 
        lml = as.numeric(logZ_EPa + logZ_EPb + logZ_EPc + logZ_EPd), 
        dlml = dlml))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
