\name{tuneHPs.fn}
\alias{tuneHPs.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tuneHPs.fn(object)
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
    C = getCovarFun(object)
    relist.template <- getHP(C)
    currTheta = NA
    library("memoise")
    shared_function <- memoize(function(theta) {
        setHP(C) <- relist(theta, relist.template)
        object@covarFun <- C
        object@K <- kernelMatrix(getKernel(C), object@X)
        ep <- EP(object)
        object@tau_loc <- ep$tau_loc
        object@nu_loc <- ep$nu_loc
        object@lml <- ep$lml
        object@dlml <- ep$dlml
        currTheta <- theta
        return(object)
    })
    objective_function = function(theta) {
        object <- shared_function(theta)
        object@lml
    }
    gradient_function = function(theta) {
        object <- shared_function(theta)
        object@dlml
    }
    theta0 = unlist(getHP(C))
    optObj = optimx(par = theta0, fn = objective_function, gr = gradient_function, 
        method = "CG", hessian = FALSE, control = list(maximize = TRUE, 
            starttests = FALSE, kkt = FALSE))
    return(setHP(C) <- relist(coef(optObj), relist.template))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
