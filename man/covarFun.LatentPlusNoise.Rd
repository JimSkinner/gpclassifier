\name{covarFun.LatentPlusNoise}
\alias{covarFun.LatentPlusNoise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
covarFun.LatentPlusNoise(x, dk = NA, hp = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{dk}{
%%     ~~Describe \code{dk} here~~
}
  \item{hp}{
%%     ~~Describe \code{hp} here~~
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
function (x, dk = NA, hp = NA) 
{
    if (is(x, "CovarFun")) {
        k_f = x@k
        dk_f = x@dk
        hp_f = x@hp
    }
    else if (is(x, "function") & is(dk, "function") & is(hp, 
        "list")) {
        k_f = x
        dk_f = dk
        hp_f = hp
    }
    else {
        stop(paste("Wrong method signature; must be (CovarFun) or (function,", 
            "function, list)"))
    }
    hp = c(hp_f, lsf = 0, lsn = 0)
    k = function(.Object, x, y) {
        exp(.Object@hp$lsf) * k_f(.Object, x, y) + exp(.Object@hp$lsn) * 
            all(x == y)
    }
    dk = function(.Object, x, y) {
        c(dk_f(.Object, x, y), lsf = exp(.Object@hp$lsf) * k_f(.Object, 
            x, y), lsn = exp(.Object@hp$lsn) * all(x == y))
    }
    return(covarFun(k, dk, hp))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
