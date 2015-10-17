\name{covarFun.SE}
\alias{covarFun.SE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
covarFun.SE(ll = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ll}{
%%     ~~Describe \code{ll} here~~
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
function (ll = 0) 
{
    k = function(.Object, x, y) {
        exp(-0.5 * sum(exp(-2 * .Object@hp$ll) * (x - y)^2))
    }
    dk = function(.Object, x, y) {
        if (length(.Object@hp$ll) == 1) {
            list(ll = .Object@k(.Object, x, y) * exp(-2 * .Object@hp$ll) * 
                crossprod(x - y))
        }
        else {
            list(ll = .Object@k(.Object, x, y) * exp(-2 * .Object@hp$ll) * 
                (x - y)^2)
        }
    }
    hp = list(ll = ll)
    return(covarFun(k, dk, hp))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
