\name{CovarFun-class}
\Rdversion{1.1}
\docType{class}
\alias{CovarFun-class}
\alias{CovarFun}
\alias{getHP,CovarFun-method}
\alias{getKernel,CovarFun-method}
\alias{getKernelGrad,CovarFun-method}
\alias{initialize,CovarFun-method}
\alias{setHP<-,CovarFun-method}

\title{Class \code{"CovarFun"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{CovarFun(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{k}:}{Object of class \code{"function"} ~~ }
    \item{\code{dk}:}{Object of class \code{"function"} ~~ }
    \item{\code{hp}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{getHP}{\code{signature(object = "CovarFun")}: ... }
    \item{getKernel}{\code{signature(object = "CovarFun")}: ... }
    \item{getKernelGrad}{\code{signature(object = "CovarFun")}: ... }
    \item{initialize}{\code{signature(.Object = "CovarFun")}: ... }
    \item{setHP<-}{\code{signature(object = "CovarFun")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("CovarFun")
}
\keyword{classes}
\usage{
CovarFun(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
