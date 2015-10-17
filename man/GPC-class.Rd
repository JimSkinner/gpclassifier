\name{GPC-class}
\Rdversion{1.1}
\docType{class}
\alias{GPC-class}
\alias{GPC}
\alias{EP,GPC-method}
\alias{fitted,GPC-method}
\alias{getCovarFun,GPC-method}
\alias{getDLml,GPC-method}
\alias{getHP,GPC-method}
\alias{getK,GPC-method}
\alias{getLml,GPC-method}
\alias{initialize,GPC-method}
\alias{predict,GPC-method}
\alias{tuneHPs,GPC-method}
\alias{update<-,GPC,list-method}

\title{Class \code{"GPC"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{GPC(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{covarFun}:}{Object of class \code{"CovarFun"} ~~ }
    \item{\code{X}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{Y}:}{Object of class \code{"logical"} ~~ }
    \item{\code{K}:}{Object of class \code{"kernelMatrix"} ~~ }
    \item{\code{lml}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dlml}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{nu_loc}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{tau_loc}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{tol}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{EP}{\code{signature(object = "GPC")}: ... }
    \item{fitted}{\code{signature(object = "GPC")}: ... }
    \item{getCovarFun}{\code{signature(object = "GPC")}: ... }
    \item{getDLml}{\code{signature(object = "GPC")}: ... }
    \item{getHP}{\code{signature(object = "GPC")}: ... }
    \item{getK}{\code{signature(object = "GPC")}: ... }
    \item{getLml}{\code{signature(object = "GPC")}: ... }
    \item{initialize}{\code{signature(.Object = "GPC")}: ... }
    \item{predict}{\code{signature(object = "GPC")}: ... }
    \item{tuneHPs}{\code{signature(object = "GPC")}: ... }
    \item{update<-}{\code{signature(object = "GPC", value = "list")}: ... }
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
showClass("GPC")
}
\keyword{classes}
\usage{
GPC(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
