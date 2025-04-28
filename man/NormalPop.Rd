\name{NormalPop}
\Rdversion{1.1}
\docType{class}
\alias{NormalPop-class}
\alias{NormalPop}

\title{Class \code{"NormalPop"}}
\description{

  A normal distribution for initial values.

}

\section{Public Fields}{

  \describe{
    \item{mu}{Population mean.}
    \item{sigma}{Population standard deviation.}
  }
}
\section{Active Fields}{
  \describe{
    \item{pvec}{The vector \code{c(mu,log(sigma))}.}
  }
}

\section{Methods}{

  \describe{
    \item{initialize}{\code{signature(name,mu,sigma)}:  Constructor,
    called by \code{$new}.}
  }
}
\section{Constructor}{

The usual expression for building this distribution is
\code{NormalPop$new(\var{name},\var{mean},\var{sd})}.

}
\section{Superclasses and Subclasses}{

  Is a subclass of
  \code{\linkS4class{PopulationModel}}

}
\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{PopulationModel}}.

}

\examples{
apop <- NormalPop$new("standard",0,1)
apop$pvec
apop$draw(3)
apop$lprob(c(0,0),c(-1,0,1),c(1,2,1)/4)
apop$cdf(c(-Inf,-1,0,1,Inf))

}
\keyword{classes}
