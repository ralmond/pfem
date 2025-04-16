\name{NormalPop}
\Rdversion{1.1}
\docType{class}
\alias{NormalPop-class}

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

}
\section{Superclasses and Subclasses}{

  \code{\link{PopulationModel}}

}
\author{Russell G. Almond}
\seealso{

  \code{\link{}}.

}

\examples{

}
\keyword{classes}
