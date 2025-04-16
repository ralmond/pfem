\name{NormalScore}
\Rdversion{1.1}
\docType{class}
\alias{NormalScore-class}

\title{Class \code{"NormalScore"}}
\description{

  An evidence model where the observation has a normal error from the
  latent variable.
  
}
\section{Public Fields}{

  \describe{
    \item{bias}{An offset added to theta to get the mean.}
    \item{se}{The standard deviation of the error.}
  }
}
\section{Active Fields}{
  \describe{
    \item{pvec}{The vector \code{c(bias,log(se))}.}
  }
}

\section{Methods}{

  \describe{
    \item{initialize}{\code{signature(name,bias,se)}:  The constructor,
      called from \code{$new}.}
  }
}
\section{Constructor}{

}
\section{Superclasses and Subclasses}{

  \code{\link{EvidenceModel}}

}
\author{Russell G. Almond}
\seealso{

  \code{\link{EvidenceModel}}

}

\examples{

}
\keyword{classes}
