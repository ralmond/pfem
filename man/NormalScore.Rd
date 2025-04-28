\name{NormalScore}
\Rdversion{1.1}
\docType{class}
\alias{NormalScore-class}
\alias{NormalScore}

\title{Class \code{"NormalScore"}}
\description{

  An evidence model where the observation has a normal error from the
  latent variable.  It was developed for testing, but may be useful for
  "true score" type models

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

Generally built by calling
\code{NormalScore$new(\var{name},\var{bias},\var{se})}.

}
\section{Superclasses and Subclasses}{

  Its superclass is:
  \code{\linkS4class{EvidenceModel}}

}
\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{EvidenceModel}}

}

\examples{
em <- NormalScore$new("ns",1,1/25)
theta <- -2:2
em$draw(theta)
em$llike(-1:1,theta)
em$lprob(c(0,0),0, theta,(1:5)/15)
}
\keyword{classes}
