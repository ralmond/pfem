\name{SpurtGrowth}
\Rdversion{1.1}
\docType{class}
\alias{SpurtGrowth-model}

\title{Class \code{"SpurtGrowth"}}
\description{

This growth model randomly switches between two Brownian Growth models,
one with a slow growth and one with a fast growth with a certain probability.

}
\section{Public Fields}{

  \describe{
    \item{name}{Identifier for the model.}
    \item{gain0}{The lower gain rate.}
    \item{gain1}{The higher gain rate.}
    \item{p}{The probability of the higher gain.}
    \item{inovSD}{The standard deviation of the innovation.}
  }
}
\section{Active Fields}{
  \describe{
    \item{pvec}{The parameter vector,
      \code{c(gain0,gain1,logit(p),log(inovSD))}.}
  }
}

\section{Methods}{

  \describe{
    \item{initialize}{\code{signature(name, gain0, gain1, p inovSD)}:
      The basic constructors.}
  }
}
\section{Constructor}{

}
\section{Superclasses and Subclasses}{

  \code{\link{GrowthModel}}
}
\author{Russell G. Almond}
\seealso{

  \code{\link{GrowthModel}}.
  \code{\link{BrownianGrowth}}.

}
\examples{

}
\keyword{classes}
