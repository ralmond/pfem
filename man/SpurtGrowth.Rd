\name{SpurtGrowth}
\Rdversion{1.1}
\docType{class}
\alias{SpurtGrowth-class}
\alias{SpurtGrowth}

\title{Class \code{"SpurtGrowth"}}
\description{

This growth model randomly switches between two Brownian Growth models,
one with a slow growth and one with a fast growth with a certain probability.

This is a mixture distribution:
  \deqn{F(\theta_1|\theta_0) = p N(\theta_0+b_1\Delta T, \sigma\sqrt{\Delta T})
  + (1-p) N(\theta_0+b_0\Delta T, \sigma\sqrt{\Delta T}),} where \eqn{\sigma} is
  \code{inovSD}, and \eqn{b_0} and \eqn{b_1} are \code{gain0} and \code{gain1}.


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
    \item{initialize}{\code{signature(name, gain0, gain1, p, inovSD)}:
      The basic constructors.}
  }
}
\section{Constructor}{

  The recommended way to build this growth model is
  \code{SpurtGrowth$new(\var{name},\var{gain0},\var{gain1},\var{p},\var{inovSD})}.


}
\section{Superclasses and Subclasses}{

  This is a subclass of
  \code{\linkS4class{GrowthModel}}
}
\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{GrowthModel}}.
  \code{\linkS4class{BrownianGrowth}}.

}
\examples{

agr <- SpurtGrowth$new("S",0,1,.5,1)
agr$pvec
agr$draw(c(-1,0,1),1)
theta0 <- rep(0,3)
theta1 <- theta0+1:3
weights <- (1:3)/9
agr$lprob(c(0,1,-Inf,0),theta0,theta1,weights,1)

}
\keyword{classes}
