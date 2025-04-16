\name{BrownianGrowth}
\Rdversion{1.1}
\docType{class}
\alias{BrownianGrowth-class}

\title{Class \code{"BrownianGrowth"}}
\description{
  
  A non-stationary Brownian motion process, where the gain and the
  variance are linear with the elapsed time.

  \deqn{\theta(t+1) \sim N(\theta(t) + b \Delta t, \sigma \sqrt{\Delta
      t}),}
  where \eqn{b} is the \code{gain} and \eqn{\sigma} is the
  \code{inovSD}. 

}
\section{Public Fields}{

  \describe{
    \item{gain}{The change per unit time.}
    \item{inovSD}{The standard deviation per unit time.}
  }
}
\section{Active Fields}{
  \describe{
    \item{pvec}{The probablity vector, \code{c(gain,log(inovSD))}.
  }
}
\section{Methods}{

  \describe{
    \item{initialize}{\code{signature(name,gain,inovSD)}:  Called by
    \code{$new}.}
  }
}
\section{Constructor}{

}
\section{Superclasses and Subclasses}{

  Subclass of \code{\link{GrowthModel}}

}
\author{Russell G. Almond}
\seealso{

  \code{\link{GrowthModel}}.

}

\examples{

}
\keyword{classes}
