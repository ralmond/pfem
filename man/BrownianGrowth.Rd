\name{BrownianGrowth}
\Rdversion{1.1}
\docType{class}
\alias{BrownianGrowth-class}
\alias{BrownianGrowth}

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
    \item{pvec}{The probablity vector, \code{c(gain,log(inovSD))}.}
  }
}
\section{Methods}{

  \describe{
    \item{initialize}{\code{signature(name,gain,inovSD)}:  Called by
    \code{$new}.}
  }
}
\section{Constructor}{

To create a new object of this type call:
\code{BrownianGrowth$new(\var{name},\var{gain},\var{inovSD})}
}
\section{Superclasses and Subclasses}{

  Subclass of \code{\linkS4class{GrowthModel}}

}
\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{GrowthModel}}, \code{\linkS4class{HMM}}

}

\examples{

bgm <- BrownianGrowth$new("amodel",.1,.05)
bgm
bgm$pvec
bgm$draw(c(-1,0,1),1)
bgm$lprob(bgm$pvec,c(-1,0,1),c(-.9,.05,.9),c(.4,.3,.3),1)


}
\keyword{classes}
