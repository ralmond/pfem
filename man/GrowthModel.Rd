\name{GrowthModel}
\Rdversion{1.1}
\docType{class}
\alias{GrowthModel-class}
\alias{GrowthModel}

\title{Class \code{"GrowthModel"}}
\description{

  An abstract class representing the distribution of the latent variable
  given its value at the previous occasion.

}
\section{Public Fields}{

  \describe{
    \item{name}{Identifier for the model.}
    \item{convergence}{If the maximizer was run, did it converge.}
    \item{lp}{If the maximizer was run, the log posterior.}
  }
}
\section{Active Fields}{
  \describe{
    \item{pvec}{Returns a vectorized version of the parameters.}
  }
}

\section{Methods}{

  \describe{
    \item{$draw}{\code{signature(theta,deltaT,variants=list())}:
      Draws a random response.}
    \item{$lprob}{\code{signature(par=self$pvec,theta0,theta1,weights,
	deltaT, variants=list())}: Calculated the log likelihood of the
      estimated growth given argument parameters.}
    \item{$mstep}{\code{signature(theta0,theta1,weights,deltaT,
	variants=list(), its=3,control=list())}:  Runs the optimizer for
      \code{its} steps.}
    \item{$print}{\code{signature(...)}: Prints the object.}
    \item{$toString}{\code{signature(...)}: Creates a string
      represenation of the object.}
  }

  The arguments to these functions are as follows:
  \describe{
    \item{theta}{The value of the latent variable at the previous
      occasion, vectorized over particles.}
    \item{theta0,theta1}{The value of the latent variable at the
      previous and current occasions, possibly a vector
      over particles.}
    \item{deltaT}{The elapsed time between the previous and current
      occasion.}
    \item{variants}{A list of additional measurement-occassion-subject
      specific covariates.}
    \item{Y}{The observed value.}
    \item{par}{The vectorized parameter.}
    \item{weights}{The weights for the particles.}
    \item{its}{The number of iterations of the optimizer to run.}
    \item{control}{Other control parameters passed to
      \code{\link{optim}}.}
  }

}
\section{Superclasses and Subclasses}{

  Its superclass is \code{\link[R6]{R6Class}}.

  Generally, a subclass needs to implment the \code{$pvec} active field
  and the \code{$draw}, \code{$lprob}, \code{$toString} and
  \code{$initialize} methods.

  Existing implementations are:
  \code{\linkS4class{BrownianGrowth}}, \code{\linkS4class{SpurtGrowth}}

}
\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{HMM}}
  \code{\linkS4class{BrownianGrowth}}, \code{\linkS4class{SpurtGrowth}}

}
\examples{
MyBrownianGrowth <- R6Class(
    classname="MyBrownianGrowth",
    inherit=GrowthModel,
    public=list(
        initialize = function(name,gain,inovSD) {
          self$name <- name
          self$gain <- gain
          self$inovSD <- inovSD
        },
        gain=0,
        inovSD=.1,
        draw = function(theta,deltaT,variants=list()){
          rnorm(length(theta),theta+self$gain*deltaT,
                self$inovSD*sqrt(deltaT))
        },
        lprob = function(par=self$pvec,theta0,theta1,weights,deltaT,
                         variants=list()) {
          mu <- theta0+par[1]*deltaT
          sig <- exp(par[2])*sqrt(deltaT)
          sum(dnorm(theta1,mu,sig,log=TRUE)*weights)
        },
        toString=function(digits=2,...){
          paste0("<BrownianGrowth: ", self$name, " ( ",
                 round(self$gain,digits=digits),
                 ", ",round(self$inovSD,digits=digits)," )>")
        }

    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$gain,log(self$inovSD)))
          self$gain <- value[1]
          self$inovSD <- exp(value[2])
        }
    )
)

}
\keyword{classes}
