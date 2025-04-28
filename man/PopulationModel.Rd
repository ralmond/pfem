\name{PopulationModel}
\Rdversion{1.1}
\docType{class}
\alias{PopulationModel-class}
\alias{PopulationModel}

\title{Class \code{"PopulationModel"}}
\description{

  An abstract class describing the initial distribution of a population.

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
    \item{$draw}{\code{signature(theta,invariants=list())}:
      Draws a random starting position.}
    \item{$cdf}{\code{signature(theta,invariants=list())}:
      Calculates the cumulative distribution function.}
    \item{$lprob}{\code{signature(par=self$pvec,theta,weights,
	invariants=list())}: Calculated the log likelihood of the
      estimated growth given argument parameters.}
    \item{$mstep}{\code{signature(theta,weights,invariants=list(),
	its=3,control=list())}:  Runs the optimizer for
      \code{its} steps.}
    \item{$print}{\code{signature(...)}: Prints the object.}
    \item{$toString}{\code{signature(...)}: Creates a string
      represenation of the object.}
  }

  The arguments to these functions are as follows:
  \describe{
    \item{theta}{The value of the latent variable at the initial time,
      vectorized over particles.}
    \item{invariants}{A list of additional subject-specific time
      invariant covariates.}
    \item{par}{The vectorized parameter.}
    \item{weights}{The weights for the particles.}
    \item{its}{The number of iterations of the optimizer to run.}
    \item{control}{Other control parameters passed to
      \code{\link{optim}}.}
  }
}
\section{Constructor}{

  This is an abstract class, usually just used for subclasses.

}
\section{Superclasses and Subclasses}{

  This is an instance of \code{\link[R6]{R6Class}}.

  Generally, a subclass needs to implment the \code{$pvec} active field
  and the \code{$draw}, \code{$lprob}, \code{$cdf}, \code{$toString} and
  \code{$initialize} methods.

  Existing subclasses:  \code{\linkS4class{NormalPop}}.

}
\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{NormalPop}}.
  \code{\linkS4class{HMM}}, \code{\linkS4class{IRTF}}

}
\examples{

MyNormalPop <- R6Class(
    classname = "MyNormalPop",
    inherit=PopulationModel,
    public=list(
        initialize = function(name,mu,sigma) {
          self$name <- name
          self$mu <- mu
          self$sigma <- sigma
        },
        mu=0,
        sigma=1,
        draw = function(npart,invariants=list()) {
          mu <- self$mu
          sigma <- self$sigma
          rnorm(npart,mu,sigma)
        },
        lprob = function(par=self$pvec,theta,weights,invariants=list()) {
          mu <- par[1]
          sigma <- exp(par[2])
          sum(dnorm(theta,mu,sigma,log=TRUE)*weights)
        },
        cdf = function(theta,invariants=list()) {
          pnorm(theta,self$mu,self$sigma)
        },
        toString=function(digits=2,...){
          paste0("<NormalPopulation: ",
          self$name, " ( ",round(self$mu,digits=digits),
                 ", ",round(self$sigma,digits=digits)," )>")
        }
    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$mu,log(self$sigma)))
          self$mu <- value[1]
          self$sigma <- exp(value[2])
        }
    )
)

}
\keyword{classes}
