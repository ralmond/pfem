\name{PopulationModel}
\Rdversion{1.1}
\docType{class}
\alias{PopulationModel-class}

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
      \code{\it} steps.} 
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

}
\section{Superclasses and Subclasses}{

  \code{\link{NormalPop}}.

}
\author{Russell G. Almond}
\seealso{

  \code{\link{NormalPop}}.

}

\examples{

}
\keyword{classes}
