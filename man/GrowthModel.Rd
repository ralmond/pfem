\name{GrowthModel}
\Rdversion{1.1}
\docType{class}
\alias{GrowthModel-class}

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
      \code{\it} steps.} 
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
\section{Constructor}{

}
\section{Superclasses and Subclasses}{

  \code{\link{BrownianGrowth}}, \code{\link{SpurtGrowth}}

}
\author{Russell G. Almond}
\seealso{

  \code{\link{BrownianGrowth}}, \code{\link{SpurtGrowth}}, 

}
\examples{

}
\keyword{classes}
