\name{EvidenceModel}
\Rdversion{1.1}
\docType{class}
\alias{EvidenceModel-class}
\alias{EvidenceModel}

\title{Class \code{"EvidenceModel"}}
\description{

  This abstract class represents the likelihood of an observation given
  the value of the latent variable.

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
    \item{$draw}{\code{signature(theta,context=list())}:  Draws a random
      response.}
    \item{$llike}{\code{signature(Y,theta,context=list())}: Calculates the log
      likelihood of the observations given the latent variable and
      current parameters.}
    \item{$lprob}{\code{signature(par=self$pvec,Y,theta, weights,
	context=list())}: Calculated the log likelihood of the
      observations given estimated theta and argument parameters.}
    \item{$mstep}{\code{signature(Y,theta,weights,context=list(),
	its=3,control=list())}:  Runs the optimizer for \code{its} steps.}
    \item{$print}{\code{signature(...)}: Prints the object.}
    \item{$toString}{\code{signature(...)}: Creates a string
      represenation of the object.}
  }

  The arguments to these functions are as follows:
  \describe{
    \item{theta}{The value of the latent variable, possibly a vector
      over particles.}
    \item{context}{A list of additional measurement-occassion-subject
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

  This is an \code{\link[R6]{R6Class}}

  Generally, a subclass needs to implment the \code{$pvec} active field
  and the \code{$draw}, \code{$lprob}, \code{$llike}, \code{$toString} and
  \code{$initialize} methods.

  Currently implemented subclasses include:
  \code{\linkS4class{GradedResponse}}, \code{\linkS4class{NormalScore}}

}

\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{GradedResponse}}, \code{\linkS4class{NormalScore}},
  \code{\linkS4class{HMM}}, \code{\linkS4class{IRTF}}

}
\examples{
MyNormalScore <- R6Class(
    classname="MyNormalScore",
    inherit=EvidenceModel,
    public=list(
        initialize = function(name, bias, se) {
          self$name <- name
          self$bias <- bias
          self$se <- se
        },
        bias=0,
        se=01,
        draw = function(theta,context=list()) {
          rnorm(length(theta),theta+self$bias,self$se)
        },
        llike = function(Y,theta,context=list()) {
          dnorm(Y,theta+self$bias,self$se,log=TRUE)
        },
        lprob = function(par=self$pvec,Y,theta,weights,context=list()) {
          sum(dnorm(Y,theta+par[1],exp(par[2]),log=TRUE)*weights)
        },
        toString=function(digits=2,...){
          paste0("<NS: ", self$name, " ( ",
                 round(self$bias,digits=digits),
                 ", ",paste(round(self$se,digits=digits),
                            collapse = ", "), " )>")
        }

    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$bias,log(self$se)))
          self$bias <- value[1]
          self$se <- exp(value[2])
        }
    )
)


}
\keyword{classes}
