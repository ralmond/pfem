\name{GradedResponse}
\Rdversion{1.1}
\docType{class}
\alias{GradedResponse-class}

\title{Class \code{"GradedResponse"}}
\description{

  An implementation of Samejimia's (1969) graded response model.  Two
  parameters, \code{a}, a discrimination and \code{b} a vector of
  difficulties, where
  \deqn{Pr(Y\ge k|\theta, a, b) = \mathrm{logit}^{-1}(a\theta -b[k])}.
  So the possible values for Y are integers from 0 to \code{length(b)}.

}
\section{Public Fields}{

  \describe{
    \item{a}{A discrimination parameter, should be positive.}
    \item{b}{A vector of difficulties, should be increasing.}
  }
}
\section{Active Fields}{
  \describe{
    \item{pvec}{The parameter vector is \code{c(log(a),b)}.}
  }
}

\section{Methods}{

  \describe{
    \item{cuts}{\code{signature(theta,a=self$a,b=self$b)}:
      Generates the cut points (cumulative probabilities).}
  }
}
\section{Constructor}{

}
\section{Superclasses and Subclasses}{

  Instance of \code{\link{EvidenceModel}}

}
\author{Russell G. Almond}
\seealso{

  \code{\link{EvidenceModel}}.

}
\references{

  Samejima, F. (1969). Estimation of latent ability using a response
  pattern of graded scores. _Psychometrika Monograph No. 17_, 34(4), (Part
  2).
  \url{https://www.psychometricsociety.org/sites/main/files/file-attachments/mn17.pdf}

}
\examples{

gr <- GradedResponse$new("agr",1,c(-1,1))
gr
gr$pvec
gr$draw(c(-1,0,1))
gr$llike(c(0,1,0),c(-1,0,1))
gr$lprob(gr$pvec,c(0,1,0),c(-1,0,1),c(.4,.3,.3))




}
\keyword{classes}
