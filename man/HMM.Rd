\name{HMM}
\Rdversion{1.1}
\docType{class}
\alias{HMM-class}

\title{Class \code{"HMM"} -- Hiddent Markov Model}
\description{

  This object contains both the description of the Hidden Markov Model,
  and after the particle filter or simultor is run contains the
  information from the procedure.

  
}
\section{Public Fields}{

  \describe{
    \item{name}{A character string identifying the model.}
    \item{popModels}{A list of \code{\link{PopulationModel}}s.}
    \item{groups}{A vector of length \code{nsubjects} giving indexes of
      which population model is used for each subject.  If all subjects
      use the same population model, this can be a vector of length one.}
    \item{growthModels}{A list of \code{\link{GrowthModel}}s.}
    \item{actions}{A \code{nsubjects} by \code{maxocc} matrix of indexes
      into the \code{growthModels} indicating which model is used
      by each subject on each occasion.  The matrix can have one row if
      all subjects use the same action, and/or one column if all
      occassions have the same growth model. }
    \item{evidenceModels}{A list of \code{\link{EvidenceModel}}s.}
    \item{tasks}{A \code{nsubjects} by \code{maxocc} matrix of indexes
      into the \code{evidenceModels} indicating which evidence model is used
      by each subject on each occasion.  The matrix can have one row if
      all subjects use the same set of tasks, and/or one column if all
      occassions have the same evidence models. }
    \item{data}{A \code{nsubjects} by \code{maxocc} matrix of data
      values.  Note that the \code{$nsubjects} and \code{$maxocc} active
      fields are defined by the dimensions of this matrix.}
    \item{theta}{This is a \code{npart} by \code{nsubjects} by
      \code{maxocc} array of theta values.  This is generally populated
      via a call to \code{\link{simulate}}, or
      \code{\link{particleFilter}}.} 
    \item{thetaNames}{The name to use for \code{theta} in data tables.}
    \item{npart}{The number of particle to use.}
    \item{lweights}{A \code{npart} by \code{nsubject} matrix of
      unnormalized log weights.  Generally populated
      \code{\link{particleFilter}}.  See \code{$weights}
    \item{cltype}{The type argument of
      \code{\link[parallel]{makeCluster}}.  The value \code{"FORK"}
      is preferred for *nix, if you have windows your only choice is "PSOCK".}
    \item{clspec}{The first argument to
      \code{\link[parallel]{makeCluster}}, which is usually the number
      of cores to use.
    \item{clargs}{Additional arguments to \code{\link[parallel]{makeCluster}}}
    \item{stopClusterOnError}{A logicial value.  If \code{TRUE}, the
      cluster will be closed on an error.}
   }
}
\section{Active Fields}{
  \describe{
    \item{deltaT}{}
    \item{times}{}
    \item{nsubjects}{Returns the number of subjects base on the number
      of rows of the \code{$data} field.}
    \item{maxocc}{Returns the number of maximum number of occasions
      based on the number of columns of the \code{$data} field.}
    \item{weights}{A \code{npart} by \code{nsubject} matrix of weights
      on the probability scale, normalized so that the sum across
      particles is one.}
  }
}
\section{Public Methods}{

  \describe{
    \item{initialize}{\code{signature(name,popModels,growthModels,
	evidenceModels)}:  Called through \code{$new}.}
    \item{delT}{\code{signature(subj,it)}: Returns the time increment
      for a specific subject and occasion.}
    \item{group}{\code{signature(subj)}: Returns the index for the
      population model for \code{subj}. }
    \item{action}{\code{signature(subj,it)}:  Returns the index for the
      growth model to use on \code{subj} on occasion \code{it}.}
    \item{task}{\code{signature(subj,it)}:  Returns the index for the
      evidence  model to use on \code{subj} on occasion \code{it}.}
    \item{toString}{\code{signature(...)}:  Provides a string
      representation for printing.}
    \item{print}{\code{signature(...)}:  Prints a representation.}
  }
}
\section{S3 Methods}{

  \describe{
    \item{\code{\link{particleFilter}}}{\code{signature(hmm, npart=hmm$npart,
	  seed=NULL, debug=FALSE, weightLog=FALSE)}:  Runs the particle
	filter. }
    \item{\code{\link{simulate}}}{\code{signature(hmm, nsim=hmm$nsubjects,
	  seed=NULL,mocc=2L,...,debug=FALSE)}: Generates simulated data
	from the model.}
    \item{\code{\link{longResults}}}{\code{signature(hmm)}:  Extracts
      the result of a particle filter or simulation as a long form data
      frame.}
  }
}
\section{Constructor}{

}
\section{Representing Time}{

}
\section{Missing and Ragged Data}{

}
\section{Parallel Computations}{


}
\section{Superclasses and Subclasses}{


}
\author{Russell G. Almond}
\seealso{

  \code{\link{}}.

}

\examples{

}
\keyword{classes}
