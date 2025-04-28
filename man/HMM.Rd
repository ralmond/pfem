\name{HMM}
\Rdversion{1.1}
\docType{class}
\alias{HMM-class}
\alias{HMM}

\title{Class \code{"HMM"} -- Hiddent Markov Model}
\description{

  This object contains both the description of the Hidden Markov Model,
  and after the particle filter or simultor is run contains the
  information from the procedure.


}
\section{Public Fields}{

  \describe{
    \item{name}{A character string identifying the model.}
    \item{popModels}{A list of \code{\linkS4class{PopulationModel}}s.}
    \item{groups}{A vector of length \code{nsubjects} giving indexes of
      which population model is used for each subject.  If all subjects
      use the same population model, this can be a vector of length one.}
    \item{growthModels}{A list of \code{\linkS4class{GrowthModel}}s.}
    \item{actions}{A \code{nsubjects} by \code{maxocc} matrix of indexes
      into the \code{growthModels} indicating which model is used
      by each subject on each occasion.  The matrix can have one row if
      all subjects use the same action, and/or one column if all
      occassions have the same growth model. }
    \item{evidenceModels}{A list of \code{\linkS4class{EvidenceModel}}s.}
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
      \code{\link{particleFilter}}.  See \code{$weights}.}
    \item{cltype}{The type argument of
      \code{\link[parallel]{makeCluster}}.  The value \code{"FORK"}
      is preferred for *nix, if you have windows your only choice is "PSOCK".}
    \item{clspec}{The first argument to
      \code{\link[parallel]{makeCluster}}, which is usually the number
      of cores to use.}
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
      evidence model to use on \code{subj} on occasion \code{it}.}
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

  To create a new HMM, use the expression
  \code{HMM$new(\var{name},\var{popModels},\var{growthModels},
  	\var{evidenceModels})}.

  The \code{popModels} argument should be a list of
  \code{\linkS4class{PopulationModel}}.  If there is more than one population,
  then the \code{$groups} field should be set to a vector of
  indexes for the population models of the same length as the number
  of subjects.

  The \code{growthModels} argument should be a list of
  \code{\linkS4class{GrowthModel}}.  If there is more than one growth model,
  then the \code{$actions} field should be set to a \code{$nsubjects}
  by \code{$maxocc} matrix of indexes into the growth models.  If all
  subjects have the same sequence of actions, then a matrix with one
  row can be used instead.  If all subjects use the same growth model
  for each occasion, then the the matrix can have a single column.

  The \code{evidenceModels} argument should be a lit of
  \code{\linkS4class{EvidenceModel}}.  If there is more than one evidence
  model, then the \code{$tasks} field should be set to a
  \code{$nsubjects} by \code{$maxocc} matrix of indexes into the
  evidence models.  If all subjects receive the same set of tasks,
  then the matrix can have a single row.  If the tasks are the same for
  every time period, then it can have a single column.

  Thus, in general, the \code{$groups}, \code{$actions}, and
  \code{$tasks} fields must be set after construction.  In addition,
  the observation times must be set by setting either \code{$deltaT} or
  \code{$times}.  Finally, for the particle filter, the
  \code{$data} field should be set as well.

}
\section{Representing Time}{

  To represent a continuous time process, there are two kinds of time.
  The first is a discrete index of the measurement occasions.  This
  index starts at 0, representing the initial value of \eqn{theta}$
  before any observations take place.  The second is the time of each
  observation.  It is probably best practice to represent these as elapsed
  time since Time 0, as well as ensuring they are all in the same units.

  Time for a single subject can be represented in one of two ways:
  (\code{$times}) as a vector of observation times or
  (\code{$deltaT}) as a vector of differences between observation times.
  Setting one of the two fields will automatically set the other.  These
  are stored in a matrix with rows corresponding to subjects and columns
  to measurement occastions.

  The \code{$times} matrix always has one more column than the
  \code{$deltaT} matrix.  The number of columns in \code{$deltaT},
  \code{$data}, \code{$actions} and \code{$tasks} should align.
}
\section{Missing and Ragged Data}{

The treatment of missing data in the particle filter depends on
which kinds of data are missing.

Missing data in the \code{$data} and \code{$task} fields indicates
that no observation was made at that opportunity.  A value for \eqn{\theta}
is imputed, but there is no update to the likelihood.

Missing data in the \code{$deltaT} or \code{$action} fields will stop
the filter for that observation.  No estimates will be made for that individual
for later measurement occasions.

Note that in the special case that \code{$deltaT} is zero means that
multiple observations were made at the same time point.  The growthModel is
not run.

Missing values elsewhere are unlikely to produce interpretable results.

}
\section{Parallel Computations}{

The \code{pfem} package uses the
\code{\link[parallel]{parallel-package}}. In particular,
looping over individuals is done using
\code{\link[parallel]{parLapply}} or one of its variants.  First
a cluster of worker processes is created and the filter for each
individual is assigned to one of the workers.

The \code{\link{particleFilter}} and \code{\link{simulate.HMM}} functions
call the \code{\link[parallel]{makeCluster}} with the arguments
\code{($clspec,$cltype,!!!$clargs)}.  See the docs for
\code{\link[parallel]{makeCluster}} for running across multiple machines.

In the simplest use case,
\code{$clspec} is the number of cores to use.  I recomend
using \code{parallel::detectCores()-2} (leaving one core for RStudio
and one for your operating system).  If unspecified it defaults to
\code{\link[base]{options}("mc.cores")} or 2 if that option is unspecified.

The cluster is automatically closed when the filter/simulator exits.
If \code{$stopClusterOnError=FALSE}, then the cluster will be left
open if the function encounters an error.

Finally, as sometimes an error will happen inside the parallelized function
both \code{\link{particleFilter}} and \code{\link{simulate.HMM}} have
a \code{debug} optional argument.  If that is set to true, then
\code{\link[base]{lapply}} will be used instead of
\code{\link[parallel]{parLapply}} to facilitate debugging.

}
\section{Superclasses and Subclasses}{

  This is an instance of \code{\link[R6]{R6Class}}.

}
\author{Russell G. Almond}
\seealso{

  \code{\linkS4class{PopulationModel}}, \code{\linkS4class{GrowthModel}},
  \code{\linkS4class{EvidenceModel}},
  \code{\link{particleFilter}}, \code{\link{simulate.HMM}},
  \code{\link{longResults}}

}

\examples{

make_pops <- function() {
  list(NormalPop$new("P1",0,1),
       NormalPop$new("P2",1,1))
}

make_growths <- function() {
  list(BrownianGrowth$new("G1",1,1),
       BrownianGrowth$new("G2",2,1))
}

make_ems <- function() {
  list(NormalScore$new("I1",0,.5),
       NormalScore$new("I2",.25,.25))
}

hmm<-HMM$new("stimulator",
             list(NormalPop$new("P1",-1/2,1),
                  NormalPop$new("P2",1/2,1)),
             list(BrownianGrowth$new("G1",1,1/100)),
             list(NormalScore$new("N0",0,1)))
hmm$clspec <- 1
hmm$groups <- 1:2
hmm$deltaT<-rep(1,3)
hmm$data <- outer(c(-1,1)/2,1:3,"+")
particleFilter(hmm,5)
longResults(hmm)

}
\keyword{classes}
