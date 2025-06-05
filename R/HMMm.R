
HMMm <- R6Class(
    "HMMm",
    inherit=HMM,
    public=list(
      name="<HMMm>",
      data=matrix(NA_integer_,1L,1L,1L),
      thetaNames="theta",
      theta = numeric(),
      lweights = numeric(),
      toString=function(...) {
        paste0("<HMMm: ",self$name,": ",
               self$nsubjects, " x ",
               self$macocc, " >")
      },
      drawGrowth = function (subj,it) {
        if (self$delT(subj,it)==0) self$theta[,,subj,it]
        else
          self$growthModels[[self$action(subj,it)]]$draw(self$theta[,,subj,it],
                                                       self$delT(subj,it))
      },
      evalEvidence= function (subj,it) {
        task <- self$task(subj,it)
        Y <- self$data[,subj,it]
        if (any(is.na(Y)) || is.na(task)) return(rep(0,dim(self$theta)[1:2]))
        else {
          self$evidenceModels[[task]]$llike(Y,self$theta[,,subj,it+1L])
        }
      }
    ),
    private=list(
        DeltaT = numeric(),
        Times = numeric(),
        Nsubjects = 1,
        Maxtimes = 1
    ),
    active=list(
        thetadim=function(){length(self$thetaNames)},
        datadim=function(){dim(self$data)[1]},
        nsubjects = function() {dim(self$data)[2]},
        maxocc = function(value) {dim(self$data)[3]},
    )
)


particleFilter.HMMm <- function (hmm, npart=hmm$npart, seed=NULL,
                                debug=FALSE, weightLog=FALSE) {
  ## Setup Clusters
  psapply <- sapply
  if (!debug) {
    clust <- inject(makeCluster(hmm$clspec,hmm$cltype,!!!hmm$clargs))
    stopOnExit <- hmm$stopClusterOnError
    withr::defer({if (stopOnExit) stopCluster(clust)})
    psapply <- function(...) parSapply(clust,...)
    if (!missing(seed)) {
      clusterSetRNGStream(clust,seed)
      mc.reset.stream()
    }
  }
  if(isTRUE(weightLog)) hmm$weightLog <- list()

  hmm$npart <- npart
  hmm$lweights <- matrix(0,npart,hmm$nsubjects)
  hmm$theta <- array(NA_real_,c(npart,hmm$thetadim,hmm$nsubjects,hmm$maxocc+1L))
  hmm$theta[,,,1L] <- psapply(1L:hmm$nsubjects, \(subj) {
    hmm$drawPop(subj,npart)
  })

  for (it in 1L:hmm$maxocc) {
    if (debug) cat("Time: ",it,".\n")
    hmm$theta[,,,it+1L] <- psapply(1L:hmm$nsubjects, \(subj) {
      hmm$drawGrowth(subj,it)
    })

    hmm$lweights <- hmm$lweights +
      psapply(1L:hmm$nsubjects, \(subj) {
        hmm$evalEvidence(subj,it)
      })
    if(isTRUE(weightLog))
      hmm$weightLog <- c(hmm$weightLog,hmm$lweights)
  }
  stopOnExit <- TRUE
  invisible(hmm)
}


# Can reuse stats::simulate
#simulate <- function(object,nsim=hmm$nsubjects,seed=NULL,...,mocc=hmm$maxocc)
#  UseMethod("simulate")

simulate.HMMm <- function(object,nsim=object$nsubjects,seed=NULL,mocc=2L,dataDim=1L,...,
                         debug=FALSE) {
  psapply <- sapply
  if (!isTRUE(debug)) {
    clust <- inject(makeCluster(object$clspec,object$cltype,!!!object$clargs))
    psapply <- function(...) parSapply(clust,...)
    stopOnExit <- object$stopClusterOnError
    withr::defer({if (stopOnExit) stopCluster(clust)})
    if (!missing(seed)) {
        clusterSetRNGStream(clust,seed)
        mc.reset.stream()
    }
  }
  if (missing(mocc)) {
    if(length(object$deltaT)>1L)
      mocc <- length(object$deltaT)
  }

  object$npart <- 1L
  object$theta <- array(NA_real_,c(object$npart,object$thetadim,nsim,mocc+1L))
  object$data <- array(NA_integer_,c(dataDim,nsim,mocc))
  object$theta[,,,1L] <- psapply(1L:object$nsubjects, \(subj) {
    object$popModels[[object$group(subj)]]$draw(object$npart)
  })


  for (it in 1L:mocc) {
    if (debug) print("Time: ",it)
    object$theta[,,,it+1L] <- psapply(1L:object$nsubjects, \(subj) {
      object$growthModels[[object$action(subj,it)]]$draw(object$theta[,subj,it],
                                           object$delT(subj,it))
    })

    object$data[,,it] <- psapply(1L:object$nsubjects, \(subj) {
      task <- object$task(subj,it)
      if (is.na(task)) return(NA_integer_)
      else {
        object$evidenceModels[[task]]$draw(object$theta[,subj,it+1L])
      }
    })

  }
  stopOnExit <- TRUE
  invisible(object)
}

