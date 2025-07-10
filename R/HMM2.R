## Theta dimensions -- time,nDims,particle,subject
## Data dimensions -- time,nDat,subject
## Weight dimension -- particle,subject
## Weight dimension -- particle,subject
## actions (it, subj)
## tasks(it,subj)
## time (it,subj)

HMM2 <- R6Class(
  "HMM",
  inherits=HMM,
  public=list(
    name="<HMM2>",
    data=matrix(NA_integer_,1L,2L,1L),
    thetaNames=c("theta1","theta2"),
    dataNames=c("Obs1","Obs2"),
    npart=1L,
    theta = numeric(),
    lweights = numeric(),
    initialize=function(name,thetaNames,dataNames,popModels,growthModels,evidenceModels) {
      super$initialize(name,popModels,growthModels,evidenceModels)
      self$thetaNames <- thetaNames
      self$dataNames <- dataNames
      }
  ),
  active=list(
    nDims=function() {length(self$thetaNames)}
    nDat=function(){length(self$dataNames)}
  )
)

particleFilter.HMM2 <- function (hmm, npart=hmm$npart, seed=NULL,
                                debug=FALSE, weightLog=FALSE) {
  ## Setup Clusters
  plapply <- lapply
  if (!debug) {
    clust <- inject(makeCluster(hmm$clspec,hmm$cltype,!!!hmm$clargs))
    stopOnExit <- hmm$stopClusterOnError
    withr::defer({if (stopOnExit) stopCluster(clust)})
    plapply <- function(...) parLapply(clust,...)
    if (!missing(seed)) {
      clusterSetRNGStream(clust,seed)
      mc.reset.stream()
    }
  }


  hmm$npart <- npart
  hmm$lweights <- matrix(0,npart,hmm$nsubjects)
  ## <<Put persons last???>>
  hmm$theta <- array(NA_real_,c(hmm$nDims,npart,hmm$nsubjects,hmm$maxocc+1L))
  hmm$weightLog <- list()

  filtres <- plapply(1L:hmm$nsubjects, \(subj) {
    stheta <- array(NA_real_,c(hmm$nDims,npart,hmm$maxocc+1L))
    sweights <- rep(0,npart)
    stheta[,,1L] <- hmm$drawPop(subj,npart)
    swl <- list()

    for (it in 1L:hmm$maxocc) {
      stheta[it+1L,,] <- hmm$drawGrowth(it,subj,stheta[it,,])

      sweights <- hmm$lweights + hmm$evalEvidence(it,subj,stheta[it+1L,,],
                                                  hmm$data[it,,subj])
      if(isTRUE(weightLog))
        swl <- c(swl,sweights)
    }
    list(theta=stheta,weights=sweights,wl=swl)
  })
  for (subj in 1L:hmm$nsubjects) {
    hmm$theta[,,,subj] <- filtres[[subj]]$theta
    hmm$lweights[,subj] <- filtres[[subj]]$weights
    hmm$weightLog <- c(hmm$weightLog,filtres[[subj]]$wl)
  }
  stopOnExit <- TRUE
  invisible(hmm)
}


# Can reuse stats::simulate
#simulate <- function(object,nsim=hmm$nsubjects,seed=NULL,...,mocc=hmm$maxocc)
#  UseMethod("simulate")

simulate.HMM2 <- function(object,nsim=object$nsubjects,seed=NULL,mocc=2L,...,
                         debug=FALSE) {
  plapply <- lapply
  if (!isTRUE(debug)) {
    clust <- inject(makeCluster(object$clspec,object$cltype,!!!object$clargs))
    plapply <- function(...) parLapply(clust,...)
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
  object$theta <- array(NA_real_,c(object$nDims,object$npart,nsim,mocc+1L))
  object$data <- array(NA_real_,c(object$nDat,nsim,mocc))
  simres <- plapply(1L:object$nsubjects, \(subj) {
    stheta <- array(NA_real_,c(object$nDims,mocc+1L))
    sdata <- array(NA_real_,c(object$nDat,mocc))

    stheta[,1L] object$drawPop(subj,object$npart)

    for (it in 1L:mocc) {
      stheta[,it+1L] <- object$drawGrowth(subj,it)                                                         object$delT(subj,it))

      task <- object$task(subj,it)
      if (!is.na(task)) {
       sdata[,it] <- object$evidenceModels[[task]]$draw(stheta[,it+1L])
      }
    }
    list(theta=stheta,data=sdata)
  })

  for (isim in 1:nsim) {
    object$theta[,1L,isim,] <- simres[[isim]]$theta
    object$data[,isim,] <- simres[[isim]]$data
  }




  stopOnExit <- TRUE
  invisible(object)
}

