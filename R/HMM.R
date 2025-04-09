name <- function(obj) {UseMethod("name")}
name.R6 <- function(obj){obj$name}

HMM <- R6Class(
    "HMM",
    public=list(
      name="<HMM>",
        popModels=list(),
        groups = 1L,
        growthModels=list(),
        actions = matrix(1L,1L,1L),
        evidenceModels=list(),
        tasks=matrix(1L,1L,1L),
        data=matrix(NA_integer_,1L,1L),
        thetaNames="theta",
        npart=1L,
        theta = numeric(),
        lweights = numeric(),
        cltype=ifelse(.Platform$OS.type=="windows","PSOCK","FORK"),
        clspec=getOption("mc.cores",2L),
        clargs=list(),
        stopClusterOnError=TRUE,
        initialize=function(name,popModels,growthModels,evidenceModels) {
          self$name <- name
          self$popModels <- popModels
          if (!is.list(popModels)) self$popModels <- list(popModels)
          self$growthModels <- growthModels
          if (!is.list(growthModels)) self$growthModels <- list(growthModels)
          self$evidenceModels <- evidenceModels
          if (!is.list(evidenceModels)) self$evidenceModels <-
                                          list(evidenceModels)
        },
        delT = function(subj,it) {
          if (nrow(private$DeltaT)==1L) subj <- 1L
          if (ncol(private$DeltaT)==1L) it <- 1L
          private$DeltaT[subj,it]
        },
        group = function (subj) {
          if (length(self$groups)==1L) return(self$groups)
          self$groups[subj]
        },
        action = function(subj,it) {
          if (!is.matrix(self$actions)) {
            if (length(self$actions)==1L)
              self$action <- matrix(self$actions,1L,1L)
            else if (length(self$actions)==self$maxocc)
              self$action <- matrix(self$action,1L,self$maxocc)
            else if (length(self$actions)==self$nsubjects)
              self$action <- matrix(self$actions,self$nsubjects,1L)
            else stop ("Unexpected size for action matrix.  Should be ",
                       self$nsubjects, "x", self$maxocc,".")
          }
          if (nrow(self$actions)==1L) subj <- 1L
          if (ncol(self$actions)==1L) it <- 1L
          self$actions[subj,it]
        },
        task = function(subj,it) {
          if (!is.matrix(self$tasks)) {
            if (length(self$tasks)==1L)
              self$task <- matrix(self$tasks,1L,1L)
            else if (length(self$tasks)==self$maxocc)
              self$tasks <- matrix(self$task,1L,self$maxocc)
            else if (length(self$tasks)==self$nsubjects)
              self$tasks <- matrix(self$tasks,self$nsubjects,1L)
            else stop ("Unexpected size for task matrix.  Should be ",
                       self$nsubjects, "x", self$maxocc,".")
          }
          if (nrow(self$tasks)==1L) subj <- 1L
          if (ncol(self$tasks)==1L) it <- 1L
          self$tasks[subj,it]
        },
        toString=function(...) {
          paste0("<HMM: ",self$name,": ",
                 self$nsubjects, " x ",
                 self$macocc, " >")
        },
        print=function(...) {
          print(self$toString(...),...)
        }
    ),
    private=list(
        DeltaT = numeric(),
        Times = numeric(),
        Nsubjects = 1,
        Maxtimes = 1
    ),
    active=list(
        deltaT = function(value) {
          if (missing(value)) return(private$DeltaT)
          if (!is.matrix(value)) {
            value <- matrix(value,1L,length(value))
          }
          private$DeltaT <- value
          private$Times <- t(apply(cbind(0,value),1,cumsum))
        },
        times = function(value) {
          if (missing(value)) return(private$Times)
          if (!is.matrix(value))
            value <- matrix(value,1L,length(value))
          private$Times <- value
          private$DeltaT <- t(apply(value,1,diff))
        },
        nsubjects = function() {nrow(self$data)},
        maxocc = function(value) {ncol(self$data)},
        weights = function() {
          weights <- exp(self$lweights)
          if (length(weights)==0L)
            return(weights)
          if(!is.matrix(weights))
            return(weights/sum(weights))
          return(sweep(weights,2,colSums(weights),"/"))
        }
    )
)


particleFilter <- function (hmm, npart=hmm$npart) {
  UseMethod("particleFilter")
}


particleFilter.HMM <- function (hmm, npart=hmm$npart, seed=NULL,
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
  hmm$theta <- array(NA_real_,c(npart,hmm$nsubjects,hmm$maxocc+1L))
  hmm$theta[,,1L] <- psapply(1L:hmm$nsubjects, \(subj) {
    hmm$popModels[[hmm$group(subj)]]$draw(npart)
  })

  for (it in 1L:hmm$maxocc) {
    hmm$theta[,,it+1L] <- psapply(1L:hmm$nsubjects, \(subj) {
      hmm$growthModels[[hmm$action(subj,it)]]$draw(hmm$theta[,subj,it],
                                           hmm$delT(subj,it))
    })

    hmm$lweights <- hmm$lweights +
      psapply(1L:hmm$nsubjects, \(subj) {
        task <- hmm$task(subj,it)
        Y <- hmm$data[subj,it]
        if (is.na(Y) || is.na(task)) return(0)
        else {
          hmm$evidenceModels[[task]]$llike(Y,hmm$theta[,subj,it+1L])
        }
      })
    if(isTRUE(weightLog))
      hmm$weightLog <- c(hmm$weightLog,hmm$lweights)
  }
  stopOnExit <- TRUE
  invisible(hmm)
}

longResults <- function (hmm) UseMethod("longResults")

longResults.HMM <- function (hmm) {
  npp <- hmm$npart*hmm$nsubjects
  nrw <- npp*(hmm$maxocc+1L)
  ptheta <- as.vector(hmm$theta)

  Y <- rep(as.vector(cbind(NA,hmm$data)),each=hmm$npart)
  tasks <- hmm$tasks
  if (!is.matrix(tasks) || ncol(tasks) < hmm$maxocc)
    tasks <- matrix(as.vector(tasks),1L,hmm$maxocc)
  tasks <- cbind(NA,tasks)
  if (nrow(tasks) > 1L)
    tasks <- rep(as.vector(tasks),each=hmm$npart)
  else
    tasks <- rep(as.vector(tasks),each=npp)
  if (nrow(hmm$times) > 1L)
    alltimes <- rep(as.vector(hmm$times),each=hmm$npart)
  else
    alltimes <- rep(as.vector(hmm$times),each=npp)

  subj<-rep(rep(1:hmm$nsubjects,each=hmm$npart),hmm$maxocc+1L)
  occ<-rep(0L:hmm$maxocc,each=npp)
  weights <- as.vector(hmm$weights)
  if (length(hmm$weights)==0L)
    weights <- rep(NA,npp)
  weights <-rep(weights,hmm$maxocc+1L)

  result <-data.frame(
      subj=subj,
      occ=occ,
      time=alltimes,
      tasks=tasks,
      Y=Y,
      weights=weights,
      ptheta)
  names(result) <- c("subj","occ","time","tasks","Y","weights",hmm$thetaNames)
  result
}

# Can reuse stats::simulate
#simulate <- function(object,nsim=hmm$nsubjects,seed=NULL,...,mocc=hmm$maxocc)
#  UseMethod("simulate")

simulate.HMM <- function(hmm,nsim=hmm$nsubjects,seed=NULL,mocc=2L,...,
                         debug=FALSE) {
  psapply <- sapply
  if (!isTRUE(debug)) {
    clust <- inject(makeCluster(hmm$clspec,hmm$cltype,!!!hmm$clargs))
    psapply <- function(...) parSapply(clust,...)
    stopOnExit <- hmm$stopClusterOnError
    withr::defer({if (stopOnExit) stopCluster(clust)})
    if (!missing(seed)) {
        clusterSetRNGStream(clust,seed)
        mc.reset.stream()
    }
  }
  if (missing(mocc)) {
    if(length(hmm$deltaT)>1L)
      mocc <- length(hmm$deltaT)
  }

  hmm$npart <- 1L
  hmm$theta <- array(NA_real_,c(hmm$npart,nsim,mocc+1L))
  hmm$data <- array(NA_integer_,c(nsim,mocc))
  hmm$theta[,,1L] <- psapply(1L:hmm$nsubjects, \(subj) {
    hmm$popModels[[hmm$group(subj)]]$draw(hmm$npart)
  })


  for (it in 1L:mocc) {
    hmm$theta[,,it+1L] <- psapply(1L:hmm$nsubjects, \(subj) {
      hmm$growthModels[[hmm$action(subj,it)]]$draw(hmm$theta[,subj,it],
                                           hmm$delT(subj,it))
    })

    hmm$data[,it] <- psapply(1L:hmm$nsubjects, \(subj) {
      task <- hmm$task(subj,it)
      if (is.na(task)) return(NA_integer_)
      else {
        hmm$evidenceModels[[task]]$draw(hmm$theta[,subj,it+1L])
      }
    })

  }
  stopOnExit <- TRUE
  invisible(hmm)
}

