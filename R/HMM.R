# !diagnostics suppress=self,private,super

name <- function(obj) {UseMethod("name")}
name.R6 <- function(obj){obj$name}
## Theta dimensions -- time,particle,subject
## Data dimensions -- time,subject
## Weight dimension -- particle,subject
## actions (it, subj)
## tasks(it,subj)
## time(it,subj)

HMM <- R6Class(
    "HMM",
    public=list(
      name="<HMM>",
      popModels=list(),
      groups = 1L,
      growthModels=list(),
      ## actions (it, subj)
      actions = matrix(1L,1L,1L),
      evidenceModels=list(),
      ## tasks(it,subj)
      tasks=matrix(1L,1L,1L),
      ## data(it,subj)
      data=matrix(NA_integer_,1L,1L),
      thetaNames="theta",
      npart=1L,
      ## theta(it,part,subj)
      theta = numeric(),
      ## weights (part,subj)
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
      delT = function(it,subj) {
        if (ncol(private$DeltaT)==1L) subj <- 1L
        if (nrow(private$DeltaT)==1L) it <- 1L
        private$DeltaT[it,subj]
      },
      group = function (subj) {
        if (length(self$groups)==1L) return(self$groups)
        self$groups[subj]
      },
      action = function(it,subj) {
        if (!is.matrix(self$actions)) {
          if (length(self$actions)==1L)
            self$action <- matrix(self$actions,1L,1L)
          else if (length(self$actions)==self$maxocc)
            self$action <- matrix(self$action,self$maxocc,1L)
          else if (length(self$actions)==self$nsubjects)
            self$action <- matrix(self$actions,1L,self$nsubjects)
          else stop ("Unexpected size for action matrix.  Should be ",
                     self$maxocc, "x", self$nsubjects,".")
        }
        if (ncol(self$actions)==1L) subj <- 1L
        if (nrow(self$actions)==1L) it <- 1L
        self$actions[it,subj]
      },
      task = function(it,subj) {
        if (!is.matrix(self$tasks)) {
          if (length(self$tasks)==1L)
            self$task <- matrix(self$tasks,1L,1L)
          else if (length(self$tasks)==self$maxocc)
            self$tasks <- matrix(self$task,self$maxocc,1L)
          else if (length(self$tasks)==self$nsubjects)
            self$tasks <- matrix(self$tasks,1L,self$nsubjects)
          else stop ("Unexpected size for task matrix.  Should be ",
                     self$maxocc, "x", self$nsubjects,".")
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
      },
      drawPop = function (subj,npart) {
        self$popModels[[self$group(subj)]]$draw(npart)
      },
      drawGrowth = function (it,subj,theta0=self$theta[it,,subj]) {
        if (self$delT(it,subj)==0) theta0
        else
          self$growthModels[[self$action(it,subj)]]$draw(theta0,
                                                       self$delT(it,subj))
      },
      evalEvidence= function (it,subj,theta1=self$theta[it+1L,,subj],
                              Y=self$data[it,subj]) {
        task <- self$task(it,subj)
        if (is.na(Y) || is.na(task)) return(rep(0,length(theta1)))
        else {
          self$evidenceModels[[task]]$llike(Y,theta1)
        }
      },
      drawEvidence= function (it,subj,theta1=self$theta[it+1L,,subj]) {
        task <- self$task(it,subj)
        if (is.na(task)) return(rep(NA_integer_,length(theta1)))
        else {
          self$evidenceModels[[task]]$draw(theta1)
        }
      }
    ),
    private=list(
        ##  deltaT (it,subj)
        DeltaT = numeric(),
        ## time (it, subj)
        Times = numeric(),
        Nsubjects = 1,
        Maxtimes = 1
    ),
    active=list(
        deltaT = function(value) {
          if (missing(value)) return(private$DeltaT)
          if (!is.matrix(value)) {
            value <- matrix(value,length(value),1L)
          }
          private$DeltaT <- value
          private$Times <- t(apply(rbind(0,value),2,cumsum))
        },
        times = function(value) {
          if (missing(value)) return(private$Times)
          if (!is.matrix(value))
            value <- matrix(value,length(value),1L)
          private$Times <- value
          private$DeltaT <- t(apply(value,2,diff))
        },
        nsubjects = function() {ncol(self$data)},
        maxocc = function(value) {nrow(self$data)},
        ## weights(particle,subject)
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


particleFilter <- function (hmm, npart=hmm$npart, seed=NULL,
                            debug=FALSE,weightLog=FALSE) {
  UseMethod("particleFilter")
}


particleFilter.HMM <- function (hmm, npart=hmm$npart, seed=NULL,
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
  if(isTRUE(weightLog)) hmm$weightLog <- list()

  hmm$npart <- npart
  hmm$lweights <- matrix(0,npart,hmm$nsubjects)
  hmm$theta <- array(NA_real_,c(hmm$maxocc+1L,npart,hmm$nsubjects))

  filtres <- plapply(1L:hmm$nsubjects, \(subj) {
    stheta <- array(NA_real_,c(hmm$maxocc+1L,npart))
    sweights <- rep(0,npart)
    stheta[1L,] <- hmm$drawPop(subj,npart)
    swl <- list()

    for (it in 1L:hmm$maxocc) {
      stheta[it+1L,] <- hmm$drawGrowth(subj,it,stheta[it,])

      sweights <- hmm$lweights +
        hmm$evalEvidence(subj,it,stheta[it+1L,],hmm$data[it,subj])
      if(isTRUE(weightLog))
        swl <- c(swl,sweights)
    }
    list(theta=stheta,weights=sweights,wl=swl)
  })
  for (subj in 1L:hmm$nsubjects) {
    hmm$theta[,,subj] <- filtres[[subj]]$theta
    hmm$lweights[,subj] <- filtres[[subj]]$weights
    hmm$weightLog <- c(hmm$weightLog,filtres[[subj]]$wl)
  }
  stopOnExit <- TRUE
  invisible(hmm)
}


# Can reuse stats::simulate
#simulate <- function(object,nsim=hmm$nsubjects,seed=NULL,...,mocc=hmm$maxocc)
#  UseMethod("simulate")

simulate.HMM <- function(object,nsim=object$nsubjects,seed=NULL,mocc=2L,...,
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
  object$theta <- array(NA_real_,c(mocc+1L,object$npart,nsim))
  object$data <- array(NA_integer_,c(mocc,nsim))
  filtres <- plapply(1L:nsim, \(subj) {
    stheta <- rep(NA_real_,mocc+1L)
    sdata <- rep(NA_integer_,mocc)

    stheta[1L] <- object$drawPop(subj,1L)
    for (it in 1L:mocc) {
      stheta[it+1L] <- object$drawGrowth(it,subj,stheta[it])
      sdata[it] <- object$drawEvidence(it,subj,stheta[it+1L])
    }
    list(theta=stheta,data=sdata)
  })
  for (subj in 1L:nsim) {
    object$theta[,subj] <- filtres[[subj]]$theta
    object$data[,subj] <- filtres[[subj]]$data
  }

  stopOnExit <- TRUE
  invisible(object)
}

