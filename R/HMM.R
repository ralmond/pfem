HMM <- R6Class(
    "HMM",
    public=list(
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
        weights = numeric(),
        initialize=function(popModels,growthModels,evidenceModels) {
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
          if (!is.matrix(self$task)) {
            if (length(self$task)==1L)
              self$task <- matrix(self$task,1L,1L)
            else if (length(self$task)==self$maxocc)
              self$task <- matrix(self$task,1L,self$maxocc)
            else if (length(self$task)==self$nsubjects)
              self$task <- matrix(self$task,self$nsubjects,1L)
            else stop ("Unexpected size for task matrix.  Should be ",
                       self$nsubjects, "x", self$maxocc,".")
          }
          if (nrow(self$task)==1L) subj <- 1L
          if (ncol(self$task)==1L) it <- 1L
          self$task[subj,it]
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
        maxocc = function(value) {ncol(self$data)}
    )
)


particleFilter <- function (hmm, npart=hmm$npart) {
  UseMethod("particleFilter")
}


#' Title
#'
#' @param hmm
#' @param npart
#'
#' @return
#' @export
#'
#' @examples
particleFilter.HMM <- function (hmm, npart=hmm$npart) {
  hmm$npart <- npart
  hmm$weights <- matrix(0,npart,hmm$nsubjects)
  hmm$theta <- array(NA_real_,c(npart,hmm$subjects,hmm$maxocc))
  hmm$theta[,,1L] <- parallel::parSapply(1L:hmm$nsubjects, \(subj) {
    hmm$popModels[[hmm$group(subj)]]$draw(npart)
  })


  for (it in 1L:hmm$maxocc) {
    hmm$theta[,,it+1L] <- parallel::parSapply(1L:hmm$nsubject, \(subj) {
      hmm$growthModels[[hmm$action(subj,it)]]$draw(hmm$theta[,subj,it],
                                           hmm$delT(subj,it))
    })

    hmm$weights <- hmm$weights +
      parallel::parSapply(1L:hmm$nsubjects, \(subj) {
        task <- hmm$task(subj,it)
        Y <- hmm$data[subj,it]
        if (is.na(Y) || is.na(task)) return(0)
        else {
          hmm$evidenceModels[[task]]$llike(Y,hmm$theta[,subj,it+1L])
        }
      })

  }
  invisible(hmm)
}

longResults <- function (hmm) UseMethod("longResults")

longResults.HMM <- function (hmm) {
  npp <- hmm$npart*hmm$nsubjects
  nrw <- npp*(hmm$maxocc+1L)
  ptheta <- as.vector(hmm$theta)

  Y <- rep(as.vector(cbind(NA,hmm$data)),each=hmm$npart)

  if (nrow(hmm$times) > 1L)
    alltimes <- rep(as.vector(hmm$times),each=hmm$npart)
  else
    alltimes <- rep(as.vector(hmm$times),each=npp)

  result <-data.frame(
      subj=rep(rep(1:hmm$nsubjects,each=hmm$npart),hmm$maxocc+1L),
      occ=rep(0L:hmm$maxcocc,each=npp),
      time=alltimes,
      Y=Y,
      weights=rep(as.vector(hmm$weights),hmm$maxocc+1L),
      ptheta)
  names(result) <- c("subj","occ","time","weights",hmm$thetaNames)
  result
}

# Can reuse stats::simulate
#simulate <- function(object,nsim=hmm$nsubjects,seed=NULL,...,mocc=hmm$maxocc)
#  UseMethod("simulate")

simulate.HMM <- function(hmm,nsim=hmm$nsubjects,seed=NULL,mocc=hmm$maxocc,...) {
  if (!missing(seed)) {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  hmm$npart <- 1L
  hmm$theta <- array(NA_real_,c(hmm$npart,nsim,hmm$maxocc))
  hmm$data <- array(NA_integer_,c(nsim,hmm$maxocc))
  hmm$theta[,,1L] <- parallel::parSapply(1L:hmm$nsubjects, \(subj) {
    hmm$popModels[[hmm$group(subj)]]$draw(hmm$npart)
  })


  for (it in 1L:hmm$maxocc) {
    hmm$theta[,,it+1L] <- parallel::parSapply(1L:hmm$nsubject, \(subj) {
      hmm$growthModels[[hmm$action(subj,it)]]$draw(hmm$theta[,subj,it],
                                           hmm$delT(subj,it))
    })

    hmm$data[,it] <- parallel::parSapply(1L:hmm$nsubjects, \(subj) {
      task <- hmm$task(subj,it)
      if (is.na(task)) return(NA_integer_)
      else {
        hmm$evidenceModels[[task]]$draw(hmm$theta[,subj,it+1L])
      }
    })

  }
  invisible(hmm)

}

