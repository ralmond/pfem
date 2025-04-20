IRTF <- R6Class(
    "IRTF",
    public=list(
        name="IRT Filter",
        thetaNames="theta",
        qpoints = qnorm(((0:15)+.5)/16),
        qbounds = qnorm((0:16)/16),
        tstar = numeric(),
        wfun = NULL,
        popModels=list(),
        groups=1L,
        evidenceModels = list(),
        tasks = matrix(1L,1L,1L),
        data=matrix(NA_integer_,1L,1L),
        times=matrix(1,1L,1L),
        lweights=numeric(),
        cltype=ifelse(.Platform$OS.type=="windows","PSOCK","FORK"),
        clspec=getOption("mc.cores",2L),
        clargs=list(),
        stopClusterOnError=TRUE,
        initialize=function(name,Nquad=16,
                            popModels=list(),
                            evidenceModels=list()) {
          self$name <- name
          self$qpoints <- qnorm(((1:Nquad)-.5)/Nquad)
          self$qbounds <- qnorm((0:Nquad)/Nquad)
          self$popModels <- popModels
          if (!is.list(popModels))
            self$popModels <- list(popModels)
          self$evidenceModels <- evidenceModels
          if (!is.list(evidenceModels))
            self$evidenceModels <- list(evidenceModels)
        },
        time = function(subj,it) {
          if (nrow(self$times)==1L) subj <- 1L
          self$times[subj,it+1L]
        },
        group = function (subj) {
          if (length(self$groups)==1L) return(self$groups)
          self$groups[subj]
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
          paste0("<IRTF: ",self$name,": ",
                 self$nsubjects, " x ",
                 self$macocc, " >")
        },
        print=function(...) {
          print(self$toString(...),...)
        },
        evalEvidence= function (subj,it) {
          task <- self$task(subj,it)
          Y <- self$data[subj,it]
          if (is.na(Y) || is.na(task)) return(0)
          else {
            self$evidenceModels[[task]]$llike(Y,self$qpoints)
          }
        }
    ),
    active=list(
        nsubjects = function() {nrow(self$data)},
        maxocc = function(value) {ncol(self$data)},
        weights = function() {
          weights <- exp(self$lweights)
          if (length(weights)==0L)
            return(weights)
          if(!is.array(weights))
            return(weights/sum(weights))
          aperm(sweep(weights,c(2,3),apply(weights,c(2,3),sum),"/"),
                c(1,3,2))
        }
    )
)

irtf <- function (irf,tstar=irf$tstar,wfun=irf$wfun,
                  debug=FALSE) {
  ## Setup Clusters
  plapply <- lapply
  if (!debug) {
    clust <- inject(makeCluster(irf$clspec,irf$cltype,!!!irf$clargs))
    stopOnExit <- irf$stopClusterOnError
    withr::defer({if (stopOnExit) stopCluster(clust)})
    plapply <- function(...) parLapply(clust,...)
  }

  irf$tstar <- tstar
  irf$wfun <- wfun
  Nsubj <- irf$nsubjects
  Nt <- length(tstar)
  Qcount <- length(irf$qpoints)
  theta <- irf$qpoints
  irf$lweights <- array(0,c(Qcount,Nt,Nsubj))

  wlist <- plapply(1:Nsubj, \(subj) {
    lwts <- matrix(0,Qcount,Nt)

    if (length(irf$popModels) > 0L) {
      wts <- diff(irf$popModels[[irf$group(subj)]]$cdf(irf$qbounds))
      lwts <- outer(wts, wfun(irf$time(subj,0L),tstar),"*")
    }

    for (it in 1L:irf$maxocc) {
      task <- irf$task(subj,it)
      tim <- irf$time(subj,it)
      if (is.na(tim)) break
      Y <- irf$data[subj,it]
      if (!is.na(Y) && !is.na(task))
        lwts <- lwts + outer(irf$evalEvidence(subj,it),
                             wfun(tim,tstar),"*")

    }
    lwts
  })
  for (subj in 1:Nsubj)
    irf$lweights[,,subj] <- wlist[[subj]]

    stopOnExit <- TRUE
  invisible(irf)
}




