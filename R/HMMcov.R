# HMM <- R6Class(
#     "PFEM",
#     public=list(
#         popModels=list(),
#         groups = 1L,
#         growthModels=list(),
#         actions = matrix(1L,1L,1L),
#         evidenceModels=list(),
#         tasks=matrix(1L,1L,1L),
#         data=matrix(NA_integer_,1L,1L),
#         Ktheta=1L,
#         thetaNames="theta",
#         npart=1L,
#         theta = numeric(),
#         weights = numeric(),
#         initialize=function(popModels,growthModels,evidenceModels) {
#           self$popModels <- popModels
#           if (!is.list(popModels)) self$popModels <- list(popModels)
#           self$growthModels <- growthModels
#           if (!is.list(growthModels)) self$growthModels <- list(growthModels)
#           self$evidenceModels <- evidenceModels
#           if (!is.list(evidenceModels)) self$evidenceModels <-
#                                           list(evidenceModels)
#         },
#         Zinv=function(subj) {
#           if (nrow(private$Zinvariants)==0L) return(NULL)
#           private$Zinvariants[subj,]
#         },
#         Zvar=function(subj,it) {
#           if (nrow(private$Zvariants)==0L) return(NULL)
#           private$Zvariants[private$Zvariants$subj==subj & private$Zvariants$occ==it,]
#         },
#         ContX=function(subj,it) {
#           if (nrow(private$Contexts)==0L) return(NULL)
#           private$Contexts[private$Contexts$subj==subj & private$Contexts$occ==it,]
#         },
#         delT = function(subj,it) {
#           if (nrow(private$DeltaT)==1L) subj <- 1L
#           if (ncol(private$DeltaT)==1L) it <- 1L
#           private$DeltaT[subj,it]
#         },
#         group = function (subj) {
#           if (length(self$groups)==1L) return(self$groups)
#           self$groups[subj]
#         },
#         action = function(subj,it) {
#           if (!is.matrix(self$actions)) {
#             if (length(self$actions)==1L)
#               self$action <- matrix(self$actions,1L,1L)
#             else if (length(self$actions)==self$maxocc)
#               self$action <- matrix(self$action,1L,self$maxocc)
#             else if (length(self$actions)==self$nsubjects)
#               self$action <- matrix(self$actions,self$nsubjects,1L)
#             else stop ("Unexpected size for action matrix.  Should be ",
#                        self$nsubjects, "x", self$maxocc,".")
#           }
#           if (nrow(self$actions)==1L) subj <- 1L
#           if (ncol(self$actions)==1L) it <- 1L
#           self$actions[subj,it]
#         },
#         task = function(subj,it) {
#           if (!is.matrix(self$task)) {
#             if (length(self$task)==1L)
#               self$task <- matrix(self$task,1L,1L)
#             else if (length(self$task)==self$maxocc)
#               self$task <- matrix(self$task,1L,self$maxocc)
#             else if (length(self$task)==self$nsubjects)
#               self$task <- matrix(self$task,self$nsubjects,1L)
#             else stop ("Unexpected size for task matrix.  Should be ",
#                        self$nsubjects, "x", self$maxocc,".")
#           }
#           if (nrow(self$task)==1L) subj <- 1L
#           if (ncol(self$task)==1L) it <- 1L
#           self$task[subj,it]
#         }
#     ),
#     private=list(
#         Zinvariants=data.frame(),
#         Zvariants=data.frame(),
#         Contexts=data.frame(),
#         DeltaT = numeric(),
#         Times = numeric(),
#         Nsubjects = 1,
#         Maxtimes = 1
#     ),
#     active=list(
#         deltaT = function(value) {
#           if (missing(value)) return(private$DeltaT)
#           if (!is.matrix(value)) {
#             value <- matrix(value,1L,length(value))
#           }
#           private$DeltaT <- value
#           private$Times <- t(apply(cbind(0,value),1,cumsum))
#         },
#         times = function(value) {
#           if (missing(value)) return(private$Times)
#           if (!is.matrix(value))
#             value <- matrix(value,1L,length(value))
#           private$Times <- value
#           private$DeltaT <- t(apply(value,1,diff))
#         },
#         nsubjects = function() {nrow(self$data)},
#         maxocc = function(value) {ncol(self$data)}
#  #        invariants = function(value) {
#  # #         if (missing(value)) {return (private$Zinvariants)}
#  # #         if (!is.data.frame(value)) {
#  # #           stop("Invariant contants must be a data frame.")
#  # #         }
#  #          if (nrow(value) > 0L && nrow(value) != self$nsubjects) {
#  #            stop("Number of rows for time invariant constants must match number of subjects.")
#  #          }
#  #          private$Zinvarants <- value
#  #        },
#  #        variants= function(value) {
#  #          if (missing(value)) return (private$Zvariants)
#  #          if (!is.data.frame(value))
#  #            stop("Time varying contants must be a data frame.")
#  #          if (nrow(value) > 0L) {
#  #            if (is.null(value$subj) || is.null(val$occ)) {
#  #              stop("Time invariant constants must have `subj` and `occ` columns.")
#  #            }
#  #          }
#  #          private$Zvarants <- value
#  #        },
#  #        contexts= function(value) {
#  #          if (missing(value)) return (private$Contexts)
#  #          if (!is.data.frame(value))
#  #            stop("Time varying contants must be a data frame.")
#  #          if (nrow(value) > 0L) {
#  #            if (is.null(value$subj) || is.null(val$occ)) {
#  #              stop("Time invariant constants must have `subj` and `occ` columns.")
#  #            }
#  #          }
#  #          private$Contexts <- value
#  #        }
#     )
# )
#
#
#
# particleFilter <- function (hmm, npart=hmm$npart) {
#   hmm$npart <- npart
#   hmm$weights <- matrix(0,npart,hmm$nsubjects)
#   theta <- array(NA_real_,c(npart,hmm$Ktheta,hmm$subjects,hmm$maxocc))
#   theta[,,,1L] <- parSapply(1L:hmm$nsubjects, \(subj) {
#     popModels[[hmm$group(subj)]]$draw(npart,hmm$Zinv(subj))
#   })
#
#
#   for (it in 1L:hmm$maxocc) {
#     theta[,,it+1L] <- parSapply(1L:hmm$nsubject, \(subj) {
#       growthModels[[hmm$action(isub,it)]]$draw(theta[,,subj,it],
#                                            hmm$delT(subj,it),
#                                            hmm$Zvar(subj,it))
#     })
#
#     weights <- weights +
#       parSapply(1L:hmm$nsubjects, \(subj) {
#         task <- hmm$task(subj,it)
#         Y <- hmm$data[isub,it]
#         if (is.na(Y) || is.na(task)) return(0)
#         else {
#           evidenceModels[[task]]$llike(Y,theta[,,subj,it+1L],
#                                        hmm$ContX(subj,it))
#         }
#       })
#
#   }
#   invisible(hmm)
# }
#
# longResults <- function (hmm) {
#   npp <- hmm$npart*hmm$nsubjects
#   nrw <- npp*(hmm$maxocc+1L)
#   ptheta <- aperm(hmm$theta,c(1,2,4,3))
#   dim(ptheta) <- c(nrw,hmm$Ktheta)
#
#   if (nrow(hmm$times) > 1L)
#     alltimes <- rep(as.vector(hmm$times),each=hmm$npart)
#   else
#     alltimes <- rep(as.vector(hmm$times),each=npp)
#
#   result <-data.frame(
#       subj=rep(rep(1:hmm$nsubjects,each=hmm$npart),hmm$maxocc+1L),
#       occ=rep(0L:hmm$maxcocc,each=npp),
#       time=alltimes,
#       weights=rep(as.vector(hmm$weights),hmm$maxocc+1L),
#       ptheta)
#   names(result) <- c("subj","occ","time","weights",hmm$thetaNames)
#   result
# }
