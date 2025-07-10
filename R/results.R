longResults <- function (obj) UseMethod("longResults")

## Theta dimensions -- time,particle,subject
## Data dimensions -- time,subject
## Weight dimension -- particle,subject
## actions (it, subj)
## tasks(it,subj)
## times(it,subj)

longResults.HMM <- function (obj) {
  npp <- obj$npart*obj$nsubjects
  nrw <- npp*(obj$maxocc+1L)
  ptheta <- as.vector(obj$theta)

  Y <- rep(as.vector(rbind(NA,obj$data)),each=obj$npart)
  tasks <- obj$tasks
  if (!is.matrix(tasks) || nrow(tasks) < obj$maxocc)
    tasks <- matrix(as.vector(tasks),obj$maxocc,1L)
  tasks <- rbind(NA,tasks)
  if (ncol(tasks) > 1L)
    tasks <- sapply(1:ncol(tasks), \(subj)
                    rep(as.vector(tasks[,subj]),obj$npart))
  else
    tasks <- rep(as.vector(tasks),npp)
  if (ncol(obj$times) > 1L)
    alltimes <- sapply(1:ncol(obj$times), \(subj)
                       rep(as.vector(obj$times[,subj]),obj$npart))
  else
    alltimes <- rep(as.vector(obj$times),npp)

  subj<-rep(1:obj$nsubjects,each=obj$npart*(obj$maxocc+1L))
  occ<-rep(0L:obj$maxocc,npp)
  weights <- as.vector(obj$weights)
  if (length(obj$weights)==0L)
    weights <- rep(NA,npp)
  weights <-rep(weights,each=(obj$maxocc+1L))

  result <-data.frame(
    occ=occ,
    subj=as.factor(subj),
    particle=rep(rep(1L:obj$npart,each=(obj$maxocc+1L)),obj$nsubjects),
    time=alltimes,
    tasks=tasks,
    Y=Y,
    weights=weights,
    ptheta)
  names(result) <- c("occ","subj","particle","time","tasks","Y","weights",obj$thetaNames)
  result
}
## Theta dimensions -- time,nDims,particle,subject
## Data dimensions -- time,nDat,subject
## Weight dimension -- particle,subject
## Weight dimension -- particle,subject
## actions (it, subj)
## tasks(it,subj)
## time (it,subj)

longResults.HMM2 <- function (obj) {
  npp <- obj$npart*obj$nsubjects
  nrw <- npp*(obj$maxocc+1L)

  nDims <- obj$nDims
  nDat <- obj$nDat

  ptheta <- sapply(1:nDims, \(idim) as.vector(obj$theta[,idim,,]))
  colnames(ptheta) <- obj$thetaNames
  Y <- sapply(1:nDat,\(idat)
              rep(as.vector(rbind(NA,obj$data[,idat,])),each=obj$npart))
  colnames(Y) <- obj$dataNames
  tasks <- obj$tasks
  if (!is.matrix(tasks) || nrow(tasks) < obj$maxocc)
    tasks <- matrix(as.vector(tasks),obj$maxocc,1L)
  tasks <- rbind(NA,tasks)
  if (ncol(tasks) > 1L)
    tasks <- sapply(1:ncol(tasks), \(subj)
                    rep(as.vector(tasks[,subj]),obj$npart))
  else
    tasks <- rep(as.vector(tasks),npp)
  if (ncol(obj$times) > 1L)
    alltimes <- sapply(1:ncol(obj$times), \(subj)
                       rep(as.vector(obj$times[,subj]),obj$npart))
  else
    alltimes <- rep(as.vector(obj$times),npp)

  subj<-rep(1:obj$nsubjects,each=obj$npart*(obj$maxocc+1L))
  occ<-rep(0L:obj$maxocc,npp)
  weights <- as.vector(obj$weights)
  if (length(obj$weights)==0L)
    weights <- rep(NA,npp)
  weights <-rep(weights,each=(obj$maxocc+1L))

  result <-data.frame(
    occ=occ,
    subj=as.factor(subj),
    particle=rep(rep(1L:obj$npart,each=(obj$maxocc+1L)),obj$nsubjects),
    time=alltimes,
    tasks=tasks,
    Y,
    weights=weights,
    ptheta)
  result
}

longResults.IRTF <- function (obj) {
  Nquad <- length(obj$qpoints)
  Ntimes <- length(obj$tstar)
  Nsubj <- obj$nsubjects
  ptheta <- rep(obj$qpoints,Ntimes*Nsubj)

  Y <- rep(as.vector(obj$data),each=Nquad)
  tasks <- obj$tasks
  if (!is.matrix(tasks) || ncol(tasks) < obj$maxocc)
    tasks <- matrix(as.vector(tasks),1L,obj$maxocc)
  if (nrow(tasks) > 1L)
    tasks <- rep(as.vector(tasks),each=Nquad)
  else
    tasks <- rep(rep(as.vector(tasks),each=Nquad),Nsubj)
  alltimes <- rep(rep(obj$tstar,each=Nquad),each=Nsubj)

  subj<-rep(rep(1:Nsubj,each=Nquad),Ntimes)
  occ<-rep(0L:(Ntimes-1),each=Nquad*Nsubj)
  weights <- as.vector(obj$weights)
  if (length(obj$weights)==0L)
    weights <- rep(NA,Nquad*Ntimes*Nsubj)

  result <-data.frame(
    occ=occ,
    subj=as.factor(subj),
    qindex=rep(1L:Nquad,Nsubj*Ntimes),
    time=alltimes,
    tasks=NA,
    Y=NA,
    weights=weights,
    ptheta)
  names(result) <- c("occ","subj","particle","time","tasks","Y","weights",obj$thetaNames)
  result
}





avePart <- function (restab) {
  dplyr::group_by(restab,restab$occ,restab$subj) |>
    dplyr::summarize(time=min(.data$time),tasks=min(.data$tasks),
                     Y=min(.data$Y),
                     theta_bar=wtd.mean(.data$theta,.data$weights),
                     theta_sd=sqrt(wtd.var(.data$theta,.data$weights,
                                           normwt=TRUE)))
}

wtd.loess.part <- function (time,theta,weights,part)
  wtd.loess.noiter(time,theta,weights)$y[part==1]

part_loess <- function(tab) {
  nsubj <- max(tab$subj)
  mocc <- max(tab$occ)+1L
  smooths <- dplyr::group_by(tab,tab$subj) |>
    dplyr::group_map(~wtd.loess.part(.x$time,.x$theta,.x$weights,.x$particle))
  smooths <- matrix(unlist(smooths),nsubj,mocc,byrow=TRUE)
  smooths
}


col2matrix <- function(tab,col,fill=NA,minocc=1) {
  mocc <- max(tab$occ)
  nsubj <- max(tab$subj)
  nana <- rep(fill,mocc)
  dplyr::filter(tab,tab$occ>=minocc) |>
    dplyr::group_by(.data$subj) |>
    dplyr::group_map(~ c(.x[[col]],nana)[1:mocc]) |>
    unlist() |> matrix(nsubj,mocc,byrow=TRUE)
}

getDeltaT <- function(tab,fill=1) {
  mocc <- max(tab$occ)
  nsubj <- max(tab$subj)
  nana <- rep(fill,mocc)
  dplyr::group_by(tab,tab$subj) |>
    dplyr::group_map(~ c(fill,diff(.x$time),nana)[1:mocc]) |>
    unlist() |> matrix(nsubj,mocc,byrow=TRUE)
}
