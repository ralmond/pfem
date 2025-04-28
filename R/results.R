longResults <- function (obj) UseMethod("longResults")

longResults.HMM <- function (obj) {
  npp <- obj$npart*obj$nsubjects
  nrw <- npp*(obj$maxocc+1L)
  ptheta <- as.vector(obj$theta)

  Y <- rep(as.vector(cbind(NA,obj$data)),each=obj$npart)
  tasks <- obj$tasks
  if (!is.matrix(tasks) || ncol(tasks) < obj$maxocc)
    tasks <- matrix(as.vector(tasks),1L,obj$maxocc)
  tasks <- cbind(NA,tasks)
  if (nrow(tasks) > 1L)
    tasks <- rep(as.vector(tasks),each=obj$npart)
  else
    tasks <- rep(as.vector(tasks),each=npp)
  if (nrow(obj$times) > 1L)
    alltimes <- rep(as.vector(obj$times),each=obj$npart)
  else
    alltimes <- rep(as.vector(obj$times),each=npp)

  subj<-rep(rep(1:obj$nsubjects,each=obj$npart),obj$maxocc+1L)
  occ<-rep(0L:obj$maxocc,each=npp)
  weights <- as.vector(obj$weights)
  if (length(obj$weights)==0L)
    weights <- rep(NA,npp)
  weights <-rep(weights,obj$maxocc+1L)

  result <-data.frame(
    occ=occ,
    subj=as.factor(subj),
    particle=rep(1L:obj$npart,obj$nsubjects*(obj$maxocc+1L)),
    time=alltimes,
    tasks=tasks,
    Y=Y,
    weights=weights,
    ptheta)
  names(result) <- c("occ","subj","particle","time","tasks","Y","weights",obj$thetaNames)
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
  dplyr::group_by(restab,dplyr::matches("occ"),dplyr::matches("subj")) |>
    dplyr::summarize(time=min(time),tasks=min(tasks),
                     Y=min(Y),
                     theta_bar=wtd.mean(theta,weights),
                     theta_sd=sqrt(wtd.var(theta,weights,normwt=TRUE)))
}

wtd.loess.part <- function (time,theta,weights,part)
  wtd.loess.noiter(time,theta,weights)$y[part==1]

part_loess <- function(tab) {
  nsubj <- max(tab$subj)
  mocc <- max(tab$occ)+1L
  smooths <- dplyr::group_by(tab,dplyr::matches("subj")) |>
    dplyr::group_map(~wtd.loess.part(.x$time,.x$theta,.x$weights,.x$particle))
  smooths <- matrix(unlist(smooths),nsubj,mocc,byrow=TRUE)
  smooths
}


col2matrix <- function(tab,col,fill=NA,minocc=1) {
  mocc <- max(tab$occ)
  nsubj <- max(tab$subj)
  nana <- rep(fill,mocc)
  dplyr::filter(tab,occ>=minocc) |>
    dplyr::group_by(subj) |>
    dplyr::group_map(~ c(.x[[col]],nana)[1:mocc]) |>
    unlist() |> matrix(nsubj,mocc,byrow=TRUE)
}

getDeltaT <- function(tab,fill=1) {
  mocc <- max(tab$occ)
  nsubj <- max(tab$subj)
  nana <- rep(fill,mocc)
  dplyr::group_by(tab,subj) |>
    dplyr::group_map(~ c(fill,diff(.x$time),nana)[1:mocc]) |>
    unlist() |> matrix(nsubj,mocc,byrow=TRUE)
}
