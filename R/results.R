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
    occ=occ,
    subj=as.factor(subj),
    particle=rep(1L:hmm$npart,hmm$nsubjects*(hmm$maxocc+1L)),
    time=alltimes,
    tasks=tasks,
    Y=Y,
    weights=weights,
    ptheta)
  names(result) <- c("occ","subj","particle","time","tasks","Y","weights",hmm$thetaNames)
  result
}

longResults.IRTF <- function (hmm) {
  Nquad <- length(hmm$qpoints)
  Ntimes <- length(hmm$tstar)
  Nsubj <- hmm$nsubjects
  ptheta <- rep(hmm$qpoints,Ntimes*Nsubj)

  Y <- rep(as.vector(hmm$data),each=Nquad)
  tasks <- hmm$tasks
  if (!is.matrix(tasks) || ncol(tasks) < hmm$maxocc)
    tasks <- matrix(as.vector(tasks),1L,hmm$maxocc)
  if (nrow(tasks) > 1L)
    tasks <- rep(as.vector(tasks),each=Nquad)
  else
    tasks <- rep(rep(as.vector(tasks),each=Nquad),Nsubj)
  alltimes <- rep(rep(hmm$tstar,each=Nquad),each=Nsubj)

  subj<-rep(rep(1:Nsubj,each=Nquad),Ntimes)
  occ<-rep(0L:(Ntimes-1),each=Nquad*Nsubj)
  weights <- as.vector(hmm$weights)
  if (length(hmm$weights)==0L)
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
  names(result) <- c("occ","subj","particle","time","tasks","Y","weights",hmm$thetaNames)
  result
}





avePart <- function (restab) {
  dplyr::group_by(restab,occ,subj) |>
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
  smooths <- dplr::group_by(tab,subj) |>
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
