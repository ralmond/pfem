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
    subj=as.factor(subj),
    occ=occ,
    particle=rep(1L:hmm$npart,hmm$nsubjects*(hmm$maxocc+1L)),
    time=alltimes,
    tasks=tasks,
    Y=Y,
    weights=weights,
    ptheta)
  names(result) <- c("subj","occ","particle","time","tasks","Y","weights",hmm$thetaNames)
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
