\name{wtd.stats}
\alias{wtd.Ecdf}
\alias{wtd.loess.noiter}
\alias{wtd.mean}
\alias{wtd.quantile}
\alias{wtd.rank}
\alias{wtd.table}
\alias{wtd.var}
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{

  These are inherited from \code{Hmisc::\link[Hmisc]{wtd.stats}}.  Note that for
  particle filter or irtf output, \code{normwt=TRUE} produces the correct result.

}
\usage{
wtd.Ecdf(x, weights = NULL, type = c("i/n", "(i-1)/(n-1)", "i/(n+1)"),
    normwt = FALSE, na.rm = TRUE)
wtd.loess.noiter(x, y, weights = rep(1, n), span = 2/3, degree = 1,
    cell = 0.13333, type = c("all", "ordered all", "evaluate"),
    evaluation = 100, na.rm = TRUE)
wtd.mean(x, weights = NULL, normwt = "ignored", na.rm = TRUE)
wtd.quantile(x, weights = NULL, probs = c(0, 0.25, 0.5, 0.75, 1),
       type = c("quantile", "(i-1)/(n-1)", "i/(n+1)", "i/n"),
       normwt = FALSE, na.rm = TRUE)
wtd.rank(x, weights = NULL, normwt = FALSE, na.rm = TRUE)
wtd.table(x, weights = NULL, type = c("list", "table"), normwt = FALSE,
    na.rm = TRUE)
wtd.var(x, weights = NULL, normwt = FALSE, na.rm = TRUE, method = c("unbiased",
    "ML"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,y}{A data vector}
  \item{weights}{A vector of weights}
  \item{probs}{A vector of probabilities for computing quantiles.}
  \item{type}{See \code{Hmisc::\link[Hmisc]{wtd.stats}} for interpretation.}
  \item{normwt}{Use \code{normwt=TRUE} with the filter weights.}
  \item{na.rm}{Logical value, should \code{NA}s be removed.}
  \item{span}{Parameter of loess.}
  \item{degree}{Parameter of loess.}
  \item{evaluation}{Parameter of loess.}
  \item{method}{Algorithm for calculating variances.}

}
\value{
  Depends on the function, see \code{Hmisc::\link[Hmisc]{wtd.stats}}
}
\author{Imported from \code{Hmisc} package, authors:  Frank Harrell and
Benjamin Tyner.

This docfile written by Russell Almond.
}
\seealso{
  See \code{Hmisc::\link[Hmisc]{wtd.stats}}
  \code{\link{avePart}}
}
\examples{
make_pops <- function() {
  list(NormalPop$new("P1",0,1),
       NormalPop$new("P2",1,1))
}

make_growths <- function() {
  list(BrownianGrowth$new("G1",1,1),
       BrownianGrowth$new("G2",2,1))
}

make_ems <- function() {
  list(NormalScore$new("I1",0,.5),
       NormalScore$new("I2",.25,.25))
}

hmm<-HMM$new("stimulator",
             list(NormalPop$new("P1",-1/2,1),
                  NormalPop$new("P2",1/2,1)),
             list(BrownianGrowth$new("G1",1,1/100)),
             list(NormalScore$new("N0",0,1)))
hmm$clspec <- 1
hmm$groups <- 1:2
hmm$deltaT<-rep(1,3)
hmm$data <- outer(c(-1,1)/2,1:3,"+")
particleFilter(hmm,5)
longResults(hmm) |>
  dplyr::group_by(occ,subj) |>
  dplyr::summarize(theta_bar=wtd.mean(theta,weights),
                   theta_sd=sqrt(wtd.var(theta,weights,normwt=TRUE)),
                   theta_lb=wts.quantile(theta,weights,probs=.025),
                   theta_ub=wts.quantile(theta,weights,probs=.925))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
