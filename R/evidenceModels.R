
EvidenceModel <- R6Class(
    classname="EvidenceModel",
    public = list(
        convergence=NA,
        lp=NA,
        draw = function(theta,context=list()) {
          stop("Draw not implemented for ", class(self))
        },
        llike = function(Y,theta,context=list()) {
          stop("Lprob not implemented for ", class(self))
        },
        lprob = function(par=self$pvec,Y,theta,weights,context=list()) {
          stop("Lprob not implemented for ", class(self))
        },
        mstep = function(Y,theta,weights,context=list(),its=3,
                         control=list()) {
          control$maxits <- its
          result <- optim(self$pvec,\(pv)
                          self$lprob(pv,Y,theta,weights,context),
                          control=control)
          if (result$convergence > 1)
            warning("Convergence issues with population model ",
                    self$Name, "\n", result$message)
          self$pvec <- result$par
          self$lp <- result$value
          self$convergence <- result$convergence
        }
    ),
    active=list(
        pvec = function(value) {
          stop("Pvec active field not implemented for this class")
        }
    )
)


logit <- function (p) log(p/(1-p))/1.7
invlogit <- function (x) 1/(1+exp(-1.7*x))

cuts2probs <- function (cuts) {
  if (!is.matrix(cuts)) cuts <- matrix(cuts,1L,length(cuts))
  t(apply(cbind(1,cuts,0),1,diff))
}


GradedResponse <- R6Class(
    classname="GradedResponse",
    inherit=EvidenceModel,
    public=list(
        initialize = function(name, a, b) {
          self$name <- name
          self$a <- a
          self$b <- b
        },
        a=1,
        b=0,
        cuts = function(theta,a=self$a,b=self$b) {
          invlogit(outer(a*theta,b,"-"))
        },
        draw = function(theta,context=list()) {
          rowSums(sweep(self$cuts(theta),1,runif(length(theta)),"<"))
        },
        llike = function(Y,theta,context=list()) {
          probs <- cuts2probs(self$cuts(theta))
          probs[,Y]
        },
        lprob = function(par=self$pvec,Y,theta,weights,context=list()) {
          probs <- cuts2probs(self$cuts(theta,a=exp(par[1]),b=par[-1]))
          sum(sapply(1L:length(theta),\(i) {
            weights[i]*log(probs[i,Y[i]])
          }))
        }
    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(log(self$a),self$b))
          self$a <- exp(value[1])
          self$b <- value[-1]
        }
    )
)

