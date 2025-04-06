
EvidenceModel <- R6Class(
    classname="EvidenceModel",
    public = list(
        convergence=NA,
        lp=NA,
        name="<EvidenceModel>",
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
        },
        print=function(...) {
          print(self$toString(...),...)
        },
        toString=function(...) {
          "<EvidenceModel>"
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
  -t(apply(cbind(1,cuts,0),1,diff))
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
          log(probs[,Y+1L])
        },
        lprob = function(par=self$pvec,Y,theta,weights,context=list()) {
          probs <- cuts2probs(self$cuts(theta,a=exp(par[1]),b=par[-1]))
          sum(sapply(1L:length(theta),\(i) {
            weights[i]*log(probs[i,Y[i]+1L])
          }))
        },
        toString=function(digits=2,...){
          paste0("<GR: ", self$name, " ( ",
                 round(self$a,digits=digits),
                 ", ",paste(round(self$b,digits=digits),
                            collapse = ", "), " )>")
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

NormalScore <- R6Class(
    classname="NormalScore",
    inherit=EvidenceModel,
    public=list(
        initialize = function(name, bias, se) {
          self$name <- name
          self$bias <- bias
          self$se <- se
        },
        bias=0,
        se=01,
        draw = function(theta,context=list()) {
          rnorm(length(theta),theta+self$bias,self$se)
        },
        llike = function(Y,theta,context=list()) {
          dnorm(Y,theta+self$bias,self$se,log=TRUE)
        },
        lprob = function(par=self$pvec,Y,theta,weights,context=list()) {
          sum(dnorm(Y,theta+par[1],exp(par[2]),log=TRUE)*weights)
        },
        toString=function(digits=2,...){
          paste0("<NS: ", self$name, " ( ",
                 round(self$bias,digits=digits),
                 ", ",paste(round(self$se,digits=digits),
                            collapse = ", "), " )>")
        }

    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$bias,log(self$se)))
          self$bias <- value[1]
          self$se <- exp(value[2])
        }
    )
)

