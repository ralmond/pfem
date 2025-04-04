PopulationModel <- R6Class(
    classname="PopulationModel",
    public = list(
        covergence=NA,
        lp=NA,
        name="<PopulationModel>",
        draw = function(npart,invariants=list()) {
          stop("Draw not implemented for ", class(self))
        },
        lprob = function(par=self$pvec,theta,weights,invariants=list()) {
          stop("Lprob not implemented for ", class(self))
        },
        mstep = function(theta,weights,invariants=list(),its=3,
                         control=list()) {
          control$maxits <- its
          result <- optim(self$pvec,\(pv)
                          self$lprob(pv,theta,weights,invariants),
                          control=control)
          if (result$convergence > 1)
            warning("Convergence issues with population model ",
                    self$name, "\n", result$message)
          self$lp <- result$value
          self$pvec <- result$par
          self$convergence <- result$convergence
        },
        print=function(...) {
          print(self$toString(...),...)
        },
        toString=function(...) {
          "<PopulationModel>"
        }
    ),
    active=list(
        pvec = function(value) {
          stop("Pvec active field not implemented for this class")
        }
    )
)


NormalPop <- R6Class(
    classname = "NormalPop",
    inherit=PopulationModel,
    public=list(
        initialize = function(name,mu,sigma) {
          self$name <- name
          self$mu <- mu
          self$sigma <- sigma
        },
        mu=0,
        sigma=1,
        draw = function(npart,invariants=list()) {
          mu <- self$mu
          sigma <- self$sigma
          rnorm(npart,mu,sigma)
        },
        lprob = function(par=self$pvec,theta,weights,invariants=list()) {
          mu <- par[1]
          sigma <- exp(par[2])
          sum(dnorm(theta,mu,sigma,log=TRUE)*weights)
        },
        toString=function(digits=2,...){
          paste0("<NormalPopulation: ",
          self$name, " ( ",round(self$mu,digits=digits),
                 ", ",round(self$sigma,digits=digits)," )>")
        }
    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$mu,log(self$sigma)))
          self$mu <- value[1]
          self$sigma <- exp(value[2])
        }
    )
)



