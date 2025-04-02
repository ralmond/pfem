PopulationModel <- R6Class(
    classname="PopulationModel",
    public = list(
        covergence=NA,
        lprob=NA,
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
                    self$Name, "\n", result$message)
          self$convergence <- result$convergence
          self$pvec <- result$par
        }
    ),
    private= list(
        namE=character(),
        imod=0
    ),
    active=list(
        name = function (value) {
          if (missing(value)) {
            if (length(private$namE) ==0L)
              private$namE <- paste0(class(self)[1],populationModel$iMod)
            private$namE
          } else {
            private$namE <- value
          }
        },
        iMod = function() {
          private$iMod <- 1+ private$iMod
          private$iMod
        },
        pvec = function(value) {
          stop("Pvec active field not implemented for this class")
        }
    )
)
              

NormalPop <- R6Class(
    classname = "NormalPop",
    inherit=PopulationModel,
    public=list(
        initialize = function(mu,sigma) {
          self$mu <- mu
          self$sigma <- sigma
        },
        mu=0,
        sigma=1,
        draw = function(npart,invariants=list())
          rnorm(npart,self$mu,self$simga)
    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$mu,log(self$sigma)))
          self$mu <- value[1]
          self$sigma <- exp(value[2])
        }
    )
)


                   
