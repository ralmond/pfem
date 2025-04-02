GrowthModel <- R6Class(
    classname="GrowthModel",
    public = list(
        covergence=NA,
        lprob=NA,
        draw = function(theta,deltaT,variants=list()) {
          stop("Draw not implemented for ", class(self))
        },
        lprob = function(par=self$pvec,theta,weights,deltaT,variants=list()) {
          stop("Lprob not implemented for ", class(self))
        },
        mstep = function(theta,weights,deltaT,variants=list(),its=3,
                         control=list()) {
          control$maxits <- its
          result <- optim(self$pvec,\(pv)
                          self$lprob(pv,theta,weights,deltaT,variants),
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
              
BrownianGrowth <- R6Class(
    classname="BrownianGrowth",
    inherit=GrowthModel,
    public=list(
        initialize = function(gain,inovSD) {
          self$gain <- gain
          self$inovSD <- inovSD
        },
        gain=0,
        inovSD=.1,
        draw = function(theta,deltaT,variants=list()){
          rnorm(length(theta),theta+self$gain*deltaT,
                self$inovSD*sqrt(deltaT))
        }
    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$gain,log(self$inovSD)))
          self$gain <- value[1]
          self$inovSD <- exp(value[2])
        }
    )
)



