GrowthModel <- R6Class(
    classname="GrowthModel",
    public = list(
        covergence=NA,
        lp=NA,
        name="<GrowthModel>",
        draw = function(theta,deltaT,variants=list()) {
          stop("Draw not implemented for ", class(self))
        },
        lprob = function(par=self$pvec,theta0,theta1,weights,deltaT,
                         variants=list()) {
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
                    self$name, "\n", result$message)
          self$pvec <- result$par
          self$lp <- result$value
          self$convergence <- result$convergence
        },
        print=function(...) {
          print(self$toString(...),...)
        },
        toString=function(...) {
          "<GrowthModel>"
        }

    ),
    active=list(
        pvec = function(value) {
          stop("Pvec active field not implemented for this class")
        }
    )
)

BrownianGrowth <- R6Class(
    classname="BrownianGrowth",
    inherit=GrowthModel,
    public=list(
        initialize = function(name,gain,inovSD) {
          self$name <- name
          self$gain <- gain
          self$inovSD <- inovSD
        },
        gain=0,
        inovSD=.1,
        draw = function(theta,deltaT,variants=list()){
          rnorm(length(theta),theta+self$gain*deltaT,
                self$inovSD*sqrt(deltaT))
        },
        lprob = function(par=self$pvec,theta0,theta1,weights,deltaT,
                         variants=list()) {
          mu <- theta0+par[1]*deltaT
          sig <- exp(par[2])*sqrt(deltaT)
          sum(dnorm(theta1,mu,sig,log=TRUE)*weights)
        },
        toString=function(digits=2,...){
          paste0("<BrownianGrowth: ", self$name, " ( ",
                 round(self$gain,digits=digits),
                 ", ",round(self$inovSD,digits=digits)," )>")
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


SpurtGrowth <- R6Class(
    classname="SpurtGrowth",
    inherit=GrowthModel,
    public=list(
        initialize = function(name,gain0,gain1,p,inovSD) {
          self$name <- name
          self$gain0 <- gain0
          self$gain1 <- gain1
          self$p <- p
          self$inovSD <- inovSD
        },
        gain0=0,
        gain1=1,
        p=.5,
        inovSD=.1,
        draw = function(theta,deltaT,variants=list()){
          mu <- ifelse(rbinom(length(theta),size=1, p=self$p),
                        self$gain1,self$gain0)*deltaT+
                theta
          sig <- self$inovSD*sqrt(deltaT)
          rnorm(length(theta),mu,sig)
        },
        lprob = function(par=self$pvec,theta0,theta1,weights,deltaT,
                         variants=list()) {
          mu0 <- theta0+par[1]*deltaT
          mu1 <- theta0+par[2]*deltaT
          p <- invlogit(par[3])
          sig <- exp(par[4])*sqrt(deltaT)
          sum(log((1-p)*dnorm(theta1,mu0,sig) +
                  p*dnorm(theta1,mu1,sig)
                  )*weights)
        },
        toString=function(digits=2,...){
          paste0("<SpurtGrowth: ", self$name, " ( ",
                 round(self$gain0,digits=digits), ",",
                 round(self$gain1,digits=digits), ",",
                 round(self$p,digits=digits), ",",
                 round(self$inovSD,digits=digits)," )>")
        }

    ),
    active=list(
        pvec = function(value) {
          if (missing(value)) return(c(self$gain0,self$gain1,
                                       logit(self$p),log(self$inovSD)))
          self$gain0 <- value[1]
          self$gain1 <- value[2]
          self$p <- invlogit(value[3])
          self$inovSD <- exp(value[4])
        }
    )
)



