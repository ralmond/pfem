squareWindow <- function(window=1,...,symmetric=TRUE) {
  if (symmetric)
    function(t,tstar) {as.numeric(abs(tstar-t) <= window)}
  else
    function(t,tstar) {as.numeric((tstar-t) <= window & tstar >= t)}
}

linearWindow <- function(window=1,slope=1,...,symmetric=TRUE) {
  if (symmetric)
    function(t,tstar) {ifelse(abs(tstar-t) > window,0,1-slope*abs(tstar-t)/window)}
  else
    function(t,tstar) {ifelse((tstar-t) > window | t>tstar,0,
                              1-slope*(tstar-t)/window)}
}

exponentialWindow <- function(window=1,lambda=.8,...,symmetric=TRUE) {
  if (symmetric)
    function(t,tstar) {lambda^abs(tstar-t)}
  else
    function(t,tstar) {ifelse(t>tstar,0,lambda^(tstar-t))}
}

gaussianWindow <- function(window=1,sigma=.5,power=1L,...,symmetric=TRUE) {
  power <- floor(power)
  if (symmetric)
    function(t,tstar) {ifelse(abs(tstar-t) > window,0,
                              exp(-1/2*((tstar-t)/sigma/window)^(2*power)))}
  else
    function(t,tstar) {ifelse((tstar-t) > window | t>tstar,0,
                              exp(-1/2*((tstar-t)/sigma/window)^(2*power)))}
}

parzenWindow <- function(window=1,...,symmetric=TRUE) {
  if (symmetric)
    function (t,tstar) { dt <- abs(tstar-t)/window
      ifelse(dt>1,0,ifelse(dt < 1/2, 1-6*dt^2*(1-dt),
                           2*(1-dt)^3))
    }
  else
    function (t,tstar) { dt <- (tstar-t)/window
      ifelse(dt>1|dt<0,0,ifelse(dt < 1/2, 1-6*dt^2*(1-dt),
                           2*(1-dt)^3))
    }
}

welchWindow <- function(window=1,slope=1,...,symmetric=TRUE) {
  if (symmetric)
    function(t,tstar) {ifelse(abs(tstar-t) > window,0,
                              1-slope*(abs(tstar-t)/window)^2)}
  else
    function(t,tstar) {ifelse((tstar-t) > window | t>tstar,0,
                              1-slope*((tstar-t)/window)^2)}
}

sineWindow <- function(window=1,power=1,...,symmetric=TRUE) {
  if (symmetric)
    function(t,tstar) {ifelse(abs(tstar-t) > window,0,
                              cos(pi*(tstar-t)/window/2)^power)}
  else
    function(t,tstar) {ifelse((tstar-t) > window | t>tstar,0,
                              cos(pi*(tstar-t)/window/2)^power)}
}

HannConstants <- c(.5,.5)
HammingConstants <- c(25,21)/46
BlackmanConsts <- function (alpha=.16)
  c((1-alpha)/2,1/2,alpha/2)
BlackmanExact <- c(7938,9240,1430)/18608
NuttallConstants <-c(0.355768, 0.487396, 0.144232, 0.012604)
BlackmanNuttallConstants <- c(0.3635819, 0.4891775, 0.1365995, 0.0106411)
BlackmanHarrisConstants <-c(0.35875, 0.48829, 0.14128, 0.01168)

FlattopConstants <- c(a0=0.21557895, a1=0.41663158, a2=0.277263158,
                      a3=0.083578947, a4=0.006947368)


cosineSumWindow <- function(window=1,
                            a=list(HannConstants,HammingConstants,
                                   BlackmanConsts(.16),BlackmanExact,
                                   NuttallConstants,BlackmanNuttallConstants,
                                   BlackmanHarrisConstants,FlattopConstants),
                            ...,symmetric=TRUE) {
  if (is.list(a)) a <- a[[1L]]
  if (symmetric)
    function(t,tstar) {
      dt <- (tstar-t)/window
      result <- a[1L]
      for (k in 2L:length(a)) {
        result <- result + (-1)^(k-1L)*a[k]*cos(pi*(k-1)*(1-dt))
      }
      ifelse(abs(dt)>1,0,result)
    } else
      function(t,tstar) {
      dt <- (tstar-t)/window
      result <- a[1L]
      for (k in 2L:length(a)) {
        result <- result + (-1)^(k-1L)*a[k]*cos(pi*(k-1)*(1-dt))
      }
      ifelse(dt>1|dt<0,0,result)
    }
}

tukeyWindow <- function(window=1,alpha=1/2,...,symmetric=TRUE) {
  if (symmetric)
    function(t,tstar) {
      dt <- abs(tstar-t)/window
      ifelse(dt > 1,0,
             ifelse(1-alpha/2 > dt, 1, 1/2*(1+cos(2*pi*(dt+1-alpha/2)/alpha))))
      }
    else
      function(t,tstar) {
        dt <- (tstar-t)/window
      ifelse(dt > 1 | dt < 0,0,
             ifelse(1-alpha/2 > dt, 1, 1/2*(1+cos(2*pi*(dt+1-alpha/2)/alpha))))
      }
}

plankWindow <- function(window=1,epsilon=1/4,...,symmetric=TRUE) {
  if (symmetric)
    function (t, tstar) {
      dt <- abs(tstar-t)
      nn <- window - dt
      ifelse(dt > window, 0,
        ifelse(dt <= (1-2*epsilon)*window,1,
               1/(1+exp(2*epsilon*window/nn -
                        2*epsilon*window/(2*epsilon*window - nn)
               ))))
    }
  else
    function (t, tstar) {
      dt <- tstar-t
      nn <- window - dt
      ifelse(dt > window | dt < 0, 0,
        ifelse(dt <= (1-2*epsilon)*window,1,
               1/(1+exp(2*epsilon*window/nn -
                        2*epsilon*window/(2*epsilon*window - nn)
               ))))
    }
}

