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

gaussianWindow <- function(window=1,sigma=1,...,symmetric=TRUE) {
  if (symmetric)
    function(t,tstar) {ifelse(abs(tstar-t) > window,0,
                              exp(-1/2*((tstar-t)/sigma/window)^2))}
  else
    function(t,tstar) {ifelse((tstar-t) > window | t>tstar,0,
                              exp(-1/2*((tstar-t)/sigma/window)^2))}
}
