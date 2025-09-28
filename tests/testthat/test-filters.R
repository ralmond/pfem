test_that("square window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- squareWindow()
  expect_equal(w(t,tstar)>0,abs(t)<=1)
  w2 <- squareWindow(window=2)
  expect_equal(w2(2*t,tstar)>0,abs(t)<=1)
  wleft <- squareWindow(symmetric=FALSE)
  expect_equal(wleft(t,tstar)>0,-1<=t & t<=0)

})


test_that("linear window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- linearWindow()
  expect_equal(w(t,tstar)>0,abs(t)<1)
  w2 <- linearWindow(window=2)
  expect_equal(w2(2*t,tstar)>0,abs(t)<1)
  wleft <- linearWindow(symmetric=FALSE)
  expect_equal(wleft(t,tstar)>0,-1<t & t<=0)
  ws2 <- linearWindow(slope=1/2)
  expect_equal(ws2(t,tstar)>0,abs(t)<=1)

})


test_that("exponential window", {
  tstar <- 0
  t <- c(-2,-1,0,1,2)
  w5 <- exponentialWindow(lambda=.5)
  expect_equal(w5(t,tstar),c(.25,.5,1,.5,.25))
  w5left <- exponentialWindow(lambda=.5,symmetric=FALSE)
  expect_equal(w5left(t,tstar),c(.25,.5,1,0,0))
  w2 <- exponentialWindow(lambda=.2)
  expect_equal(w2(t,tstar),c(.04,.2,1,.2,.04))

})

test_that("Gaussian window", {
  tstar <- 0
  t <- c(-1.1,-1,0,1,1.1)
  w <- gaussianWindow(sigma=.5)
  wmin <- sqrt(2*pi)*.5*dnorm(-1,sd=.5)
  expect_equal(w(t,tstar),c(0,wmin,1,wmin,0),tolerance=.0001)
  w2 <- gaussianWindow(window=2,sigma=.5)
  expect_equal(w2(2*t,tstar),c(0,wmin,1,wmin,0),tolerance=.0001)
  wleft <- gaussianWindow(sigma=.5,symmetric=FALSE)
  expect_equal(wleft(t,tstar),c(0,wmin,1,0,0),tolerance=.0001)
  ws4 <- gaussianWindow(sigma=.4,symmetric=FALSE)
  wmin4 <-sqrt(2*pi)*.4*dnorm(-1,sd=.4)
  expect_equal(ws4(t,tstar),c(0,wmin4,1,0,0),tolerance=.0001)

})


test_that("Gaussian 2 window", {
  tstar <- 0
  t <- c(-1.1,-1,0,1,1.1)
  w <- gaussianWindow(power=2,sigma=.5)
  wmin <- sqrt(2*pi)*dnorm((1/.5)^2)
  expect_equal(w(t,tstar),c(0,wmin,1,wmin,0),tolerance=.0001)
  w2 <- gaussianWindow(window=2,power=2,sigma=.5)
  expect_equal(w2(2*t,tstar),c(0,wmin,1,wmin,0),tolerance=.0001)
  wleft <- gaussianWindow(power=2,sigma=.5,symmetric=FALSE)
  expect_equal(wleft(t,tstar),c(0,wmin,1,0,0),tolerance=.0001)
  ws4 <- gaussianWindow(power=2,sigma=.4,symmetric=FALSE)
  wmin4 <-sqrt(2*pi)*dnorm((1/.4)^2)
  expect_equal(ws4(t,tstar),c(0,wmin4,1,0,0),tolerance=.0001)

})

test_that("Parzen window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- parzenWindow()
  expect_equal(w(t,tstar)>0,abs(t)<1)
  w2 <- parzenWindow(window=2)
  expect_equal(w2(2*t,tstar)>0,abs(t)<1)
  wleft <- parzenWindow(symmetric=FALSE)
  expect_equal(wleft(t,tstar)>0,-1<t & t<=0)

})


test_that("Welch window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- welchWindow()
  expect_equal(w(t,tstar)>0,abs(t)<1)
  w2 <- welchWindow(window=2)
  expect_equal(w2(2*t,tstar)>0,abs(t)<1)
  wleft <- welchWindow(symmetric=FALSE)
  expect_equal(wleft(t,tstar)>0,-1<t & t<=0)

})

test_that("sine window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- sineWindow()
  expect_equal(w(t,tstar)>0.0001,abs(t)<1)
  w2 <- sineWindow(window=2)
  expect_equal(w2(2*t,tstar)>0.0001,abs(t)<1)
  wleft <- sineWindow(symmetric=FALSE)
  expect_equal(wleft(t,tstar)>0.0001,-1<t & t<=0)

})


test_that("sine 2 window", {
  t <- seq(-1.1,1.1,.1)
  tstar <- 0
  sw <- sineWindow(power=2)
  hw <- cosineSumWindow(a=HannConstants)
  expect_equal(sw(t,tstar),hw(t,tstar),tolerance=.0001)

})


test_that("cosine sum window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- cosineSumWindow(a=HannConstants)
  expect_equal(w(t,tstar)>0.0001,abs(t)<1)
  w2 <- cosineSumWindow(window=2,a=HannConstants)
  expect_equal(w2(2*t,tstar)>0.0001,abs(t)<1)
  wleft <- cosineSumWindow(a=HannConstants,symmetric=FALSE)
  expect_equal(wleft(t,tstar)>0.0001,-1<t & t<=0)

  expect_equal(cosineSumWindow(a=BlackmanConsts(.16))(t,tstar)>0.0001,
               abs(t)<1)

    for (cons in c("HannConstants",
                   "NuttallConstants",
                   "BlackmanHarrisConstants")) {
      w <- cosineSumWindow(a=get(cons))
      expect_equal(w(t,tstar)>0.0001,abs(t)<1,label=cons)
    }
    for (cons in c("HammingConstants",
                   "BlackmanExact",
                   "BlackmanNuttallConstants")) {
      w <- cosineSumWindow(a=get(cons))
      expect_equal(w(t,tstar)>0.0001,abs(t)<=1,label=cons)
    }
    w <- cosineSumWindow(a=FlattopConstants)
    expect_equal(w(c(-1.1,1.1),0),c(0,0))
    tt <- seq(-1.0,1.0,.1)
    expect_equal(w(tt,0)<0,abs(tt)>=.5)
    expect_equal(w(tt,0)>0,abs(tt)<.5)
    expect_equal(w(0,0),1,tolerance=.0001)
})


test_that("Tukey window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- tukeyWindow()
  expect_equal(w(t,tstar)==1,abs(t)<=.75)
  w2 <- tukeyWindow(window=2)
  expect_equal(w2(2*t,tstar)==1,abs(t)<=.75)
  wleft <- tukeyWindow(symmetric=FALSE)
  expect_equal(wleft(t,tstar)==1,t>=-.75 & t<=0)
  w4 <- tukeyWindow(alpha=.4)
  expect_equal(w4(t,tstar)==1,abs(t)<=.8)

})

test_that("Plank window", {
  tstar <- 0
  t <- seq(-1.1,1.1,.1)
  w <- plankWindow()
  expect_equal(w(t,tstar)==1,abs(t)<=1/2)
  w2 <- plankWindow(window=2)
  expect_equal(w2(2*t,tstar)==1,abs(t)<=1/2)
  wleft <- plankWindow(symmetric=FALSE)
  expect_equal(wleft(t,tstar)==1,t>=-.5 & t<=0)
  w4 <- plankWindow(epsilon=.2)
  expect_equal(w4(t,tstar)==1,abs(t)<=.61)

})

