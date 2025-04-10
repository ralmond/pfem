test_that("BrownianGrowth initialize",{
  agr <- BrownianGrowth$new("B",.1,.04)
  expect_equal(agr$name,"B")
  expect_equal(agr$gain,.1)
  expect_equal(agr$inovSD,.04)
})

test_that("BrownianGrowth pvec", {
  agr <- BrownianGrowth$new("B",.1,.04)
  expect_equal(agr$pvec,c(.1,log(.04)))
  agr$pvec <- c(0,0)
  expect_equal(agr$gain,0)
  expect_equal(agr$inovSD,1)
})

test_that("BrownianGrowth draw", {
  tol <- .0001
  agr <- BrownianGrowth$new("B",0,tol/25)

  ## Theta
  samp <- agr$draw(c(-1,0,1),1)
  expect_equal(samp,c(-1,0,1),tolerance=tol)

  ## deltaT --mean
  agr$gain <- 1
  samp <- agr$draw(rep(0,3),1:3)
  expect_equal(samp,1:3,tolerance=tol)

  ## deltaT -- SD
  agr$gain <- 0
  agr$inovSD <- 1
  samp1 <- agr$draw(rep(0,100),1)
  expect_lt(sum(samp1^2),qchisq(.999,100))
  samp2 <- agr$draw(rep(0,100),4)
  expect_lt(sum(samp2^2),4*qchisq(.999,100))

})

test_that("BrownianGrowth lprob", {

  agr <- BrownianGrowth$new("B",1,1)
  theta0 <- rep(0,3)
  theta1 <- theta0+1:3
  weights <- 1:3
  expect_equal(agr$lprob(c(0,0),theta0,theta1,weights,1),
               sum(dnorm(theta1,log=TRUE)*weights))
  expect_equal(agr$lprob(c(1,log(2)),theta0,theta1,weights,1),
               sum(dnorm(theta1,1,2,log=TRUE)*weights))



})


test_that("SpurtGrowth initialize",{
  agr <- SpurtGrowth$new("S",0,.1,.25,.04)
  expect_equal(agr$name,"S")
  expect_equal(agr$gain0,0)
  expect_equal(agr$gain1,.1)
  expect_equal(agr$p,.25)
  expect_equal(agr$inovSD,.04)

})

test_that("SpurtGrowth pvec", {
  agr <- SpurtGrowth$new("S",0,1,.5,1)
  expect_equal(agr$pvec,c(0,1,0,0))
  agr$pvec<-c(-1,0,log(.25/.75),log(.04))
  expect_equal(agr$gain0,-1)
  expect_equal(agr$gain1,0)
  expect_equal(agr$p,.25)
  expect_equal(agr$inovSD,.04)
})

test_that("SpurtGrowth draw", {
  tol <- .0001
  agr <- SpurtGrowth$new("B",0,1,0,tol/25)

  ## Theta
  samp <- agr$draw(c(-1,0,1),1)
  expect_equal(samp,c(-1,0,1),tolerance=tol)

  ## p
  agr$p <- 1
  samp <- agr$draw(c(-1,0,1),1)
  expect_equal(samp,c(0,1,2),tolerance=tol)

  agr$p <- .5
  samp <- agr$draw(rep(0,100),1)
  z <- (sum(samp)-50)/sqrt(100*.5*.5)
  expect_lt(abs(z),qnorm(.999))

  ## deltaT --mean
  agr$gain0 <- 1
  agr$gain1 <- 0
  agr$p <- 0
  samp <- agr$draw(rep(0,3),1:3)
  expect_equal(samp,1:3,tolerance=tol)

  agr$gain0 <- 0
  agr$gain1 <- 1
  agr$p <- 1
  samp <- agr$draw(rep(0,3),1:3)
  expect_equal(samp,1:3,tolerance=tol)

  ## deltaT -- SD
  agr$gain0 <- 0
  agr$gain1 <- 0
  agr$inovSD <- 1
  samp1 <- agr$draw(rep(0,100),1)
  expect_lt(sum(samp1^2),qchisq(.999,100))
  samp2 <- agr$draw(rep(0,100),4)
  expect_lt(sum(samp2^2),4*qchisq(.999,100))

})

test_that("SpurtGrowth lprob", {
  agr <- SpurtGrowth$new("B",0,1,.5,1)
  theta0 <- rep(0,3)
  theta1 <- theta0+1:3
  weights <- 1:3
  expect_equal(agr$lprob(c(0,1,-Inf,0),theta0,theta1,weights,1),
               sum(dnorm(theta1,log=TRUE)*weights))
  expect_equal(agr$lprob(c(0,1,Inf,0),theta0,theta1,weights,1),
               sum(dnorm(theta1,1,1,log=TRUE)*weights))
  expect_equal(agr$lprob(c(0,1,0,0),theta0,theta1,weights,1),
               sum(log(dnorm(theta1)/2 +
                       dnorm(theta1,1,1)/2)*weights))

})

