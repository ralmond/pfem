test_that("NormalPop initialize",{
  apop <- NormalPop$new("Pop",3,1/2)
  expect_s3_class(apop,"PopulationModel")
  expect_equal(apop$name,"Pop")
  expect_equal(apop$mu,3)
  expect_equal(apop$sigma,1/2)

})

test_that("NormalPop pvec", {
  apop <- NormalPop$new("Pop",0,1)
  expect_equal(apop$pvec,c(0,0))
  apop$pvec <- c(1,log(2))
  expect_equal(apop$mu,1)
  expect_equal(apop$sigma,2)
})

test_that("NormalPop draw", {
  apop <- NormalPop$new("Pop",0,1)
  z <- apop$draw(100)
  expect_lt(sum(z^2),qchisq(.99,100))
})

test_that("NormalPop lprob", {
  apop <- NormalPop$new("Pop",0,1)
  theta <- 1:3
  weights <- 1:3
  expect_equal(apop$lprob(c(0,0),theta,weights),
               sum(dnorm(theta,log=TRUE)*weights))
  expect_equal(apop$lprob(c(1,log(2)),theta,weights),
               sum(dnorm(theta,1,2,log=TRUE)*weights))

})


test_that("NormalPop cdf",{
  apop <- NormalPop$new("standard",0,1)
  expect_equal(apop$cdf(c(-Inf,0,Inf)),c(0,.5,1))
})

test_that("NormalPop mstep",{

})
