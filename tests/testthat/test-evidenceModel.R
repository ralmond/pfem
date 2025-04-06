test_that("GradedResponse initialize",{
  em <- GradedResponse$new("gr",1,c(-1,1))
  expect_s3_class(em,"EvidenceModel")
  expect_equal(em$name,"gr")
  expect_equal(em$a,1)
  expect_equal(em$b,c(-1,1))

})

test_that("GradedResponse pvec", {
  em <- GradedResponse$new("gr",1,c(-1,1))
  expect_equal(em$pvec,c(0,-1,1))
  em$pvec <- c(log(2),0,1)
  expect_equal(em$a,2)
  expect_equal(em$b,0:1)
})

test_that("GradedResponse cuts", {
  em <- GradedResponse$new("gr",1,c(-1,1))
  theta <- -2:2
  expect_equal(logit(em$cuts(theta)),
               cbind(theta+1,theta-1))
  expect_equal(logit(em$cuts(theta,1/2,c(0,1))),
               cbind(theta/2,theta/2-1))

})


test_that("GradedResponse draw", {
  em <- GradedResponse$new("gr",1,c(-1,1))
  theta <- -2:2
  Y <- em$draw(theta)
  expect_length(Y,length(theta))
  expect_lte(max(Y),2)
  expect_gte(min(Y),0)
})

test_that("GradedResponse lprob", {
  em <- GradedResponse$new("gr",1,c(-1,1))
  theta <- -2:2
  weights <- c(1,1,2,3,3)
  Y <- c(0,0,1,2,2)
  probs <- cuts2probs(em$cuts(theta))
  expect_equal(em$lprob(c(0,-1,1),Y,
                        theta,weights),
               sum(weights*
                   log(c(probs[1,1],
                         probs[2,1],
                         probs[3,2],
                         probs[4,3],
                         probs[5,3]))))
  probs <- cuts2probs(em$cuts(theta,.5,c(-1,1)))
  expect_equal(em$lprob(c(log(.5),-1,1),Y,
                        theta,weights),
               sum(weights*
                   log(c(probs[1,1],
                         probs[2,1],
                         probs[3,2],
                         probs[4,3],
                         probs[5,3]))))

})

test_that("GradedResponse llike", {
  em <- GradedResponse$new("gr",1,c(-1,1))
  theta <- -2:2
  cuts <- em$cuts(theta)
  expect_equal(em$llike(0,theta),
               log(1-cuts[,1]))
  expect_equal(em$llike(1,theta),
               log(cuts[,1]-cuts[,2]))
  expect_equal(em$llike(2,theta),
               log(cuts[,2]))
})



test_that("NormalScore initialize",{
  em <- NormalScore$new("ns",0,1)
  expect_s3_class(em,"EvidenceModel")
  expect_equal(em$name,"ns")
  expect_equal(em$bias,0)
  expect_equal(em$se,1)

})

test_that("NormalScore pvec", {
  em <- NormalScore$new("ns",0,1)
  expect_equal(em$pvec,c(0,0))
  em$pvec <- c(1,log(2))
  expect_equal(em$bias,1)
  expect_equal(em$se,2)
})

test_that("NormalScore draw", {
  tol <- .000001
  em <- NormalScore$new("ns",1,tol/25)
  theta <- -2:2
  Y <- em$draw(theta)
  expect_equal(Y,theta+1,tolerance=tol)
})

test_that("NormalScore lprob", {
  em <- NormalScore$new("ns",0,1)
  theta <- -2:2
  weights <- c(1,1,2,3,3)
  expect_equal(em$lprob(c(0,0),0, theta,weights),
               sum(weights*dnorm(theta,log=TRUE)))
  expect_equal(em$lprob(c(1,log(.5)),0,theta,weights),
               sum(weights*(dnorm(0,theta+1,.5,log=TRUE))))
})

test_that("NormalScore llike", {
  em <- NormalScore$new("ns",0,1)
  theta <- -2:2
  expect_equal(em$llike(0,theta),dnorm(theta,log=TRUE))
})




