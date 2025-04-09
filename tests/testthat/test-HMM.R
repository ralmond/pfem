make_pops <- function() {
  list(NormalPop$new("P1",0,1),
       NormalPop$new("P2",1,1))
}

make_growths <- function() {
  list(BrownianGrowth$new("G1",1,1),
       BrownianGrowth$new("G2",2,1))
}

make_ems <- function() {
  list(NormalScore$new("I1",0,.5),
       NormalScore$new("I2",.25,.25))
}


test_that("HMM new", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  expect_equal(hmm$name,"test")
  expect_equal(sapply(hmm$popModels,name.R6),
               c("P1","P2"))
  expect_equal(sapply(hmm$growthModels,name.R6),
               c("G1","G2"))
  expect_equal(sapply(hmm$evidenceModels,name.R6),
               c("I1","I2"))

})

test_that("HMM delT", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  hmm$data <- matrix(0,2,3)
  hmm$deltaT <- matrix(1:6,2,3)
  expect_equal(hmm$delT(1,1),1)
  expect_equal(hmm$delT(2,2),4)

  hmm$deltaT <- matrix(1:6,1,6)
  expect_equal(hmm$delT(3,3),3)

  hmm$deltaT <- matrix(1:6,6,1)
  expect_equal(hmm$delT(2,2),2)

})

test_that("HMM group", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  hmm$data <- matrix(0,2,3)
  hmm$deltaT <- matrix(1:6,2,3)
  hmm$groups <- 1:2
  expect_equal(hmm$group(1),1)
  expect_equal(hmm$group(2),2)
  hmm$groups <- 1
  expect_equal(hmm$group(2),1)

})

test_that("HMM action", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  hmm$data <- matrix(0,2,3)
  hmm$actions <- matrix(1:6,2,3)
  expect_equal(hmm$action(1,1),1)
  expect_equal(hmm$action(2,2),4)

  hmm$actions <- matrix(1:6,1,6)
  expect_equal(hmm$action(3,2),2)

  hmm$actions <- matrix(1:6,6,1)
  expect_equal(hmm$action(4,2),4)


})

test_that("HMM task", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  hmm$data <- matrix(0,2,3)
  hmm$tasks <- matrix(1:6,2,3)
  expect_equal(hmm$task(1,1),1)
  expect_equal(hmm$task(2,2),4)

  hmm$tasks <- matrix(1:6,1,6)
  expect_equal(hmm$task(3,2),2)

  hmm$tasks <- matrix(1:6,6,1)
  expect_equal(hmm$task(4,2),4)

})

test_that("HMM deltaT", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  hmm$data <- matrix(0,2,3)
  hmm$deltaT <- rep(1,5)
  expect_equal(hmm$deltaT,matrix(1,1,5))
  expect_equal(hmm$times,matrix(0:5,1,6))

  hmm$deltaT <- matrix(1:2,2,3)
  expect_equal(hmm$delT(1,3),1)
  expect_equal(hmm$delT(2,3),2)
  expect_equal(hmm$times,matrix(c(0:3,2*(0:3)),2,4,byrow=TRUE))

})

test_that("HMM times", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  hmm$data <- matrix(0,2,3)
  hmm$times <- 0:5
  expect_equal(hmm$deltaT,matrix(1,1,5))
  expect_equal(hmm$times,matrix(0:5,1,6))

  hmm$times <- matrix(c(0:3,2*(0:3)),2,4,byrow=TRUE)
  expect_equal(hmm$delT(1,3),1)
  expect_equal(hmm$delT(2,3),2)
  expect_equal(hmm$deltaT,matrix(1:2,2,3))

})

test_that("HMM dims", {
  hmm <- HMM$new("test",make_pops(),
                 make_growths(),
                 make_ems())
  hmm$data <- matrix(0,2,3)
  expect_equal(hmm$nsubjects,2)
  expect_equal(hmm$maxocc,3)
})


test_that("HMM particleFilter", {
  tol <- .0001
  hmm<-HMM$new("stimulator",
               list(NormalPop$new("P1",-1/2,1),
                    NormalPop$new("P2",1/2,1)),
               list(BrownianGrowth$new("G1",1,tol/100)),
               list(NormalScore$new("N0",0,1)))
  hmm$clspec <- 1
  hmm$groups <- 1:2
  hmm$deltaT<-rep(1,3)
  hmm$data <- outer(c(-1,1)/2,1:3,"+")
  particleFilter(hmm,5)
  expect_equal(dim(hmm$theta),c(5,2,4))
  expect_equal(dim(hmm$weights),c(5,2))
  expect_equal(colSums(hmm$weights),c(1,1),
               tolerance=tol)
  expect_equal(apply(hmm$theta[,1L,],1,diff),
               matrix(1,3,5),tolerance=tol)
  expect_equal(apply(hmm$theta[,2L,],1,diff),
               matrix(1,3,5),tolerance=tol)
  llikes <-dnorm(hmm$theta[,1L,1L],-1/2,1,log=TRUE)
  expect_equal(hmm$lweights[,1L]/llikes,rep(3,5),
               tolerance=5*tol)
})


test_that("HMM longResults simulation", {
  tol <- .0001
  hmm<-HMM$new("stimulator",
               list(NormalPop$new("P0",0,1)),
               list(BrownianGrowth$new("G1",1,tol/100)),
               list(NormalScore$new("N0",0,tol/100)))
  hmm$deltaT<-rep(1,3)
  simulate(hmm,5,mocc=3)
  table <- longResults(hmm)
  expect_equal(dim(table),c(20,7))
  expect_equal(table$subj[1:5],1:5)
  expect_equal(table$occ,table$time)
  expect_equal(is.na(table$Y),c(rep(TRUE,5),rep(FALSE,15)))
  expect_equal(is.na(table$tasks),c(rep(TRUE,5),rep(FALSE,15)))
  expect_true(all(is.na(table$weights)))
  expect_true(all(!is.na(table$theta)))
})

test_that("HMM longResults PF", {
  tol <- .0001
  hmm<-HMM$new("stimulator",
               list(NormalPop$new("P1",-1/2,1),
                    NormalPop$new("P2",1/2,1)),
               list(BrownianGrowth$new("G1",1,tol/100)),
               list(NormalScore$new("N0",0,1)))
  hmm$clspec <- 1
  hmm$groups <- 1:2
  hmm$deltaT<-rep(1,3)
  hmm$data <- outer(c(-1,1)/2,1:3,"+")
  particleFilter(hmm,5)
  table <- longResults(hmm)
  expect_equal(dim(table),c(40,7))
  expect_equal(table$subj,rep(c(rep(1,5),rep(2,5)),4))
  expect_equal(table$occ,table$time)
  expect_equal(is.na(table$Y),c(rep(TRUE,10),rep(FALSE,30)))
  expect_equal(is.na(table$tasks),c(rep(TRUE,10),rep(FALSE,30)))
  expect_equal(sum(table$weights),8,tolerance=tol)
  expect_true(all(!is.na(table$theta)))
})



test_that("HMM simulate", {
  tol <- .0001
  hmm<-HMM$new("stimulator",
          list(NormalPop$new("P0",0,1)),
          list(BrownianGrowth$new("G1",1,tol/100)),
          list(NormalScore$new("N0",0,tol/100)))
  hmm$deltaT<-rep(1,3)
  simulate(hmm,5,mocc=3)
  expect_equal(dim(hmm$theta),c(1,5,4))
  expect_equal(dim(hmm$data),c(5,3))
  expect_equal(apply(hmm$theta[1,,],1,diff),
               matrix(1,3,5),tolerance=tol)
  expect_equal(hmm$theta[1,,-1],hmm$data,tolerance=tol)
})

