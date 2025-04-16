EM3 <- list(NormalScore$new("EM-1",-1,.001),
            NormalScore$new("EM 0",0,.001),
            NormalScore$new("EM 1",1,.001))
task3 <- matrix(c(c(2,2,2),3:1),2,3,byrow=TRUE)
data3 <- matrix(c(-1:1,1:-1),2,3,byrow=TRUE)


test_that("IRTF",{
  Afilt<-IRTF$new("testf",3,list(),EM3)
  Afilt$qpoints <- c(-1,0,1)
  Afilt$tasks <- task3
  Afilt$data <- data3
  Afilt$times <- matrix(0:3,1,4)

  expect_equal(Afilt$nsubjects,2)
  expect_equal(Afilt$maxocc,3)
  expect_equal(Afilt$time(2,2),2)
  expect_equal(Afilt$group(2),1)
  expect_equal(Afilt$task(1,2),2)
  expect_equal(Afilt$task(2,1),3)


})

test_that("irtf",{
  Afilt<-IRTF$new("testf",3,list(),EM3)
  Afilt$qpoints <- c(-1,0,1)
  Afilt$tasks <- task3
  Afilt$data <- data3
  Afilt$times <- matrix(0:3,1,4)

  irtf(Afilt,1:3,linearWindow(.5),debug=TRUE)
  expect_equal(dim(Afilt$lweights),c(3,3,2))
  expect_equal(Afilt$lweights[,1,1],rev(Afilt$lweights[,3,1]),
               tolerance=.001)
  expect_equal(Afilt$lweights[,1,2],Afilt$lweights[,3,2],
               tolerance=.001)
  expect_equal(Afilt$lweights[,2,1],Afilt$lweights[,2,2],
               tolerance=.001)
  expect_equal(Afilt$lweights[,3,2],Afilt$lweights[,2,2],
               tolerance=.001)
  expect_equal(Afilt$weights[,,1],diag(3))
  expect_equal(Afilt$weights[,,2],matrix(c(0,1,0),3,3))

})
