library(testthat)

source("script1.R");
## Test norm data
test_that("norm_data detects nonpositive values", {
  expect_error(norm_data(c(1, -1)), "All values should be positive")
  expect_error(norm_data(c(1, 0)), "All values should be positive")
})

test_that("normdata returns matrices", {
  expect_equal(norm_data(c(1, 1)), matrix(c(1, 1), nrow = 1))
  expect_equal(norm_data(matrix(c(1, 1, 1, 4), nrow = 2)),
               matrix(c(1, 1, 1, 4), nrow = 2))
})

## Test closure

test_that("closure works properly", {
  expect_equal(closure(c(1, 1)), matrix(c(0.5, 0.5), nrow = 1))
  expect_equal(closure(matrix(c(1, 1, 1, 4), nrow = 2)),
               matrix(c(0.5, 0.2, 0.5, 0.8), nrow = 2))
  expect_equal(closure(c(1, 1), k = 4), matrix(c(2, 2), nrow = 1))
})



## Test Inner product

test_that("Inner product works properly", {
  expect_equal(Inner_product(c(1.2, 1.9),c(8.1, 9.2)),0.02926,tolerance=1e-5)
  expect_equal(Inner_product(c(7.1, 6.2),c(0.3, 7.5)),-0.21815,tolerance=1e-5)
  expect_equal(Inner_product(c(12.8, 3.2),c(18.3, 73.5)),-0.96374,tolerance=1e-5)
  
  expect_error(Inner_product(c(1.2, 1.9),c(8.1, 9.2,2.3)), "x and y should have the same dimension")
  expect_error(Inner_product(matrix(c(8.1, 9.2, 9.2, 8.1),nrow=2),matrix(c(8.1, 9.2, 9.2, 8.1),nrow=2)), "wrong dimension")
  
})


## Test pertubation

test_that("closure works properly", {
  expect_equal(perturbation(matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2), c(2, 3)),matrix(c(0.4, 0.727, 0.6, 0.272), nrow=2), tolerance=1e-2)
  expect_equal(perturbation(matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2), 3 ), matrix(c(0.5, 0.8, 0.5, 0.2), nrow = 2), tolerance=1e-2)
  
})
